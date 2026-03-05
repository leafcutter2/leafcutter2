#!/usr/bin/env python


# TODO: refactor chromLst out of module-level scope — pass as argument to
# pool_junc_reads(), addlowusage(), and sort_junctions() instead of relying on module global.
chromLst = (
    [f"chr{x}" for x in range(1, 23)]
    + ["chrX", "chrY"]
    + [f"{x}" for x in range(1, 23)]
    + ["X", "Y"]
)

import argparse
import gzip
import os
import re
import shutil
import sys
import tempfile
from statistics import stdev
from datetime import datetime

from leafcutter2 import classifier as sjcf
import pandas as pd
import pyfastx
import logging
from leafcutter2 import transcript_tools as Transcript_tools
import shlex

logger = logging.getLogger(__name__)

def setup_logging(level_name: str, verbose_flag: bool):
    # Resolve level precedence: explicit level overrides verbose flag
    level = getattr(logging, level_name.upper(), None)
    if level is None:
        level = logging.DEBUG if verbose_flag else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')
    logger.debug(f"Logger initialized at level: {logging.getLevelName(level)}")


def natural_sort(l: list):
    """Natural sort a list of string/tuple, similar to bash `sort -V`

    Parameters:
    -----------
    l : list
        l can be a list of string or numerics; or a list of varing length of tuples

    Returns:
    --------
    return : a sorted list
    """

    untuple = lambda tup: "".join([str(e) for e in tup])
    convert = lambda text: int(text) if text.isdigit() else text.lower()
    alphanum_key = lambda key: [convert(c) for c in re.split("([0-9]+)", untuple(key))]
    return sorted(l, key=alphanum_key)


def cluster_intervals(E: list):
    """Clusters intervals together

    Parameters:
    -----------
    E : list
        list of tuples, e.g. [(start, end)]

    Returns: tuple
    ---------
    Eclusters : list of list
        Each element is a list of tuples (introns) clustered into 1 cluster
    cluster: list
        A list of tuples of introns
    """
    E = natural_sort(E)
    current = E[0]
    Eclusters, cluster = [], []

    i = 0
    while i < len(E):

        if overlaps(E[i], current):
            cluster.append(E[i])
        else:
            Eclusters.append(cluster)
            cluster = [E[i]]
        current = (E[i][0], max([current[1], E[i][1]]))
        i += 1

    if len(cluster) > 0:

        Eclusters.append(cluster)

    return Eclusters, E


def overlaps(A: tuple, B: tuple):
    """Checks if A and B overlaps

    Parameters:
    -----------
    A : tuple
        start and end coordinates of 1 intron
    B : tuple
        start and end coordinates of another intron

    Returns:
    --------
    return : boolean
        Indicates whether genomic ranges of A and B overlap or not.
    """

    if A[1] < B[0] or B[1] < A[0]:
        return False
    else:
        return True


def pool_junc_reads(flist, options):
    """Pool junction reads

    Parameters
    ----------
    flist : str
        The file list
    options : argparse object
        Passed arguments

    Returns:
    --------
    return
        No returns. Use side effects.

    Side-effects:
    -------------
    write introns and counts by clusters. Output file is NOT versions sorted.
    format: [chrom]:[strand] [start]:[end]:[reads]
            e.g. chr17:+ 410646:413144:3 410646:413147:62
    """


    outPrefix = options.outprefix
    rundir = options.rundir
    maxIntronLen = int(options.maxintronlen)
    checkchrom = options.checkchrom
    logger.info(f"Max Intron Length: {maxIntronLen}")
    outFile = f"{rundir}/clustering/{outPrefix}_pooled"

    if not os.path.exists(rundir):
        os.mkdir(rundir)

    # store introns in `by_chrom`, a nested dictionary
    by_chrom = {}  # { k=(chrom, strand) : v={ k=(start, end) : v=reads } }

    for libl in flist:
        # FIX: Extract only the filepath part, ignoring sample name
        lib = libl.strip().split('\t')[0].strip()
        
        if not os.path.isfile(lib):
            continue

        if options.verbose:
            logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} scanning {lib}...\n")

        if ".gz" in lib:
            F = gzip.open(lib)
        else:
            F = open(lib)

        for ln in F:

            if type(ln) == bytes:
                ln = ln.decode("utf-8")  # convert bytes to string

            lnsplit = ln.split()
            if len(lnsplit) < 6:
                logger.error(f"Error in {lib} \n")
                continue

            if len(lnsplit) == 12:  # 12 fields regtools junc file
                (
                    chrom,
                    A,
                    B,
                    dot,
                    counts,
                    strand,
                    rA,
                    rb,
                    rgb,
                    blockCount,
                    blockSize,
                    blockStarts,
                ) = lnsplit
                if int(blockCount) > 2:
                    logger.debug("ignored junction with blockCount>2")
                    continue
                Aoff, Boff = blockSize.split(",")[:2]
                A, B = int(A) + int(Aoff), int(B) - int(Boff)  # get intron

            elif len(lnsplit) == 6:
                # old leafcutter junctions
                chrom, A, B, dot, counts, strand = lnsplit
                A, B = int(A), int(B)

            if checkchrom and (chrom not in chromLst):
                continue

            A, B = int(A), int(B) + int(options.offset)

            if B - A > int(maxIntronLen):
                continue

            # sum up all the reads at the same junctions if junctions already exist
            try:
                by_chrom[(chrom, strand)][(A, B)] = (
                    int(counts) + by_chrom[(chrom, strand)][(A, B)]
                )
            except:
                try:
                    by_chrom[(chrom, strand)][(A, B)] = int(
                        counts
                    )  # when only 1 junction
                except:
                    by_chrom[(chrom, strand)] = {
                        (A, B): int(counts)
                    }  # when only 1 junction

    with open(outFile, "w") as fout:
        Ncluster = 0
        logger.info("Parsing pooled junctions...")

        for chrom in by_chrom:
            read_ks = [
                k for k, v in by_chrom[chrom].items() if v >= 3
            ]  # read-keys, require junction reads > 3
            read_ks.sort()  # sort read-keys: (start, end)

            logger.info(f"{chrom[0]}:{chrom[1]}..\n")

            if read_ks:
                clu = cluster_intervals(read_ks)[
                    0
                ]  # clusters of introns, [[(start, end),..],..]

                for cl in clu:
                    if len(cl) > 0:  # 1 if cluster has more than one intron
                        buf = f"{chrom[0]}:{chrom[1]} "  # chr:strand
                        for interval, count in [(x, by_chrom[chrom][x]) for x in cl]:
                            buf += (
                                f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                            )  # start:end:reads
                        fout.write(buf + "\n")
                        Ncluster += 1

        logger.info(f"Wrote {Ncluster} clusters..\n")


def refine_linked(clusters):
    """Re-cluster introns into clusters of linked introns

    Linked introns are introns that share either 5' or 3' splice site

    Parameters:
    -----------
    clusters : tuple
        format is [((start, end), reads)], eg.
        [((413430, 423479), 3), ((410646, 413144), 3), ((410646, 413147), 62), ((410646, 420671), 4)]

    Returns:
    --------
    return : list of list
        base element is a tuple, format: ((start, end), reads). e.g.
        [[((413430, 423479), 3)], [((410646, 413144), 3), ((410646, 413147), 62), ((410646, 420671), 4)]]
    """

    unassigned = [x for x in clusters[1:]]
    current = [clusters[0]]
    splicesites = set([current[0][0][0], current[0][0][1]])  # start & end of intron
    newClusters = []
    while len(unassigned) > 0:
        finished = False

        while not finished:
            finished = True
            torm = []
            for intron in unassigned:
                (start, end), count = intron
                if start in splicesites or end in splicesites:
                    current.append(intron)
                    splicesites.add(start)
                    splicesites.add(end)
                    finished = False
                    torm.append(intron)
            for intron in torm:
                unassigned.remove(intron)
        newClusters.append(current)
        current = []
        if len(unassigned) > 0:
            current = [unassigned[0]]
            splicesites = set([current[0][0][0], current[0][0][1]])
            if len(unassigned) == 1:  # assign the last intron
                newClusters.append(current)
            unassigned = unassigned[1:]
    return newClusters


def refine_cluster(clu: list, cutoff: float, readcutoff: int):
    """Filter introns based on cutoffs

    Parameters:
    -----------
    clu : list of tuples, a single cluster
        list of tuples, each tuple an intron of the cluster
    cutoff : float
        reads ratio cutoff, passed in from option --mincluratio
    readcutoff : int
        minimum reads cutoff, passed in from option --minreads

    Filters:
    --------
        1. compute ratio of reads for each intron in a cluster
        2. remove intron if:
            - ratio < ratio_cutoff
            - OR reads of the intron < readcutoff
        3. re-cluster with remaining introns

    Returns:
    --------
        return : list of list
        list of refined (filtered) clusters

    """

    remove = []
    dic = {}
    intervals = []

    reCLU = False  # re-cluster flag
    totN = 0

    for inter, count in clu:
        totN += count
    for inter, count in clu:
        if count / float(totN) >= cutoff and count >= readcutoff:
            intervals.append(inter)
            dic[inter] = count  # {(start, end): reads}
        else:
            reCLU = True  # any intron not passing filters will enforce reCLU

    if len(intervals) == 0:
        return []  # base case

    # Below makes sure that after trimming/filtering, the clusters are still good
    # afterwards - each clusters have linked introns that pass filters.

    Atmp, _ = cluster_intervals(intervals)
    A = []

    # A: a list of linked intron clusters
    if len(Atmp) > 0:
        for cl in Atmp:  # Atmp is a list of list
            if len(cl) == 1:  # 1 intron
                A.append(cl)
            for c in refine_linked([(x, 0) for x in cl]):  # >1 introns
                if len(c) > 0:
                    A.append([x[0] for x in c])

    if len(A) == 1:  # A has 1 cluster of introns
        rc = [(x, dic[x]) for x in A[0]]

        if len(rc) > 0:
            if reCLU:  # recompute because ratio changed after removal of some introns
                return refine_cluster(
                    [(x, dic[x]) for x in A[0]], cutoff, readcutoff
                )  # recursive
            else:
                return [[(x, dic[x]) for x in A[0]]]

    NCs = []  # As in N Clusters, here A has more than 1 clusters of introns
    for c in A:
        if len(c) > 1:  # c has more than 1 introns
            NC = refine_cluster(
                [(x, dic[x]) for x in c], cutoff, readcutoff
            )  # recursive
            NCs += NC

    return NCs


def refine_clusters(options):
    """Refine clusters.

    Refine clusters such that kept clusters that are written to file meets
    the following criteria:
        * introns a linked (share either 5' or 3' splice site)
        * minimum total cluster reads cutoff
        * minimum intron reads cutoff
        * minimum reads ratio per cluster cutoff

    However, if constitutive flag `const` is on, then non-linked introns
    are also written out, and they do not subject to cluster ratio and
    cluster reads cutoff filters.

    Parameters:
    -----------
        options : argparse object

    Returns:
    --------
        return : no returns. Use side-effects

    Side-effects:
    -------------
        write refined clusters to file - `*_refined`. Output file is NOT
        version sorted.

    """

    outPrefix = options.outprefix
    rundir = options.rundir
    minratio = float(options.mincluratio)
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    inFile = f"{rundir}/clustering/{outPrefix}_pooled"
    outFile = f"{rundir}/clustering/{outPrefix}_refined"
    fout = open(outFile, "w")

    logger.info(f"Refine clusters from {inFile}...\n")

    Ncl = 0
    for ln in open(inFile):  # pooled juncs
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        clu = []  # each cluster: [((start, end), reads),..]
        totN = 0  # total cluster reads
        chrom = ln.split()[0]
        for ex in ln.split()[1:]:  # for an exon
            A, B, N = ex.split(":")
            clu.append(((int(A), int(B)), int(N)))
            totN += int(N)

        if totN < minclureads:
            continue

        if (
            options.const
        ):  # include constitutive introns. These are clusters that only have 1 intron, hence "constitutive"
            if len(clu) == 1:
                buf = f"{chrom} "
                for interval, count in clu:
                    buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                Ncl += 1
                fout.write(buf + "\n")  # e.g. 'chr:strand start:end:reads'

        for cl in refine_linked(clu):  # only linked intron clusters
            rc = refine_cluster(cl, minratio, minreads)
            if len(rc) > 0:
                for clu in rc:
                    buf = f"{chrom} "
                    for interval, count in clu:
                        buf += f"{interval[0]}:{interval[1]}" + f":{count}" + " "
                    Ncl += 1
                    fout.write(buf + "\n")
    logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Split into {Ncl} clusters...\n")
    fout.close()


def addlowusage(options):
    """Add low usage introns to refined clusters

    Parameters:
    -----------
    options : argparse object
        pass in command options


    Returns:
    --------
    return : null
        no returns. Write files in side-effects.

    Side-effects:
    ------------
        written files:
            - [out_prefix]_lowusage_introns : file stores low usage introns (by cluster).
              Output file is version sorted.
            - [out_prefix]_clusters    : file stores all usage introns (by cluster),
              although each cluster must pass min cluster reads cutoff. Output file is
              version sorted.

    """


    logger.info(f"Add low usage introns...\n")

    outPrefix = options.outprefix
    rundir = options.rundir
    pooled = f"{rundir}/clustering/{outPrefix}_pooled"
    minclureads = int(options.minclureads)
    minreads = int(options.minreads)

    if options.cluster == None:
        refined_cluster = f"{rundir}/clustering/{outPrefix}_refined"
    else:
        refined_cluster = options.cluster

    outFile = (
        f"{rundir}/clustering/{outPrefix}_clusters"  # out file that includes noisy introns
    )
    outFile_lowusageintrons = (
        f"{rundir}/clustering/{outPrefix}_lowusage_introns"  # out file for lowusage introns
    )

    fout = open(outFile, "w")
    fout_lowusage = open(outFile_lowusageintrons, "w")

    # get 5' sites, 5' sites, and clusters of introns from refined file, see data structure below
    exons5, exons3, cluExons = {}, {}, {}
    cluN = 0  # clusterID

    # construct 5' sites, 3' sites, and clusters dict from refined
    for ln in open(refined_cluster):
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        chrom = ln.split()[0]
        cluN += 1
        cluExons[(chrom, cluN)] = []  # keys are (chrom, cluster_number)
        for exon in ln.split()[1:]:
            A, B, count = exon.split(":")  # start, end, reads
            if chrom not in exons5:
                exons5[chrom] = {}
                exons3[chrom] = {}
            exons5[chrom][int(A)] = (
                chrom,
                cluN,
            )  # 5' sites, { k=chrom, v={ k=start, v=(chrom, clusterID) } }
            exons3[chrom][int(B)] = (
                chrom,
                cluN,
            )  # 3' sites, { k=chrom, v={ k=end, v=(chrom, clusterID) } }
            cluExons[(chrom, cluN)].append(
                exon
            )  # introns, { k=(chrom, clusterID), v=['start:end:reads'] }

    # Below for loop adds back clusters (previously filtered out in refined_clusters)
    # in the pooled junc file. These previously removed introns are added to all
    # cluExons, as well as to lowusage_intron
    lowusage_intron = {}  # { k=(chrom, clusterID), v=['start:end:reads'...]}
    for ln in open(pooled):  # each cluster/line from pool_juncs file

        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        clu = []
        totN = 0
        chrom = ln.split()[0]

        if chrom in exons5:  # ensure chrom is in exons5 level-1 keys

            # Below for loop adds introns that were filtered out in refined, aka noisy_introns,
            # back to a total intron cluster dict and to a lowusage (noisy) intron cluster dict
            for exon in ln.split()[1:]:
                A, B, N = exon.split(":")  # start, end, reads

                # when 5' site in refined
                if int(A) in exons5[chrom]:
                    clu = exons5[chrom][
                        int(A)
                    ]  # set clu=(chrom, clusterID), key for cluExons
                    if exon not in cluExons[clu]:  # exon was filtered out by refined
                        cluExons[clu].append(exon)  # add it to cluExons
                        if clu not in lowusage_intron:
                            lowusage_intron[clu] = []
                        lowusage_intron[clu].append(exon)  # also add it to lowusage

                # else when 3' site in refined, perform same procedure
                elif int(B) in exons3[chrom]:  # when 3' site is in refined
                    clu = exons3[chrom][int(B)]
                    if exon not in cluExons[clu]:
                        cluExons[clu].append(exon)
                        if clu not in lowusage_intron:
                            lowusage_intron[clu] = []
                        lowusage_intron[clu].append(exon)

                # neither 5' nor 3' splice site in refined, only add cluster if intron meets minreads requirement
                # because of the minreads requirement, this intron is not noisy, thus do not add to lowusage_intron
                else:
                    if int(N) > minreads:
                        cluN += 1
                        cluExons[(chrom, cluN)] = [exon]

    # write low usage introns
    ks = natural_sort(
        lowusage_intron.keys()
    )  # e.g. { k=(chrom, clusterID), v=['start:end:reads'...]}
    for clu in ks:  # e.g. (chrom, clusterID)
        fout_lowusage.write(clu[0] + " " + " ".join(lowusage_intron[clu]) + "\n")
    fout_lowusage.close()

    # write all intron clusters
    cluLst = natural_sort(cluExons.keys())
    for clu in cluLst:
        if not options.const:  # if -C flag not set, do not write constitutive introns
            if len(cluExons[clu]) == 1:
                continue  # skip write out if only 1 intron in cluster, aka, constitutive

        # only write introns if minimum cluster reads criteria is met
        if sum([int(ex.split(":")[-1]) for ex in cluExons[clu]]) < minclureads:
            continue
        chrom = clu[0]
        buf = f"{chrom}"
        for ex in cluExons[clu]:
            buf += " " + ex
        fout.write(buf + "\n")
    fout.close()


def get_sample_name_from_line(line):
    """Extract sample name from junction file line"""
    parts = line.strip().split('\t')
    filepath = parts[0].strip()
    
    if len(parts) >= 2 and parts[1].strip():
        # Use provided sample name from second column
        return parts[1].strip(), filepath
    else:
        # Use basename approach
        basename = os.path.basename(filepath)
        
        # Remove common file extensions
        for ext in ['.junc.gz', '.bed.gz', '.junc', '.bed', '.gz']:
            if basename.endswith(ext):
                basename = basename[:-len(ext)]
                break
        
        return basename, filepath

def sort_junctions(libl, options):
    """Sort junctions by cluster

    For each intron cluster, sort introns. Write both numerator (intron reads) and
    denominator (cluster reads) into output file.

    Parameters:
    -----------
        libl : str
            A list of junction file paths
        options: argparse object
            Attributes store command line options

    Returns:
    --------
        return : no returns. Use site effect.

    Side-effects:
    -------------
        text file : '{rundir}/{outPrefix}_sortedLibs'
            store junfile names that are processed/sorted
        text file:  '{rundir}/{outPrefix}...sorted.gz'
            a series of sorted input junction files, sorted.
    """


    outPrefix = options.outprefix
    rundir = options.rundir
    checkchrom = options.checkchrom

    if options.cluster == None:  # if not providing refined clusters externally
        refined_cluster = (
            f"{rundir}/clustering/{outPrefix}_clusters"  # note refined noisy intron clusters
        )
        logger.info(f"Using {refined_cluster} as refined cluster...")
    else:
        refined_cluster = options.cluster

    runName = f"{rundir}/{outPrefix}"

    # exons:  { k=chrom : v={ k=(start, end) : v=clusterID } }
    # cluExons:  { k=clusterID : v=[(chrom, start, end)] }
    exons, cluExons = {}, {}
    cluN = 0  # clusterID

    # fill in exons, cluExons dict from `*clusters` intron cluster file
    for ln in open(
        refined_cluster
    ):  # e.g. ln = "chr10:+ 135203:179993:5 135302:179993:29"

        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string

        chrom = ln.split()[0]  # e.g. "chr10:+"
        cluN += 1
        for exon in ln.split()[1:]:  # e.g. "135203:179993:5 135302:179993:29"
            A, B, count = exon.split(":")

            if chrom not in exons:
                exons[chrom] = {}
            if (int(A), int(B)) not in exons[chrom]:
                exons[chrom][(int(A), int(B))] = [cluN]
            else:
                exons[chrom][(int(A), int(B))].append(cluN)
            if cluN not in cluExons:
                cluExons[cluN] = []
            cluExons[cluN].append((chrom, A, B))

    merges = {}  # stores junc file names as dict { k=sample_name : v=[filepath] }

    for line in libl:
        sample_name, filepath = get_sample_name_from_line(line)
        
        if not os.path.isfile(filepath):
            continue
            
        if sample_name not in merges:
            merges[sample_name] = []
        merges[sample_name].append(filepath)

    fout_runlibs = open(
        os.path.join(rundir, outPrefix) + "_sortedlibs", "w"
    )  # intermediate file to store sorted junc file names to be written

    # operate on each sample_name (library), each sample_name can have more than 1+ junc files
    for sample_name in merges:
        by_chrom = {}  # to store junctions from original unsorted junc file

        # write sorted junc file names into intermediate file
        foutName = os.path.join(
            rundir, outPrefix + "_" + sample_name + ".junc.sorted.gz"
        )  # 'test/gtex_w_clu/gtex_GTEX-1IDJU-0006-SM-CMKFK.junc.sorted.gz'
        fout_runlibs.write(foutName + "\n")  # e.g. 'test/gtex_w_clu/gtex_sortedlibs'

        if options.verbose:
            logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Sorting {sample_name}..\n")
        if len(merges[sample_name]) > 1:
            if options.verbose:
                logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} merging {' '.join(merges[sample_name])}...\n")
        else:
            pass
        fout = gzip.open(
            foutName, "wt"
        )  # e.g. 'test/gtex_w_clu/gtex_GTEX-1IDJU-0006-SM-CMKFK.junc.sorted.gz'

        # -------- Process and write junction files --------

        # write header
        fout.write(f"chrom {sample_name}\n")  # 'chrom GTEX-111VG-0526-SM-5N9BW\n'

        # -------- Gather counts from all junc files of library --------
        # store in by_chrom: { ('chr1', '+') : { (100, 300) : 5, (500, 700): 10, ... } }
        for lib in merges[sample_name]:
            if ".gz" in lib:
                F = gzip.open(lib)
            else:
                F = open(lib)

            for ln in F:  # 1 line: e.g. "chr17\t81701131\t81701534\t.\t1\t+"

                if type(ln) == bytes:
                    ln = ln.decode("utf-8")  # convert bytes to string
                lnsplit = ln.split()

                if len(lnsplit) < 6:
                    logger.error(f"Error in {lib} \n")
                    continue

                if len(lnsplit) == 12:
                    (
                        chrom,
                        A,
                        B,
                        dot,
                        counts,
                        strand,
                        rA,
                        rb,
                        rgb,
                        blockCount,
                        blockSize,
                        blockStarts,
                    ) = lnsplit
                    if int(blockCount) > 2:
                        logger.debug("ignored junction with blockCount>2")
                        continue
                    Aoff, Boff = blockSize.split(",")[:2]
                    A, B = int(A) + int(Aoff), int(B) - int(Boff)

                elif len(lnsplit) == 6:
                    # old leafcutter junctions
                    chrom, A, B, dot, counts, strand = lnsplit
                    A, B = int(A), int(B)

                A, B = int(A), int(B) + int(options.offset)  # start, end + offset

                chrom = (chrom, strand)
                if chrom not in by_chrom:
                    by_chrom[chrom] = (
                        {}
                    )  # store introns from junc file, key: ('chr1', '+')

                intron = (A, B)
                if intron in by_chrom[chrom]:  # sum up reads by intron from junc files
                    by_chrom[chrom][intron] += int(counts)
                else:
                    by_chrom[chrom][intron] = int(counts)

        # ------- Take clusters from clusters, assign reads -------
        # reads are from by_chrom (junc files)
        # For each intron cluster, write fraction for each intron (one intron per line).
        for clu in cluExons:  # cluExons: { k=cluID : v=[(chrom, start, end)...]}
            buf = []
            ks = cluExons[clu]  # eg: [('chr1:+', 827776, 829002), ..] introns of a clu
            ks.sort()  # no need to version sort within cluster

            # Step 1: sum cluster level reads from each intron
            # gather (sum) reads for each cluster in clusters, read counts are from junc file (by_chrom)
            tot = 0  # sum of total read counts per cluster
            usages = []
            for exon in ks:
                chrom, start, end = exon
                chrom = tuple(
                    chrom.split(":")
                )  # convert 'chr3:+' to ('chr3', '+') as in by_chrom
                start, end = int(start), int(end)

                if chrom not in by_chrom:
                    pass
                elif (start, end) in by_chrom[chrom]:
                    tot += by_chrom[chrom][(start, end)]

            # Step 2: append intron usage fraction to stream buffer
            for exon in ks:
                chrom, start, end = exon
                start, end = int(start), int(end)
                chrom = tuple(chrom.split(":"))
                chromID, strand = chrom  # chromID eg: 'chr3'

                intron = (
                    chromID,
                    start,
                    end + 1,
                    strand,
                )  # converting to 1-based coordinates

                if chrom not in by_chrom:
                    # if refined exon chrom is not found in junc file, write 0/cluster_total
                    buf.append(f"{chromID}:{start}:{end}:clu_{clu}_{strand} 0/{tot}\n")
                elif (start, end) in by_chrom[chrom]:
                    # if refind exon is in junc file, write exon reads / cluster_total
                    buf.append(
                        f"{chromID}:{start}:{end}:clu_{clu}_{strand} {by_chrom[chrom][(start,end)]}/{tot}\n"
                    )
                else:
                    # if refined exon is not found in junc file, write 0/cluster_total
                    buf.append(f"{chromID}:{start}:{end}:clu_{clu}_{strand} 0/{tot}\n")

            fout.write("".join(buf))
        fout.close()
    fout_runlibs.close()


def merge_files(fnames, fout, options):
    """Merge a list of files into a gzip file

    Parameters:
    -----------
    fnames : list
        list of file path. This is a list of files within a single batch
    fout : a single opened file io
        specifically, a gzip.open('w') io object
    options :
        argparse object

    Returns:
    --------
        return : no returns. Use side-effects.

    Side-effects:
    -------------
        After merging all files from `fnames` list, write out into a gzipped
    file io, as opened in `fout`.
    """

    fopen = []
    for fname in fnames:  # each file to be merged
        if fname[-3:] == ".gz":
            fopen.append(gzip.open(fname))
        else:
            fopen.append(open(fname))

    finished = False
    N = 0
    while not finished:  # cycle through files in batch
        N += 1
        if N % 50000 == 0:
            logger.debug(".")
        buf = []
        for f in fopen:  # each opened file
            ln = f.readline().decode().split()  # read 1 line
            if len(ln) == 0:  # end of line = finish
                finished = True
                break
            chrom = ln[0]  # e.g. "chrom" or "chr1:825552:829002:clu_1_+"
            data = ln[1:]  # e.g. "GTEX-1117F-0626-SM-5N9CS.leafcutter" or "0/0"
            if len(buf) == 0:
                buf.append(chrom)
            buf += data  # e.g. ['chrom', 'GTEX-111VG-0526-SM-5N9BW.leafcutter', 'GTEX-1117F-0626-SM-5N9CS.leafcutter'] for first lines, or ['chr1:825552:829002:clu_1_+', '0/0', '0/0'] for 2+ lines
            # each file the exact same chromosome coordinates, effectively we are collecting counts into columns 2 and after

        if len(buf) > 0:
            if buf[0] == "chrom":
                if options.verbose:
                    logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} merging {len(buf)-1} files")
            fout.write(" ".join(buf) + "\n")  # combining sample counts into columns
        else:
            break

    logger.info(f" done.\n")
    for fin in fopen:
        fin.close()


def merge_junctions(options):
    """Merge junctions

    Merge a list of sorted junction files into a single merged junction file.
    Each input sorted junction files must have the same introns, i.e. first
    column of each row must be the same across all files to be merged.

    Parameters:
    -----------
    options : argparse object

    Returns:
    ---------
    return : null
        No returns. Use side effect.

    Side-effects:
    -------------
        Collect previously sorted junction files. Merge junction files in batches.
    And finally, all batches are merged into a single file. Reads fractions are in
    columns.
        row1  : col1=`chrom`, col2 and beyond are input file names merged
        row2+ : col1=`intron identifier`, reads fraction from each input file
    """

    outPrefix = options.outprefix
    rundir = options.rundir

    fnameout = os.path.join(f"{rundir}/{outPrefix}")
    flist = fnameout + "_sortedlibs"  # sorted juncs file list
    lsts = []  # = flist
    for ln in open(flist):
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string
        lsts.append(ln.strip())

    if options.verbose:
        logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} \nMerging {len(lsts)} junction files...\n")

    # Change 300 if max open file is < 300
    # set up batch N per batch
    N = min([300, max([100, int(len(lsts) ** (0.5))])])

    # tmpfiles = []
    while len(lsts) > 1:  # initial list of sorted junc files

        # convert lsts (list of file paths) to clst (list of lists)
        # each sublist is a batch of upto 100 files.
        clst = []  # list of batches, each batch has up to 100 sorted junc files
        for i in range(
            0, int(len(lsts) / N) + 1
        ):  # merge in batches of max(100, len(lsts))
            lst = lsts[N * i : N * (i + 1)]
            if len(lst) > 0:
                clst.append(lst)
        lsts = (
            []
        )  # clear initial file list, now repurposed to store merged file names (temp)

        for lst in clst:  # run by batch
            if len(lst) == 0:
                continue
            tmpfile = tempfile.mktemp()
            os.mkdir(tmpfile)
            foutname = tmpfile + "/tmpmerge.gz"
            fout = gzip.open(
                foutname, "wt"
            )  # create a temp file for the batch of files to be merged

            merge_files(lst, fout, options)  # merge the batch into `fout`
            lsts.append(foutname)  # save the temp merged file name
            # tmpfiles.append(foutname) # this line is not needed.
            fout.close()

    if not options.const:
        shutil.move(lsts[0], fnameout + "_perind.counts.gz")
    else:
        shutil.move(lsts[0], fnameout + "_perind.constcounts.gz")


def get_numers(options):
    """Get numerators from merged count table

    Parameters:
    -----------
    options : argparse object

    Returns:
    --------
    return : null
        No returns. Use side-effect to write out file.

    Side-effects:
    -------------
        Take in count tables, extract numerators for each sample.

    """

    outPrefix = options.outprefix
    rundir = options.rundir

    if not options.const:
        fname = f"{rundir}/{outPrefix}_perind.counts.gz"
        fnameout = f"{rundir}/{outPrefix}_perind_numers.counts.gz"
    else:
        fname = f"{rundir}/{outPrefix}_perind.constcounts.gz"
        fnameout = f"{rundir}/{outPrefix}_perind_numers.constcounts.gz"

    input_file = gzip.open(fname, "r")
    fout = gzip.open(fnameout, "wt")
    first_line = True

    logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} \nExtracting numerators (read counts) from {fname}...")

    for l in input_file:
        if first_line:
            fout.write(
                " ".join(l.decode().strip().split(" ")[1:]) + "\n"
            )  # print the sample names
            first_line = False
        else:
            l = l.decode().strip()
            words = l.split(" ")
            fout.write(
                words[0] + " " + " ".join([g.split("/")[0] for g in words[1:]]) + "\n"
            )  # write intron and numerators

    input_file.close()
    fout.close()
    logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} done.\n")

def pick_sjc_label(df):
    '''Pick one label for each intron based on priroty'''
    sjc_priority = {'PR': 1, 'UP': 2, 'NE':3} # PR > UP > NE
    df['priority'] = df['SJClass'].map(sjc_priority)
    chosen = df.sort_values(by = 'priority', ascending = True).iloc[0]
    return chosen.drop('priority')


def merge_discordant_logics(sjc_file: str):
    '''some junctions have multiple classifications. Use conservative approach
    to merge them.

    OPTIMIZED: Uses vectorized sort + drop_duplicates instead of slow groupby + apply
    '''

    classifier = {
        # each bit represents:
        # is_GTFAnnotatedCoding?  is_GTFAnnotated?  is_LF2AnnotatedCoding?  is_ClosetoUTR?
        '0000': 'UP', '0001': 'NE', '0010': 'PR', '0011': 'PR',
        '0100': 'UP', '0101': 'NE', '0110': 'PR', '0111': 'PR',
        '1000': 'PR', '1001': 'PR', '1010': 'PR', '1011': 'PR',
        '1100': 'PR', '1101': 'PR', '1110': 'PR', '1111': 'PR'
    }

    classifer_3bits = {
        # each bit represents:
        # is_GTFAnnotated?  is_LF2AnnotatedCoding?  is_ClosetoUTR?
        '000': 'UP', '001': 'NE', '010': 'PR', '011': 'PR',
        '100': 'UP', '101': 'NE', '110': 'PR', '111': 'PR',
    }

    sjc_priority = {'PR': 1, 'UP': 2, 'NE': 3}

    sjc = pd.read_csv(sjc_file, sep="\t")

    # group dt; NOTE:ForwardSpliceJunctionClassifier has an extra Strand column, backward doesn't
    if 'Strand' in sjc.columns:
        sjc = sjc[['Gene_name', 'Intron_coord', 'Strand', 'GencodePC', 'Annot', 'Coding', 'UTR']]

        # convert Annotation, Coding, UTR status to SJ categories
        sjc['SJClass'] = sjc[['GencodePC', 'Annot', 'Coding', 'UTR'
                              ]].astype(int).astype(str).agg(''.join, axis=1).map(classifier)

        # Vectorized priority selection (replaces slow groupby + apply)
        sjc['priority'] = sjc['SJClass'].map(sjc_priority)
        sjc = sjc.sort_values(['Intron_coord', 'Strand', 'priority'])
        sjc = sjc.drop_duplicates(subset=['Intron_coord', 'Strand'], keep='first')
        sjc = sjc[['Intron_coord', 'Strand', 'SJClass', 'Gene_name']]

        # convert df to dict
        sjc = sjc.set_index(['Intron_coord', 'Strand']).to_dict(orient='index')
    else:
        sjc = sjc[['Gene_name', 'Intron_coord', 'Annot', 'Coding', 'UTR']]

        # convert Annotation, Coding, UTR status to SJ categories
        sjc['SJClass'] = sjc[['Annot', 'Coding', 'UTR'
                              ]].astype(int).astype(str).agg(''.join, axis=1).map(classifer_3bits)

        # Vectorized priority selection (replaces slow groupby + apply)
        sjc['priority'] = sjc['SJClass'].map(sjc_priority)
        sjc = sjc.sort_values(['Intron_coord', 'priority'])
        sjc = sjc.drop_duplicates(subset=['Intron_coord'], keep='first')
        sjc = sjc[['Intron_coord', 'SJClass', 'Gene_name']]

        # convert df to dict
        sjc = sjc.set_index('Intron_coord').to_dict(orient='index')

    sjc = {flatten_tuple(k): v for k, v in sjc.items()}

    # sjc is a dictionary with:
    # - keys: intron coordinates, e.g. ('chr1', 1000, 2000, '+') or ('chr1', 1000, 2000) for backward
    # - values: a dictionary e.g. {'SJClass': 'UP', 'Gene_name': 'DNMBP'}
    return sjc


def flatten_tuple(t):
    # t: tuple like ('chr1:100-200', '+')' or str 'chr1:100-200'
    if isinstance(t, tuple):
        c, ab = t[0].split(":")
    elif isinstance(t, str):
        c, ab = t.split(":")
    a, b = ab.split("-")
    a, b = int(a), int(b)
    if isinstance(t, tuple):
        s = t[1]
        return((c, a, b, s)) # e.g. ('chr1', 100, 200, '+')
    else:
        return((c, a, b))

def validate_gtf_requirements(gtf_file, options, CheckRequirementsEveryNLines=None):
    """
    Validate GTF file and auto-detect attribute names if not specified by user.
    Optionally perform periodic checks (every N lines) to exit early when all
    required features/attributes are found.
    
    Args:
        gtf_file: Path to GTF file
        options: argparse options object
        CheckRequirementsEveryNLines: int or None. If int, check every N lines
            for completion and break early when satisfied. If None, parse full file.
    
    Returns:
        options: Modified options object with auto-detected attributes
    
    Raises:
        SystemExit if required features/attributes are missing
    """
    
    # Attribute alternatives in priority order
    attribute_alternatives = {
        'gene_name': ['gene_name', 'gene_id', 'gene_symbol'],
        'transcript_name': ['transcript_name', 'transcript_id'], 
        'transcript_type': ['transcript_type', 'transcript_biotype']
    }
    
    # Required feature types
    required_features = {'transcript','gene','exon', 'start_codon', 'stop_codon', 'CDS'}
    
    # Track what we find
    found_features = set()
    found_attributes = set()
    protein_coding_found = False
    
    logger.info(f"Validating GTF file: {gtf_file}\n")
    
    # Helper: are all attribute alternatives present?
    def _attributes_satisfied():
        for _, alts in attribute_alternatives.items():
            if not any(a in found_attributes for a in alts):
                return False
        return True
    
    # Parse GTF file
    try:
        if gtf_file.endswith('.gz'):
            file_handle = gzip.open(gtf_file, 'rt')
        else:
            file_handle = open(gtf_file, 'r')
            
        line_count = 0
        for line in file_handle:
            line_count += 1
            if line_count % 100000 == 0:
                logger.info(f"  Processed {line_count} lines...\n")
                
            # Skip comments and empty lines
            if line.startswith('#') or not line.strip():
                continue
                
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
                
            feature_type = fields[2]
            attributes_str = fields[8]
            
            # Track feature types
            found_features.add(feature_type)
            
            # Parse attributes
            attributes = {}
            for attr_pair in attributes_str.split(';'):
                attr_pair = attr_pair.strip()
                if not attr_pair:
                    continue
                    
                # Handle both quoted and unquoted values
                if '"' in attr_pair:
                    # Quoted format: gene_name "WASH7P"
                    parts = attr_pair.split('"')
                    if len(parts) >= 2:
                        key = parts[0].strip()
                        value = parts[1]
                        attributes[key] = value
                else:
                    # Unquoted format: gene_name WASH7P
                    parts = attr_pair.split(None, 1)
                    if len(parts) == 2:
                        key, value = parts
                        attributes[key] = value
            
            # Track all found attributes
            found_attributes.update(attributes.keys())
            
            # Check for protein_coding
            for attr_name in ['transcript_type', 'transcript_biotype']:
                if attr_name in attributes and attributes[attr_name] == 'protein_coding':
                    protein_coding_found = True
            
            # Optional periodic early check
            if CheckRequirementsEveryNLines and (line_count % int(CheckRequirementsEveryNLines) == 0):
                if required_features.issubset(found_features) and _attributes_satisfied():
                    logger.info(f"  Early stop at {line_count} lines: requirements satisfied.\n")
                    break
        file_handle.close()
        
    except Exception as e:
        logger.error(f"Error reading GTF file {gtf_file}: {e}\n")
        exit(1)
    
    logger.info(f"GTF validation complete. Processed {line_count} lines.\n")
    
    # Validate required features
    missing_required = required_features - found_features
    
    if missing_required:
        logger.error("Error: Missing required feature types in GTF:\n")
        logger.error(f"  Required but missing: {sorted(missing_required)}\n")
        logger.error(f"  Found feature types: {sorted(found_features)}\n")
        logger.error("  Required features: gene, exon, start_codon, stop_codon, CDS\n")
        logger.error("  Note: UTR classification is determined from start/stop codons, not UTR features\n")
        raise SystemExit(1)
    
    # Check for protein_coding transcripts
    if not protein_coding_found:
        logger.warning("Warning: No 'protein_coding' transcript types found in GTF.\n")
        logger.warning("  This may affect junction classification.\n")
    
    # Auto-detect and validate attributes
    user_specified = {}
    auto_detected = {}
    
    for option_name, alternatives in attribute_alternatives.items():
        current_value = getattr(options, option_name)
        
        if current_value and current_value != "":
            # User specified - validate it exists
            user_specified[option_name] = current_value
            if current_value not in found_attributes:
                logger.error(f"Error: User-specified attribute '{current_value}' for {option_name} not found in GTF\n")
                logger.error(f"Available attributes: {sorted(found_attributes)}\n")
                raise SystemExit(1)
        else:
            # Auto-detect from alternatives
            detected = None
            for alternative in alternatives:
                if alternative in found_attributes:
                    detected = alternative
                    break
                    
            if detected:
                auto_detected[option_name] = detected
                setattr(options, option_name, detected)
            else:
                logger.error(f"Error: No suitable attribute found for {option_name}\n")
                logger.error(f"  Looked for: {alternatives}\n") 
                logger.error(f"  Available attributes: {sorted(found_attributes)}\n")
                logger.error(f"  Suggestion: Your GTF may use non-standard attribute names\n")
                raise SystemExit(1)
    
    # For non-required attributes (anything other than gene_id/transcript_id), verify that
    # the chosen attribute is actually populated in every record across the full GTF.
    # gene_id and transcript_id are required by the GTF spec and can be trusted without a
    # full scan. If any records are missing a non-required attribute, fall back to the safe
    # required alternative and warn the user.
    safe_defaults = {'gene_name': 'gene_id', 'transcript_name': 'transcript_id'}
    attrs_to_check = {
        opt: getattr(options, opt)
        for opt in ('gene_name', 'transcript_name')
        if getattr(options, opt) != safe_defaults[opt]
    }

    if attrs_to_check:
        logger.info(
            f"Performing full GTF scan to verify completeness of non-required attributes: "
            f"{list(attrs_to_check.values())}\n"
        )
        null_counts = {opt: 0 for opt in attrs_to_check}
        total_data_lines = 0
        try:
            opener = gzip.open if gtf_file.endswith('.gz') else open
            with opener(gtf_file, 'rt') as fh:
                for line in fh:
                    if line.startswith('#') or not line.strip():
                        continue
                    fields = line.strip().split('\t')
                    if len(fields) < 9:
                        continue
                    total_data_lines += 1
                    attrs_str = fields[8]
                    parsed = {}
                    for attr_pair in attrs_str.split(';'):
                        attr_pair = attr_pair.strip()
                        if not attr_pair:
                            continue
                        if '"' in attr_pair:
                            parts = attr_pair.split('"')
                            if len(parts) >= 2:
                                parsed[parts[0].strip()] = parts[1]
                        else:
                            parts = attr_pair.split(None, 1)
                            if len(parts) == 2:
                                parsed[parts[0]] = parts[1]
                    for opt, attr in attrs_to_check.items():
                        if not parsed.get(attr, '').strip():
                            null_counts[opt] += 1
        except Exception as e:
            logger.warning(f"Could not verify attribute completeness: {e}\n")

        for opt, attr in list(attrs_to_check.items()):
            if null_counts[opt] > 0:
                fallback = safe_defaults[opt]
                logger.warning(
                    f"Attribute '{attr}' missing or empty in {null_counts[opt]}/{total_data_lines} "
                    f"GTF records; falling back to required attribute '{fallback}'\n"
                )
                setattr(options, opt, fallback)
                if opt in auto_detected:
                    auto_detected[opt] = fallback
                elif opt in user_specified:
                    user_specified[opt] = fallback

    # Print summary of what was used
    if user_specified:
        logger.info("Using user-specified GTF attributes:\n")
        for option_name, value in user_specified.items():
            logger.info(f"  {option_name}: {value}\n")

    if auto_detected:
        logger.info("Auto-detected GTF attributes:\n")
        for option_name, value in auto_detected.items():
            logger.info(f"  {option_name}: {value} (auto-detected)\n")

    logger.info("GTF validation and attribute detection complete.\n")

    return options

def has_cds_feature(gtf_file: str) -> bool:
    """Detect if any CDS feature exists in the GTF (streaming)."""
    try:
        opener = gzip.open if gtf_file.endswith('.gz') else open
        with opener(gtf_file, 'rt') as fh:
            for line in fh:
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) >= 3 and parts[2] == 'CDS':
                    return True
    except Exception as e:
        logger.warning(f"Failed scanning GTF for CDS: {e}")
    return False


def validate_or_reformat_gtf(gtf_file, options):
    """
    Validate GTF; if validation fails and auto-reformat is enabled, attempt
    to reformat with Reformat_gtf and re-validate. Updates options.annot.
    """
    try:
        return validate_gtf_requirements(gtf_file, options, CheckRequirementsEveryNLines=1000)
    except SystemExit:
        if getattr(options, 'no_auto_reformat', False):
            logger.error("GTF validation failed and auto-reformat is disabled.")
            raise
        if not options.genome:
            logger.error("GTF validation failed and genome fasta is required to attempt auto-reformat.")
            raise
        base = os.path.basename(gtf_file)
        base = base[:-3] if base.endswith('.gz') else base
        base = base[:-4] if base.endswith('.gtf') else base
        reformatted = os.path.join(options.rundir, f"reformatted_{base}.gtf")
        cds_present = has_cds_feature(gtf_file)
        trans_approach = 'B' if cds_present else 'D'
        arg_str = (
            f"-i {gtf_file} -fa {options.genome} "
            f"-extra_attributes {options.gene_name},{options.transcript_name},{options.transcript_type} "
            f"-transcript_id_attribute_name {options.transcript_name or 'transcript_id'} "
            f"-gene_id_attribute_name {options.gene_name or 'gene_id'} "
            f"-infer_gene_type_approach B -infer_transcript_type_approach C "
            f"-translation_approach {trans_approach} -o {reformatted}"
        )
        logger.info(f"Attempting to reformat GTF (CDS={'yes' if cds_present else 'no'}) -> {reformatted}")
        try:
            logger.warning("Reformatting gtf is experimental; please check output GTF carefully.")
            logger.warning("Running Transcript_tools to reformat gtf with arguments:\n" + arg_str)
            # Call main with a list of args
            Transcript_tools.main(shlex.split(arg_str))
        except SystemExit as e:
            logger.error(f"Transcript_tools exited with status {e.code}")
            raise SystemExit(1)
        except Exception as e:
            logger.error(f"Transcript_tools failed: {e}")
            raise SystemExit(1)
        options.annot = reformatted
        return validate_gtf_requirements(options.annot, options, CheckRequirementsEveryNLines=1000)

def annotate_noisy(options):
    """Annotate introns

    Produces 2 files. Details in side-effects.

    Parameters:
    -----------
        options : argparse object

    Returns:
    --------
        return : no return. Use side-effects to write output files.

    Side-effects:
    -------------
        noisy counts by intron : output file.
            Same count table as the output of sorted juncs count table,
            except that each intron is annotated with `UP`, `PR`, `NE`,
            or `IN` to denote unporductive, productive, neither, or intergenic.

        noisy numerators : output result file.
            Same as noisy counts, except here only reports the numerators.
    """

    outPrefix = options.outprefix
    rundir = options.rundir
    minreadstd = float(options.minreadstd)
    fnameout = f"{rundir}/{outPrefix}"
    sjc_f = f"{rundir}/clustering/{outPrefix}_junction_classifications.txt"  # Classified junction annotations

    logger.info(
        f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotating introns with custom-classified annotations: {sjc_f} ...\n"
    )

    if sjc_f != None:

        if options.verbose:
            logger.info(
                f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loading {sjc_f} for (un/)productive splicing classification..\n"
            )
        sjc = merge_discordant_logics(sjc_f)

        if options.verbose:
            logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Loaded.\n")

        # Use the counts file source (either provided or generated)
        if options.counts_file:
            fname = options.counts_file
        else:
            if not options.const:
                fname = fnameout + "_perind.counts.gz"  # no constitutive introns
            else:
                fname = fnameout + "_perind.constcounts.gz"

        # If using provided counts file, adjust output paths to use rundir
        if options.counts_file:
            noisydiag = f"{rundir}/{outPrefix}.cluster_ratios.gz"
            numersdiag = f"{rundir}/{outPrefix}.junction_counts.gz"
        else:
            if not options.const:
                noisydiag = fname.replace(
                    "_perind.counts", ".cluster_ratios"
                )  # eg: run/out_perind.counts.classified.gz
                numersdiag = fname.replace(
                    "_perind.counts", ".junction_counts"
                )  # eg: run/out_perind_numers.counts.noise.gz
            else:
                noisydiag = fname.replace(
                    "_perind.constcounts", ".cluster_ratios.const"
                )  # eg: run/out_perind.counts.classified.gz
                numersdiag = fname.replace(
                    "_perind.constcounts", ".junction_counts.const"
                )  # eg: run/out_perind_numers.counts.noise.gz

    foutdiag = gzip.open(noisydiag, "wt")
    foutdiagnumers = gzip.open(numersdiag, "wt")

    F = gzip.open(fname)
    ln = F.readline().decode()
    foutdiag.write(ln)
    foutdiagnumers.write(ln)
    
    # Get number of samples from header line (subtract 1 for the 'chrom' column)
    num_samples = len(ln.split()) - 1

    N_introns_annotated = 0
    N_skipped_introns = 0

    # Determine key format ONCE before the loop (avoids O(n*m) complexity)
    sample_key = next(iter(sjc.keys()))
    use_strand = len(sample_key) == 4

    for ln in F:
        if type(ln) == bytes:
            ln = ln.decode("utf-8")  # convert bytes to string
        ln = ln.split()
        intron = ln[0]
        chrom, s, e, clu = intron.split(":")  # chr, start, end, clu_1_+
        strand = clu.split("_")[-1]

        # check if Strand is in sjc keys:
        if use_strand:
            intronid = chrom, int(s), int(e), strand
        else:
            intronid = chrom, int(s), int(e)

        fractions = [x.split("/") for x in ln[1:]]
        usages = [int(f[0]) / (float(f[1]) + 0.1) for f in fractions] # intron usage ratios
        reads = [int(f[0]) for f in fractions]  # numerators
        if num_samples > 1:
            sdreads = stdev(reads)  # standard deviation of read counts across samples
        else:
            sdreads = 0
            minreadstd = -1

        # remove intron if read count SD < 0.5 and usage ratios are all 0
        if sum(usages) == 0 or sdreads < minreadstd:
            N_skipped_introns += 1
            continue

        # annotate using custom classification
        if intronid in sjc:
            classification = sjc[intronid]["SJClass"]
        else:
            classification = "IN"  # IN: INtergenic

        # add class flag and write to *_perind.noise_by_intron.gz, eg: 'chr1:825552:829002:clu_1_+:F 1/14 0/25 1/33 1/14 1/33'
        foutdiag.write(intron + f":{classification}" + " " + " ".join(ln[1:]) + "\n")
        foutdiagnumers.write(
            intron
            + f":{classification}"
            + " "
            + " ".join([str(x) for x in reads])
            + "\n"
        )
        
        N_introns_annotated += 1
        if N_introns_annotated % 10000 == 0:
            logger.info(f" ... {N_introns_annotated} introns annotated.\n")

    foutdiag.close()
    foutdiagnumers.close()

    logger.info(f"{datetime.now().strftime('%Y-%m-%d %H:%M:%S')} Annotation done.\n")
    logger.info(f"Annotated {N_introns_annotated} introns.\n")
    logger.info(f"Filtered out {N_skipped_introns} introns with stdev(reads) < {minreadstd} or zero usage.\n")


def main(options, libl):
    
    if not os.path.exists(options.rundir):
        os.mkdir(options.rundir)
    os.makedirs(f"{options.rundir}/clustering/", exist_ok=True)

    # Determine perind_file source
    if options.counts_file:
        # User provided pre-existing counts file
        perind_file = options.counts_file
        logger.info(f"Using provided counts file: {perind_file}\n")
    else:
        # Run full pipeline to generate counts
        if options.cluster == None:
            pool_junc_reads(libl, options)
            refine_clusters(options)
            addlowusage(options)

        sort_junctions(libl, options)
        merge_junctions(options)
        get_numers(options)
        
        # Set perind_file to generated file
        if not options.const:
            perind_file = f"{options.rundir}/{options.outprefix}_perind.counts.gz"
        else:
            perind_file = f"{options.rundir}/{options.outprefix}_perind.constcounts.gz"

    if options.annot != None and options.genome != None:

        # Validate GTF and auto-detect attributes
        options = validate_or_reformat_gtf(options.annot, options)

        logger.info(f"Loading genome {options.genome} ...")
        fa = pyfastx.Fasta(options.genome)
        logger.info("done!")

        logger.info("Classifying splice junctions...")
        sjcf.ClassifySpliceJunction(
            perind_file=perind_file,
            gtf_annot=options.annot,
            fa = fa,
            rundir=options.rundir,
            outprefix=options.outprefix,
            max_juncs=options.max_juncs,
            keepannot=options.keepannot,
            transcript_type=options.transcript_type,
            gene_name=options.gene_name,
            transcript_name=options.transcript_name,
            verbose=options.verbose,
        )
        annotate_noisy(options)

    else:
        logger.info("Skipping annotation step...\n")



def main_cli():

    parser = argparse.ArgumentParser()
    
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        "-j",
        "--juncfiles",
        dest="juncfiles",
        type=str,
        help="a text file storing paths to junction files, one path per line",
    )

    input_group.add_argument(
        "--leafcutter1-counts-file",
        dest="counts_file",
        default=None,
        help="Pre-existing counts file (*.counts.gz or *.constcounts.gz). If provided, skip clustering and counting steps and proceed directly to classification.",
    )

    parser.add_argument(
        "-o",
        "--outprefix",
        dest="outprefix",
        default="leafcutter2",
        help="output prefix (default leafcutter2)",
    )

    parser.add_argument(
        "-q",
        "--quiet",
        dest="verbose",
        default=True,
        action="store_false",
        help="don't print status messages to stdout",
    )

    parser.add_argument(
        "-r",
        "--rundir",
        dest="rundir",
        default="./",
        help="write to directory (default ./)",
    )

    parser.add_argument(
        "-l",
        "--maxintronlen",
        dest="maxintronlen",
        default=100000,
        help="maximum intron length in bp (default 100,000bp)",
    )

    parser.add_argument(
        "-m",
        "--minclureads",
        dest="minclureads",
        default=30,
        help="minimum reads in a cluster (default 30 reads)",
    )

    parser.add_argument(
        "-M",
        "--minreads",
        dest="minreads",
        default=5,
        help="minimum reads for a junction to be considered for \
                        clustering(default 5 reads)",
    )

    parser.add_argument(
        "-D",
        "--minreadstd",
        dest="minreadstd",
        default=0.5,
        help="minimum standard deviation of reads across samples for a \
                        junction to be included in output (default 0.5)",
    )

    parser.add_argument(
        "-p",
        "--mincluratio",
        dest="mincluratio",
        default=0.001,
        help="minimum fraction of reads in a cluster that support a junction \
              (default 0.001)",
    )

    parser.add_argument(
        "-c",
        "--cluster",
        dest="cluster",
        default=None,
        help="skip intron clustering step, use pre-determined clusters",
    )

    parser.add_argument(
        "-k",
        "--checkchrom",
        dest="checkchrom",
        action="store_true",
        default=False,
        help="check that the chromosomes are well formated e.g. chr1, chr2, \
              ..., or 1, 2, ...",
    )

    parser.add_argument(
        "-C",
        "--includeconst",
        dest="const",
        action="store_true",
        default=False,
        help="also include constitutive introns",
    )

    parser.add_argument(
        "-A",
        "--annot",
        dest="annot",
        default=None,
        help="Gencode GTF annotation file, e.g. gencode.v37.annotation.gtf.gz",
    )

    parser.add_argument(
        "-G",
        "--genome",
        dest="genome",
        default=None,
        help="Genome fasta file, e.g. hg38.fa",
    )

    parser.add_argument(
        "-f",
        "--offset",
        dest="offset",
        default=0,
        help="Offset sometimes useful for off by 1 annotations. (default 0)",
    )

    parser.add_argument(
        "-T",
        "--keeptemp",
        dest="keeptemp",
        action="store_true",
        default=False,
        help="keep temporary files. (default false)",
    )
    
    parser.add_argument(
        "-L",
        "--keepleafcutter1",
        dest="keepleafcutter1",
        action="store_true",
        default=False,
        help="keep temporary LeafCutter1 files. Useful for running differential splicing \
              analysis with leafcutter's R package. (default false)",
    )
    
    parser.add_argument(
        "-P",
        "--keepannot",
        dest="keepannot",
        action="store_true",
        default=False,
        help="save parsed annotations to .pckle files. (default false)",
    )
    
    parser.add_argument(
        "-t",
        "--transcript_type",
        dest="transcript_type",
        default="",
        help="tag for transcript type in GTF file (default: auto-detect from transcript_type, transcript_biotype)",
    )

    parser.add_argument(
        "-gn",
        "--gene_name",
        dest="gene_name",
        default="",
        help="tag for gene name or ID in GTF file (default: auto-detect from gene_name, gene_id, gene_symbol)",
    )

    parser.add_argument(
        "-tn",
        "--transcript_name",
        dest="transcript_name",
        default="",
        help="tag for transcript name or ID in GTF file (default: auto-detect from transcript_name, transcript_id)",
    )

    parser.add_argument(
        "-N",
        "--max_juncs",
        dest="max_juncs",
        metavar='N',
        default=10000,
        type=int,
        help="skip solveNMD function if gene contains more than N juncs. Juncs in skipped genes are assigned Coding=False. Default 10000")

    parser.add_argument(
        "--log-level",
        dest="log_level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default INFO)",
    )

    parser.add_argument(
        "--no-auto-reformat-gtf",
        dest="no_auto_reformat",
        action="store_true",
        default=False,
        help="Disable automatic GTF reformatting when validation fails",
    )

    options = parser.parse_args()

    # Initialize logging once
    setup_logging(options.log_level, options.verbose)

    # Validate input arguments
    if not options.counts_file and not options.juncfiles:
        logger.error("Error: Either --juncfiles or --counts-file must be provided\n")
        exit(1)
    
    if options.counts_file:
        if not os.path.exists(options.counts_file):
            logger.error(f"Error: Counts file {options.counts_file} does not exist\n")
            exit(1)
        if not (options.annot and options.genome):
            logger.error("Error: When using --counts-file, both --annot and --genome are required for classification\n")
            exit(1)
    elif options.juncfiles == None:
        logger.error("Error: no junction file provided...\n")
        exit(0)

    # Get the junction file list with optional sample names
    libl = []
    if options.juncfiles:
        for line_num, line in enumerate(open(options.juncfiles), 1):
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            
            # Extract filepath for validation
            filepath = line.split('\t')[0].strip()
            
            try:
                open(filepath)
            except:
                logger.error(f"Error on line {line_num}: {filepath} does not exist... check your junction files.\n")
                exit(0)
            
            libl.append(line)  # Keep the full line (with optional sample name)

    main(options, libl)

    if not options.keeptemp:
        logger.info("Remove generated temp files... \n")
        
        sortedlibs_file = os.path.join(options.rundir, options.outprefix) + "_sortedlibs"
        
        # Check if the sortedlibs file exists before trying to process it
        if os.path.exists(sortedlibs_file):
            try:
                with open(sortedlibs_file) as f:
                    for tmp in [ln.strip() for ln in f.readlines()]:
                        if os.path.exists(tmp):
                            try:
                                os.remove(tmp)
                                if options.verbose:
                                    logger.info(f"Removed {tmp}")
                            except Exception as e:
                                logger.warning(f"Warning: Could not remove {tmp}: {e}\n")
                        else:
                            if options.verbose:
                                logger.info(f"File {tmp} already removed or doesn't exist\n")
                
                # Remove the sortedlibs file itself
                os.remove(sortedlibs_file)
            except Exception as e:
                logger.warning(f"Warning: Could not process {sortedlibs_file}: {e}\n")
        else:
            logger.warning(f"{sortedlibs_file} not found\n")
        
        # Clean up clustering files if they exist
        if options.cluster == None:
            pooled_file = f"{options.rundir}/clustering/{options.outprefix}_pooled"
            refined_file = f"{options.rundir}/clustering/{options.outprefix}_refined"
            
            for cleanup_file in [pooled_file, refined_file]:
                if os.path.exists(cleanup_file):
                    try:
                        os.remove(cleanup_file)
                        if options.verbose:
                            logger.info(f"Removed {cleanup_file}")
                    except Exception as e:
                        logger.warning(f"Warning: Could not remove {cleanup_file}: {e}\n")
        
        logger.info(f"Done.\n")
        
    if (options.annot == None) or (options.genome == None):
        if not options.const:
            shutil.move(os.path.join(options.rundir, options.outprefix) + "_perind.counts.gz", 
                        os.path.join(options.rundir, options.outprefix) + ".cluster_ratios.unclassified.gz")
            shutil.move(os.path.join(options.rundir, options.outprefix) + "_perind_numers.counts.gz", 
                        os.path.join(options.rundir, options.outprefix) + ".junction_counts.unclassified.gz")
        else:
            shutil.move(os.path.join(options.rundir, options.outprefix) + "_perind.constcounts.gz", 
                        os.path.join(options.rundir, options.outprefix) + ".cluster_ratios.unclassified.gz")
            shutil.move(os.path.join(options.rundir, options.outprefix) + "_perind_numers.constcounts.gz", 
                        os.path.join(options.rundir, options.outprefix) + ".junction_counts.unclassified.gz")
    else:
        
        if not options.keepleafcutter1:
            logger.info("Remove generated LeafCutter1 files... \n")
            
            # Only remove files if they were generated (not using provided counts file)
            if not options.counts_file:
                perind_file = os.path.join(options.rundir, options.outprefix) + "_perind.counts.gz"
                numers_file = os.path.join(options.rundir, options.outprefix) + "_perind_numers.counts.gz"
                
                if os.path.exists(perind_file):
                    os.remove(perind_file)
                    if options.verbose:
                        logger.info(f"Removed {perind_file}")
                
                if os.path.exists(numers_file):
                    os.remove(numers_file)
                    if options.verbose:
                        logger.info(f"Removed {numers_file}")
            else:
                if options.verbose:
                    logger.info("Skipping LeafCutter1 file removal - used provided counts file\n")
                    
            logger.info(f"Done.\n")
        else:
            # Only move files if they were generated (not using provided counts file)
            if not options.counts_file:
                os.makedirs(os.path.join(options.rundir, "leafcutter1_files"), exist_ok=True)
                
                perind_file = os.path.join(options.rundir, options.outprefix) + "_perind.counts.gz"
                numers_file = os.path.join(options.rundir, options.outprefix) + "_perind_numers.counts.gz"
                
                if os.path.exists(perind_file):
                    shutil.move(perind_file, 
                                os.path.join(options.rundir, "leafcutter1_files", options.outprefix) + "_perind.counts.gz")
                
                if os.path.exists(numers_file):
                    shutil.move(numers_file, 
                                os.path.join(options.rundir, "leafcutter1_files", options.outprefix) + "_perind_numers.counts.gz")
            else:
                if options.verbose:
                    logger.info("Skipping LeafCutter1 file preservation - used provided counts file\n")

# python ../scripts/leafcutter2.py -r example/outtut_dir -A annotation/chr10.gtf.gz -r output_dir -G annotation/chr10.fa.gz -j junction_files.txt -L -gn gene_id tn transcript_id


if __name__ == "__main__":
    main_cli()
