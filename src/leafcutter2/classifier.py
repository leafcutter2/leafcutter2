#!/bin/env python3
# SpliceJunctionClassified V0.2 (Updated Apr 2025)
# Written by Yang Li Nov-2023
# NOTE: 2024-08-30 This version takes in GTF file and perind_file from previous leafcutter2 steps. 
# GTF coordinates are kept as 1-based like GTF. Perind count files are BED formatted, and converted into leafcutter 1 format immediately after having read in the ClassifierClassifySpliceJunction( function. (chrom, BEDstart, BEDend) to (chrom, BEDstart, BEDend +1). 
# After solving for NMD, intron junctions are revereted back to chrom, BEDstart, BEDend before writing into text output.


import argparse
import gzip
import hashlib
import pickle
import sys
import logging
from bisect import insort
from statistics import mean, median, mode, multimode

import pyfastx
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

def tx_by_gene(gtf_annot, 
               gene_type: str, 
               transcript_type: str,
               gene_name: str, 
               transcript_name: str
              ):
    transcripts_by_gene = {}
    for dic in parse_gtf(gtf_annot, gene_type, transcript_type, gene_name, transcript_name):
        if dic['type'] == 'transcript':
            if dic['gene_name'] not in transcripts_by_gene:
                transcripts_by_gene[dic['gene_name']] = {dic['transcript_name']: []}
            else:
                transcripts_by_gene[dic['gene_name']] = transcripts_by_gene[dic['gene_name']] | {dic['transcript_name']: []}
        if dic['type'] == 'CDS':
            insort(transcripts_by_gene[dic['gene_name']][dic['transcript_name']], int((dic['start'])))
            insort(transcripts_by_gene[dic['gene_name']][dic['transcript_name']], int((dic['end'])))
    return transcripts_by_gene

def NMD_tx(gtf_annot, 
           gene_type: str, 
           transcript_type: str,
           gene_name: str, 
           transcript_name: str
          ):
    nmd_tx_by_gene = {}
    for dic in parse_gtf(gtf_annot, gene_type, transcript_type, gene_name, transcript_name):
        if dic['transcript_type'] != 'nonsense_mediated_decay':
            continue
        if dic['type'] == 'transcript':
            if dic['gene_name'] not in nmd_tx_by_gene:
                nmd_tx_by_gene[dic['gene_name']] = {dic['transcript_name']: []}
            else:
                nmd_tx_by_gene[dic['gene_name']] = nmd_tx_by_gene[dic['gene_name']] | {dic['transcript_name']: []}
        if dic['type'] == 'CDS':
            insort(nmd_tx_by_gene[dic['gene_name']][dic['transcript_name']], int((dic['start'])))
            insort(nmd_tx_by_gene[dic['gene_name']][dic['transcript_name']], int((dic['end'])))
    return nmd_tx_by_gene

def ptc_pos_from_prot(prot, sub):
    to_return = []
    start = 0
    while True:
        start = prot.find(sub, start)
        if start == -1: return to_return
        else:
            to_return.append(start)
        start += 1


def nucleotide_rule(failing_juncs, gene_name, transcripts_by_gene, strand, chrom, nmd_tx_by_gene, fa, exonLcutoff = 2000):
    distances = []
    nuc_rule = []
    unique_juncs_pre_ptc = []
    last_exon_length = {}
    for junc in failing_juncs:
        junc_break = False
        distances_to_ejc = []
        last_exon_lengths = []
        possible_transcripts = transcripts_by_gene[gene_name]
        if gene_name in nmd_tx_by_gene.keys():
            nmd_transcripts = nmd_tx_by_gene[gene_name]
        else:
            nmd_transcripts = {}

        for transcript in possible_transcripts:
            if transcript in nmd_transcripts.keys():
                nmd_transcript = True
            else:
                nmd_transcript = False

            if nmd_transcript: #remove junctions in NMD transcripts completely
                junc_break = True
                break

            s = list(possible_transcripts[transcript])
            #cast to list of integers
            if s == []:
                continue

            new_junc = False
            failing_junc = False
            for i in range(len(s)-2, 0, -2):
                junction = tuple([s[i-1], s[i]])

                #overlapping 3' or 5' splice site, just replace the junction
                if junc[0] == junction[0] or junc[1] == junction[1]:
                    s[i-1] = junc[0]
                    s[i] = junc[1]
                    new_junc = True
                #if this does not fit in our transcript, move on
                    if not (s[i-1] > s[i-2] and s[i] < s[i+1]):
                        failing_junc = True

                    break  

            #If our overlapping site configuration fails in the transcript, skip this transcript     
            if failing_junc or not new_junc:
                continue

            if strand == '-':
                s.reverse()
            allprot = Seq("")
            leftover = Seq("")

            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = int(exon_coord[1])-int(exon_coord[0])
                if exlen > exonLcutoff:
                    break

                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                we must add 1 to length"""
                endpos = (len(leftover)+exlen+1)%3
                if junc not in unique_juncs_pre_ptc:
                    unique_juncs_pre_ptc.append(junc)
                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                        leftover = Seq("")
                    else:
                        prot = seq[:-endpos].translate()
                        leftover = seq[-endpos:]                                                                                                               


                    bool_ptc = "*" in prot 
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot) 

                    allprot = allprot+prot
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()

                    bool_ptc = "*" in prot
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot) 

                    allprot = allprot+prot
                
                if bool_ptc:
                    ptc_coord = [(x + 1)*3 for x in ptc_pos]
                    if (i != len(s)-2 and i != len(s) - 4): #or (nmd_transcript and i != len(s)-4 and ptc_prot_len*3-min(ptc_coord) < 51):
                        junc_break = True
                        break
                    else:
                        if i == len(s) - 2:
                            distances_to_ejc.append(-1)
                            last_exon_lengths.append(exlen)
                            break
                        else:
                            distances_to_ejc.append(ptc_prot_len*3-min(ptc_coord))
                            break
                
            if junc_break:
                break
        if not junc_break and len(distances_to_ejc) != 0:
            if max(distances_to_ejc) > 50: #dealing with PTCs in 2nd to last exon not within 50nt of ejc for plotting
                nuc_rule.append(junc)
                distances.append(max(distances_to_ejc))
            else:
                nuc_rule.append(junc)
                distances.append(max(multimode(distances_to_ejc)))

    return nuc_rule, distances

def many_junctions(failing_juncs, gene_name, transcripts_by_gene, strand, chrom, fa, exonLcutoff = 2000):

    before = {}
    after = {}
    for junc in failing_juncs:
        junc_added = False
        possible_transcripts = transcripts_by_gene[gene_name]
        possible_exons_before = []
        possible_exons_after = []
        for transcript in possible_transcripts:
            s = list(possible_transcripts[transcript])
            #cast to list of integers
 

            new_junc = False
            failing_junc = False
            for i in range(len(s)-2, 0, -2):
                junction = tuple([s[i-1], s[i]])

                #overlapping 3' or 5' splice site, just replace the junction
                if junc[0] == junction[0] or junc[1] == junction[1]:
                    s[i-1] = junc[0]
                    s[i] = junc[1]
                    new_junc = True
                #if this does not fit in our transcript, move on
                    if not (s[i-1] > s[i-2] and s[i] < s[i+1]):
                        failing_junc = True
                    
                    break  
            #If our overlapping site configuration fails in the transcript, skip this transcript     
            if failing_junc or not new_junc:
                continue

            if strand == '-':
                s.reverse()
            allprot = Seq("")
            leftover = Seq("")

            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = int(exon_coord[1])-int(exon_coord[0])
                if exlen > exonLcutoff:
                    break

                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                we must add 1 to length"""
                endpos = (len(leftover)+exlen+1)%3

                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                        leftover = Seq("")
                    else:
                        prot = seq[:-endpos].translate()
                        leftover = seq[-endpos:]                                                                                                                


                    bool_ptc = "*" in prot 
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot)
                        ptc_coord = (min(ptc_pos) + 1)*3


                    allprot = allprot+prot
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()

                    bool_ptc = "*" in prot
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot) 
                        ptc_coord = (min(ptc_pos) + 1)*3 
                        

                    allprot = allprot+prot

                exons_after = (len(s) - 2)/2 - i/2
                exons_before = i/2
                if bool_ptc and i != len(s)-2 and not (i == len(s)-4 and ptc_prot_len*3 - ptc_coord < 55): #and 55 nt rule not satisfied
                    possible_exons_before.append(exons_before)
                    possible_exons_after.append(exons_after)
                    break
        if len(possible_exons_before) != 0:
            before[junc] = median(possible_exons_before)
            after[junc] = median(possible_exons_after)
    return before, after


def long_exon_finder(failing_juncs, gene_name, transcripts_by_gene, strand, chrom, fa, exonLcutoff = 1000):
    ptc_junctions = []
    ptc_exon_lens = []
    ptc_distances = []
    long_exons = []
    for junc in failing_juncs:
        junc_added = False
        possible_transcripts = transcripts_by_gene[gene_name]
        possible_distances = []
        possible_lens = []


        for transcript in possible_transcripts:
            s = list(possible_transcripts[transcript])
            #cast to list of integers
            
            if s == []:
                continue

            new_junc = False
            failing_junc = False
            for i in range(len(s)-2, 0, -2):
                junction = tuple([s[i-1], s[i]])

                #overlapping 3' or 5' splice site, just replace the junction
                if junc[0] == junction[0] or junc[1] == junction[1]:
                    s[i-1] = junc[0]
                    s[i] = junc[1]
                    new_junc = True
                #if this does not fit in our transcript, move on
                    if not (s[i-1] > s[i-2] and s[i] < s[i+1]):
                        failing_junc = True
                    
                    break  

            #If our overlapping site configuration fails in the transcript, skip this transcript     
            if failing_junc or not new_junc:
                continue

            if strand == '-':
                s.reverse()
            allprot = Seq("")
            leftover = Seq("")

            bool_long_exon = False
            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = int(exon_coord[1])-int(exon_coord[0])
                if exlen > exonLcutoff:
                    break

                """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                we must add 1 to length"""
                endpos = (len(leftover)+exlen+1)%3

                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                        leftover = Seq("")
                    else:
                        prot = seq[:-endpos].translate()
                        leftover = seq[-endpos:]                                                                                                               


                    bool_ptc = "*" in prot 
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot)
                        ptc_coord = (min(ptc_pos) + 1)*3 

                    allprot = allprot+prot
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                    else:
                        leftover = Seq("")
                    seq = seq.reverse_complement()
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()

                    bool_ptc = "*" in prot
                    if bool_ptc:
                        ptc_pos = ptc_pos_from_prot(prot, '*')
                        ptc_prot_len = len(prot)
                        ptc_coord = (min(ptc_pos) + 1)*3 

                    allprot = allprot+prot
                
                if bool_ptc and i != len(s)-2 and not (i == len(s)-4 and ptc_prot_len*3 - ptc_coord < 55):
            
                    possible_distances.append(ptc_prot_len*3 - ptc_coord)
                    possible_lens.append(ptc_prot_len*3)
                    break

        if len(possible_distances) != 0:
            ptc_distances.append(median(multimode(possible_distances)))
            ptc_exon_lens.append(median(multimode(possible_lens)))
            ptc_junctions.append(junc)
    return ptc_junctions, ptc_distances, ptc_exon_lens

def check_utrs(junc,strand,start_codons,stop_codons):
    '''
    checks if junction is within UTR
    '''
    all_start_coords = []
    start_codons = list(start_codons)
    stop_codons = list(stop_codons)
    for i in range(len(start_codons)):
        all_start_coords.append(start_codons[i][0])
        all_start_coords.append(start_codons[i][1])
    all_stop_coords = []
    for i in range(len(stop_codons)):
        all_stop_coords.append(stop_codons[i][0])
        all_stop_coords.append(stop_codons[i][1])
    if strand == '+':
        start = max(all_start_coords)
        end = min(all_stop_coords)
    else:
        start = max(all_stop_coords)
        end = min(all_start_coords)
    min_junc = min(junc)
    max_junc = max(junc)
    if min_junc < start or max_junc > end:
        return True
    return False

def solve_NMD(chrom, strand, junc, start_codons, stop_codons, gene_name, fa,
              verbose=False, exonLcutoff=2000, max_paths=100000, max_exon_search=100000):
    '''
    Compute whether there is a possible combination that uses the junction without
    inducing a PTC. We start with all annotated stop codon and go backwards.
    '''

    seed = []

    junc.sort()
    if strand == "-":
        junc.reverse()
        
    """Quinn Comment: Adds all 'start codons' to a nested list called seed""" 
    ##in an individual transcript    
    for c in start_codons:
        if strand == "-":
            seed.append([c[1]])
        else:
            seed.append([c[0]])

    # seed starts with just stop codon and then a possible 3'ss-5'ss junction
    # without introducing a PTC [stop_codon,3'ss, 5'ss, 3'ss, ..., start_codon]
    
    junc_pass = {}
    junc_fail = {}
    path_pass = []
    
    dic_terminus = {}

    depth = 0
    paths_explored = 0
    truncated = False

    """Quinn Comment: while our seed length is greater than 0 - which means we have charted all possible paths through
    all junctions ending in a stop codon (or there is an exon longer than 1000 bp and we have no complete paths)"""
    while len(seed) > 0:
        new_seed = []
        final_check = []
        depth += 1

        # Check path limit to prevent memory blowup
        paths_explored += len(seed)
        if paths_explored > max_paths:
            logger.warning(f"{gene_name} exceeded {max_paths} paths at depth {depth}, truncating search")
            truncated = True
            break

        if verbose:
            logger.debug(f"Depth {depth}, Seed L = {len(seed)}")
                    
        for s in seed:
            # first check that the seed paths are good        
            bool_ptc = False
            leftover = ''
            if len(s) > 0:                
                leftover = Seq("")
                allprot = Seq("")

                """Quinn Comment: loop through the exons, calculating lengths"""
                for i in range(0, len(s)-1, 2):
                    exon_coord = s[i:i+2]
                    exon_coord.sort()
                    exon_coord = tuple(exon_coord)
                    exlen = exon_coord[1]-exon_coord[0]


                    """Quinn Comment: find start position relative to named start of this exon and translate to protein"""
                    """Quinn Comment: Coordinates from PERIND file and GTF file are exon start and end coordinates, so 
                    we must add 1 to length"""
                    endpos = (len(leftover)+exlen+1)%3
                    if strand == '+':
                        seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                        if endpos == 0:
                            prot = seq.translate()
                            leftover = Seq("")
                        else:
                            prot = seq[:-endpos].translate()
                            leftover = seq[-endpos:]                                                                                                               
                        allprot = allprot+prot  
                    else:
                        seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                        if endpos > 0:
                            leftover = seq[:endpos]
                            seq = seq.reverse_complement()
                            prot = seq[:-endpos].translate()
                        else:
                            leftover = Seq("")
                            seq = seq.reverse_complement()
                            prot = seq.translate()
                        allprot = allprot+prot

                    #found a PTC in this transcript if any element but the last is a stop codon    
                    bool_ptc = "*" in allprot

            """Quinn Comment: if we found a PTC, add all intron coordinate pairs involved in the transcript to junc_fail"""        
            if bool_ptc:
                #This transcript failed
                for i in range(1, len(s)-1, 2):                                                                                                                  
                    j_coord = s[i:i+2]                                                                                                                           
                    j_coord.sort()                                                                                                                             
                    j_coord = tuple(j_coord)                                                                                                                     
                    if j_coord not in junc_fail:                                                                                                                 
                        junc_fail[j_coord] = 0                                                                                                                   
                    junc_fail[j_coord] += 1  

                continue
        
            # passed

            if len(s) > 2:
                terminus = (s[-2],s[-1],leftover)
                
                if terminus in dic_terminus:
                    dic_terminus[terminus].append(tuple(s))
                    continue
                else:
                    dic_terminus[terminus] = [tuple(s)]
            
            last_pos = s[-1]
            
        #fix for abnormally long exons Quinn's addition 4/3/2025
            next_exon_max = 0

            ptc_ahead = False
            while ptc_ahead == False and next_exon_max < max_exon_search:
                next_exon_max += 500
                endpos = (len(leftover)+next_exon_max+1)%3
                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (last_pos, last_pos+next_exon_max)))
                    if endpos == 0:
                        prot = seq.translate()
                        
                    else:
                        prot = seq[:-endpos].translate()
                                                                                                                                    
                else:
                    seq = Seq(fa.fetch(chrom, (last_pos-next_exon_max,last_pos)))+leftover
                    if endpos > 0:
                        seq = seq.reverse_complement()
                        prot = seq[:-endpos].translate()
                    else:
                        seq = seq.reverse_complement()
                        prot = seq.translate()
            
                ptc_ahead = "*" in prot
            
            """Quinn Comment: check the last position of our seed to see if it is close to a start codon, within a potential exon's length,
             and add your seed plus this start codon to final_check """
            for stop in stop_codons:                
                if strand == "+" and last_pos < stop[0] and stop[0]-last_pos < next_exon_max:
                    final_check.append(s+[stop[1]])
                elif strand == "-" and last_pos > stop[1] and last_pos-stop[1] < next_exon_max:
                    final_check.append(s+[stop[0]]) 

            """Quinn Comment: add all possible places to go from our last_pos to the seed (nested list)"""
            for j0,j1 in junc:                
                if strand == "+" and last_pos < j0 and j0-last_pos < next_exon_max:
                    new_seed.append(s+[j0,j1])
                if strand == "-" and last_pos > j1 and last_pos-j1 < next_exon_max: 
                    new_seed.append(s+[j1,j0])
                    
        """Quinn Comment: Exited from s in seed loop, now we check our final_checks of the full paths, we do not
        eliminate paths based on presence of a PTC, rather we classify full complete paths without PTCs if they exist"""
        # check that the possible final paths are good
        for s in final_check:
            leftover = Seq("")
            allprot = Seq("")
            for i in range(0, len(s)-1, 2):
                exon_coord = s[i:i+2]
                exon_coord.sort()
                exon_coord = tuple(exon_coord)
                exlen = exon_coord[1]-exon_coord[0]
                endpos = (len(leftover)+exlen+1)%3
                if strand == '+':
                    seq = leftover + Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))
                    if endpos == 0:
                        prot = seq.translate()
                    else:
                        prot = seq[:-endpos].translate()
                    leftover = seq[-endpos:]                                                                                                               
                    allprot = allprot+prot  
                else:
                    seq = Seq(fa.fetch(chrom, (exon_coord[0],exon_coord[1])))+leftover
                    if endpos > 0:
                        leftover = seq[:endpos]
                        seq = seq.reverse_complement()
                        prot = seq[:-endpos].translate()
                    else:
                        leftover = Seq("")
                        seq = seq.reverse_complement()
                        prot = seq.translate()
                    allprot = allprot+prot

                #found a PTC in this transcript if any element but the last is a stop codon    
                bool_ptc = "*" in allprot[:-1]
        
            if not bool_ptc:
                # all pass
                path_pass.append(tuple(s))
                for i in range(1, len(s), 2):
                    j_coord = s[i:i+2]
                    j_coord.sort()
                    j_coord = tuple(j_coord)
                    if j_coord not in junc_pass:
                        junc_pass[j_coord] = 0
                    junc_pass[j_coord] += 1

        # Limit new_seed size to prevent memory blowup in single iteration
        if len(new_seed) > max_paths:
            logger.warning(f"{gene_name} new_seed size ({len(new_seed)}) exceeds limit, keeping first {max_paths} paths")
            new_seed = new_seed[:max_paths]
            truncated = True

        seed = new_seed
    
    
    """Quinn Comment: OUT OF WHILE LOOP through all possible paths/seeds; 
    check all terminus' to see if they are part of a full path that has been classified as passing"""
    terminus_iterations = 0
    max_terminus_iterations = 50000

    while True:
        terminus_iterations += 1
        if terminus_iterations > max_terminus_iterations:
            logger.warning(f"{gene_name} terminus checking exceeded {max_terminus_iterations} iterations, stopping")
            break

        new_paths = []
        for terminus in dic_terminus:
            terminus_pass = False
            for path_subset in dic_terminus[terminus]:
                for path in path_pass:
                    if path[:len(path_subset)] == path_subset:
                        terminus_pass = True
                        break
            #print(terminus, terminus_pass)

            """Quinn Comment: if our terminus is part of a passing path, we want to make sure if is reflected in passing paths and
            add the associate junctions to junc_pass, only if they are not present"""
            if terminus_pass:
                for path_subset in dic_terminus[terminus]:
                    if path_subset in path_pass: continue
                    new_paths.append(path_subset)
                    path_pass.append(path_subset)
                    for i in range(1, len(path_subset), 2):
                        j_coord = list(path_subset[i:i+2])
                        j_coord.sort()
                        j_coord = tuple(j_coord)
                        if j_coord not in junc_pass:
                            junc_pass[j_coord] = 0
                            if verbose:
                                logger.debug(f"junction {j_coord} pass")
        """Quinn Comment: we could have a new path_pass added, so our while loop checks again to see if there are any new paths 
        that are now going to be passing considering our additions"""
        if len(new_paths) == 0:
            break
            
    return junc_pass,junc_fail

def parse_gtf(gtf: str, 
              gene_type: str, 
              transcript_type: str,
              gene_name: str, 
              transcript_name: str
             ):
    '''Lower level function to parse GTF file
    - gtf: str : path to GTF annotation file
    - returns: dictionary with keys: 
        chrom, source, type, start, end, strand, frame, info, gene_name, transcript_type, transcript_name, gene_type
    '''
    fields = ["chrom", "source", "type", "start", "end", ".","strand","frame", "info"]
    open_gtf = lambda x: gzip.open(x) if ".gz" in x else open(x)
    for ln in open_gtf(gtf):
        ln = ln.decode('ascii') if ".gz" in gtf else ln
        dic = {}
        if ln[0] == "#": continue
        ln = ln.strip().split('\t')
        for i in range(len(fields)):
            dic[fields[i]] = ln[i]

        # Parse all attributes first
        info_fields = [{x.split()[0]: x.split()[1].replace('"', '')} 
                    for x in dic['info'].split(';') if len(x.split()) > 1]
        info_fields = {k: v for d in info_fields for k, v in d.items()}

        # Map custom attribute names to standard keys
        attribute_mapping = {
            'gene_name': gene_name,
            'transcript_name': transcript_name, 
            'gene_type': gene_type,
            'transcript_type': transcript_type
        }

        for standard_key, custom_attr in attribute_mapping.items():
            try:
                dic[standard_key] = info_fields[custom_attr]
            except KeyError:
                dic[standard_key] = None

        yield dic
         

def parse_annotation(gtf_annot: str,
                     gene_type: str, 
                     transcript_type: str,
                     gene_name: str, 
                     transcript_name: str):
    '''
    Used `parse_gtf` to first parse a gtf file, then extract and return further
    information, including gene coordinates, intron info, and splice site to 
    gene-name dictionary.
    
    - gtf_anno: str : path to GTF annotation file
    - returns: 
        - genes_coords: gene coordinates, grouped by chromosome and strand 
        - introns_info: a dictionary with these keys: junctions (ie all introns),
                        start_codon, stop_codon, utrs, pcjunctions (ie only protein coding introns) 
        - ss2gene: a dictionary with splice site as keys, and gene name as values, eg. ('chr1', 11869): 'DDX11L1'
    '''
    genes_info = {}
    introns_info = {}
    genes_coords = {}
    ss2gene = {}

    for dic in parse_gtf(gtf_annot, gene_type, transcript_type, gene_name, transcript_name):
        chrom = dic['chrom']
        gname = dic['gene_name']
        tname = dic['transcript_name'], gname
        anntype  = dic['type']
        if dic['type'] == 'gene': 
            dic['transcript_type'] = "gene"
        if dic['transcript_type'] == "nonsense_mediated_decay" and anntype == "stop_codon": 
            continue 
        if dic['transcript_type'] != "protein_coding" and anntype == "UTR": 
            continue
        """Quinn: CONVERT TO BED FORMAT IS ERROR"""
        start, end = int(dic['start']), int(dic['end']) # GTF format coordinates
        strand = dic['strand']
    
        if (chrom, strand) not in genes_coords:
            genes_coords[(chrom, strand)] = []

        if tname not in genes_info: # tname is (transcript_name, gene_name)
            genes_info[tname] = {'exons':[],
                                 'start_codon':set(), 
                                 'stop_codon':set(), 
                                 'utrs':set(),
                                 'type':dic['transcript_type']
                                 }
        
        if anntype in ["start_codon", "stop_codon"]:
            genes_info[tname][anntype].add((start,end))
        elif anntype in ["gene"]:
            genes_coords[(chrom,strand)].append(((start,end), gname))
        elif anntype in ['exon']:
            genes_info[tname]['exons'].append((start,end))
            # Store gene info for splice sites
            ss2gene[(chrom, int(dic['start']))] = dic['gene_name'] # GTF format
            ss2gene[(chrom, int(dic['end']))] = dic['gene_name'] # GTF format

        elif anntype in ['UTR']:
            genes_info[tname]['utrs'].add((start,end))

    for gene in genes_info: # this is actually transcript level!!! (transcript_name, gene_name)
        gene_name = gene[1]
        exons = genes_info[gene]['exons']
        exons.sort()

        pc = False # Basically use transcript_type == protein_coding to set flag
        if genes_info[gene]['type'] == "protein_coding": # AGAIN here gene = (transcript_name, gene_name)
            pc = True

        junctions = set() # all junctions/introns
        pcjunctions = set() # only protein coding junctions/introns

        for i in range(len(exons)-1):
            intron = (exons[i][1], exons[i+1][0])
            junctions.add(intron)
            if pc:
                pcjunctions.add(intron)

        if gene_name not in introns_info: # here only gene_name, does not incl. transcript_name
            introns_info[gene_name] = {'junctions':set(),
                                       'start_codon':set(),
                                       'stop_codon':set(),
                                       'utrs':set(),
                                       'pcjunctions':set(),
                                 }

        introns_info[gene_name]['start_codon'] = introns_info[gene_name]['start_codon'].union(genes_info[gene]['start_codon'])
        introns_info[gene_name]['stop_codon'] = introns_info[gene_name]['stop_codon'].union(genes_info[gene]['stop_codon'])
        introns_info[gene_name]['junctions'] = introns_info[gene_name]['junctions'].union(junctions)
        introns_info[gene_name]['utrs'] = introns_info[gene_name]['utrs'].union(genes_info[gene]['utrs'])
        introns_info[gene_name]['pcjunctions'] = introns_info[gene_name]['pcjunctions'].union(pcjunctions)

    return genes_coords, introns_info, ss2gene





def get_feature(fname, feature = "exon"):
    ss2gene = {}
    if ".gz" in fname:
        F = gzip.open(fname)
    else:
        F = open(fname)
    for ln in F:
        if ln[0] == "#": continue
        ln = ln.split('\t')
        gID = ln[-1].split('gene_name "')[1].split('"')[0]

        if ln[2] != feature:
            continue
        ss2gene[(ln[0], int(ln[3]))] = gID
        ss2gene[(ln[0], int(ln[4]))] = gID
    return ss2gene


def get_overlap_stream(L1, L2, relax = 0):
    '''                                                                                                                                                                
    L1 and L2 are sorted lists of tuples with key, values                                                                                                              
    '''
    i, j = 0, 0
    while i < len(L1) and j < len(L2):
        if L1[i][0][1] < L2[j][0][0]:
            i += 1
            continue
        elif L2[j][0][1] < L1[i][0][0]:
            j += 1
            continue
        else:
            k = 0
            # hits overlapping, check all L2 that may overlap with intron                                                                                              
            while L2[j+k][0][0] <= L1[i][0][1]:
                if overlaps(L1[i][0], L2[j+k][0]):
                    yield L1[i], L2[j+k]
                k += 1
                if j+k == len(L2): break
            i += 1

def overlaps(A: tuple, B: tuple):
    '''                                                                                                                                                                
    Checks if A and B overlaps                                                                                                                                         
    A: tuple : (start, end)
    B: tuple : (start, end)
    '''
    if A[1] < B[0] or B[1] < A[0]:
        return False
    else: return True


def ClassifySpliceJunction(
    perind_file: str,
    gtf_annot: str,
    fa,
    rundir: str = ".",
    outprefix: str = "Leaf2",
    max_juncs: int = 10000,
    max_paths: int = 100000,
    keepannot: bool = False,
    gene_type: str = "gene_type",
    transcript_type: str = "transcript_type",
    gene_name: str = "gene_name",
    transcript_name: str = "transcript_name",
    verbose: bool = False):
    """
    perind_file: str : LeafCutter perind counts file, e.g. leafcutter_perind.counts.gz
    gtf_annot: str : Annotation GTF file, for example gencode.v37.annotation.gtf.gz
    fa: pyfastx object
    genome: str : Reference genome fasta file
    rundir: str : run directory, default is current directory
    outprefix: str : output prefix (default: Leaf2)
    max_juncs: int : If a gene has more than `max_juncs`, skip solv_NMD function,
                     Set Coding = False (default: 10000)
    verbose: bool : verbose mode (default: False)

    NOTE (Aug 2024): This script expects all coordinates to be 1-based. GTF
    annotations and processed intron junctions should be 1-based going into
    this function and the solve_NMD function.
    """

    # read leafcutter perind file and store junctions in dictionary: dic_junc
    # key = (chrom,strand), value = list of junctions [(start,end)]
    dic_junc = {}
    logger.info(f"Processing junction counts {perind_file}...")
    for ln in gzip.open(perind_file):
        junc_info = ln.decode('ascii').split()[0] # first column
        if junc_info == "chrom": continue # header
        
        chrom, start, end, clu_strand = junc_info.split(":")
        strand = clu_strand.split("_")[-1]
        if (chrom,strand) not in dic_junc: 
            dic_junc[(chrom,strand)] = []
        dic_junc[(chrom,strand)].append((int(start), int(end)))

    #NOTE: convert start, end coordinates in dic_junc to (leafcutter1) 1-based from 0 based (like BED)
    for chrom, strand in dic_junc:
        dic_junc[(chrom, strand)] = [(x[0] , x[1] + 1) for x in dic_junc[(chrom, strand)]]

    logger.info("Processed junctions.")
    if verbose:
        for chrstrand in dic_junc:
            logger.debug(f"{len(dic_junc[chrstrand])} jxns on {chrstrand[0]} ({chrstrand[1]}).")


    # Compute content hash for stable pickle cache key (avoids basename collisions)
    with open(gtf_annot, 'rb') as _f:
        _gtf_hash = hashlib.md5(_f.read()).hexdigest()[:12]

    # load or parse gtf annotations
    # g_coords: gene coordinates, grouped by chromosome and strand
    # g_info: a dictionary with (transcript_name, gene_name) as keys, and intron info as values
    try:
        logger.info("Loading annotations...")
        parsed_gtf = f"{rundir}/genesDict_{_gtf_hash}_SJC_annotations.pckle"
        with open(parsed_gtf, 'rb') as f:
            g_coords, g_info = pickle.load(f)
        logger.info("Loaded annotations.")
    except:
        logger.info("Parsing annotations for the first time...")
        g_coords, g_info, ss2gene = parse_annotation(gtf_annot, gene_type, transcript_type, gene_name, transcript_name)
        
    
        for chrom,strand in g_coords:
            to_remove_gcoords = set()
            to_remove_ginfo = set()
            for gene in g_coords[(chrom,strand)]:
                if gene[1] in g_info: # check if gene_name is in introns_info keys
                    if len(g_info[gene[1]]['stop_codon']) == 0: # if there are no stop codons
                        to_remove_gcoords.add(gene) # mark gene for remove
                        to_remove_ginfo.add(gene[1]) # mark gene_name for removal

            for g in to_remove_gcoords:
                g_coords[(chrom,strand)].remove(g)
            for g in to_remove_ginfo:
                g_info.pop(g)
        logger.info("Saving parsed annotations...")
        if keepannot:
            with open(parsed_gtf, 'wb') as f:
                pickle.dump((g_coords, g_info), f)


    txn2gene = f"{rundir}/txn2gene_{_gtf_hash}_SJC_annotations.pckle"
    try:
        logger.info("Loading txn2gene annotations...")
        with open(txn2gene, 'rb') as f:
            transcripts_by_gene = pickle.load(f)
    except:
        logger.info("Failed... Making txn2gene annotations...")
        transcripts_by_gene = tx_by_gene(gtf_annot, gene_type, transcript_type, gene_name, transcript_name)
        if keepannot:
            with open(txn2gene, 'wb') as f:
                pickle.dump(transcripts_by_gene, f)
    logger.info("Loaded txn2gene annotations.")

    nmd_tx2gene = f"{rundir}/nmd_txn2gene_{_gtf_hash}_SJC_annotations.pckle"
    try:
        logger.info("Loading nmd_txn2gene annotations...")
        with open(nmd_tx2gene, 'rb') as f:
            nmd_tx_by_gene = pickle.load(f)
    except:
        logger.info("Failed... Making nmd_txn2gene annotations...")
        nmd_tx_by_gene = NMD_tx(gtf_annot, gene_type, transcript_type, gene_name, transcript_name)
        if keepannot:
            with open(nmd_tx2gene, 'wb') as f:
                pickle.dump(nmd_tx_by_gene, f)
    logger.info("Loaded nmd_txn2gene annotations.")

    gene_juncs = {}
    for chrom,strand in dic_junc:
        if (chrom,strand) not in g_coords: 
            logger.error(f"Could not find {chrom} ({strand}) in annotations...")
            continue
        juncs = [(x,x) for x in dic_junc[(chrom,strand)]]
        juncs.sort()

        coords = g_coords[(chrom,strand)]
        coords.sort()
        
        # save junctions that overlapping a gene in gene_juncs dictionary: (gene_name, chrom, strand) : [junctions]
        for junc, geneinfo in get_overlap_stream(juncs,coords): 
            info = (geneinfo[1], chrom,strand)
            if info not in gene_juncs:
                gene_juncs[info] = []
            gene_juncs[info].append(junc[0])

    fout = open(f"{rundir}/clustering/{outprefix}_junction_classifications.txt",'w')
    fout.write("\t".join(["Gene_name","Intron_coord","Strand","Annot","Coding","UTR","GencodePC"])+'\n')
    lout = open(f"{rundir}/clustering/{outprefix}_long_exon_distances.txt",'w')
    lout.write("\t".join(["Gene_name","Intron_coord","PTC_position","Exon_length"])+'\n')
    eout = open(f"{rundir}/clustering/{outprefix}_exon_stats.txt",'w')
    eout.write("\t".join(["Gene_name","Intron_coord","Exons_before","Exons_after"])+'\n')
    nout = open(f"{rundir}/clustering/{outprefix}_nuc_rule_distances.txt",'w')
    nout.write("\t".join(["Gene_name","Intron_coord","ejc_distance"])+'\n')

    for gene_name, chrom, strand in gene_juncs:

        extra_utr_rule = {'junction': [], 'annotated':[], 'coding': [], 'utr': [], 'gencode': []}

        logger.debug(f"Processing {gene_name} ({chrom}:{strand})")
        
        query_juncs = gene_juncs[(gene_name,chrom,strand)] # from LeafCutter perind file
        if gene_name not in g_info: continue
        junctions = g_info[gene_name]['junctions'] # from annotation
        
        # classify all junctions in gene
        junctions = list(junctions.union(query_juncs))

        start_codons = g_info[gene_name]['start_codon'] 
        stop_codons = g_info[gene_name]['stop_codon']

        if verbose:
            logger.debug(f"LeafCutter junctions ({len(query_juncs)}) All junctions ({len(junctions)}) Start codons ({len(start_codons)}) Stop codons ({len(stop_codons)})")

        
        if len(junctions) <= max_juncs:
            try:
                junc_pass, junc_fail = solve_NMD(chrom, strand, junctions,
                                                start_codons, stop_codons,
                                                gene_name, fa,
                                                max_paths=max_paths)
                junc_fail = set(junc_fail.keys())
                junc_pass = set(junc_pass.keys())
            except:
                if verbose:
                    logger.warning("Skipping... seq fetching problem")
                continue
        else:
            if verbose:
                logger.info(f"Skipping... Too many juncs (> {str(max_juncs)})")
            junc_pass = set()
            #Don't want to run all these junctions through the extra rules
            junc_fail = set()


        failing_juncs = junc_fail.difference(junc_pass)


        old_junc_pass = junc_pass
        junc_pass = {}
        junc_pass['normal'] = old_junc_pass
        ptc_junctions, ptc_distances, ptc_exon_lens = long_exon_finder(failing_juncs, gene_name, transcripts_by_gene, strand, chrom, fa)
        exons_before, exons_after = many_junctions(failing_juncs, gene_name, transcripts_by_gene, strand, chrom, fa)
        junc_pass['nuc_rule'], ejc_distances = nucleotide_rule(failing_juncs, gene_name, transcripts_by_gene, strand, chrom, nmd_tx_by_gene, fa)
        for j in junctions:

            bool_pass = j in junc_pass['normal'] 
            gencode = j in g_info[gene_name]['pcjunctions']
            bool_fail = j in failing_juncs

            if bool_fail or bool_pass:
                tested = True
            else:
                tested = False
            annotated = j in g_info[gene_name]['junctions']
            
            if len(start_codons) != 0 and len(stop_codons) != 0:
                utr = check_utrs(j, strand, start_codons, stop_codons)
            else:
                utr = False

            #if not bool_pass and annotated:
            #print("%s %s %s junction: %s tested: %s utr: %s coding: %s annotated: %s "%(chrom, strand, gene_name, j, tested,utr, bool_pass, annotated))
            
            extra_utr_rule['junction'].append(j)
            extra_utr_rule['annotated'].append(annotated)
            extra_utr_rule['coding'].append(bool_pass)
            extra_utr_rule['utr'].append(utr)
            extra_utr_rule['gencode'].append(gencode)
            
            #NOTE: revert to leafcutter2 bed format coordiantes for junctions
        coding = [extra_utr_rule['junction'][k] for k in [i for i,bol in enumerate(extra_utr_rule['coding']) if bol]]
        if strand == '+':
            coding_5_prime = [j[0] for j in coding]
        else:
            coding_5_prime = [j[1] for j in coding]

        utr_js = []
        for i in range(len(extra_utr_rule['annotated'])):
            j = extra_utr_rule['junction'][i]
            utr = extra_utr_rule['utr'][i]
            if strand == '+':
                if utr and j[0] in coding_5_prime and j[0] > max([k[0] for k in start_codons]):
                    utr = False
            else:
                if utr and j[1] in coding_5_prime and j[1] < min([k[1] for k in start_codons]):
                    utr = False
            annotated = extra_utr_rule['annotated'][i]
            bool_pass = extra_utr_rule['coding'][i]
            gencode = extra_utr_rule['gencode'][i]
            fout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]-1}',strand,
                                    str(annotated), str(bool_pass), str(utr), str(gencode)])+'\n')
            if utr:
                utr_js.append(j)
        
        for w in range(len(ptc_junctions)):
            j = ptc_junctions[w]
            if j not in utr_js:
                if j not in g_info[gene_name]['pcjunctions']:
                    lout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]-1}',
                                        str(ptc_distances[w]), str(ptc_exon_lens[w])])+'\n')
        for j in exons_before:
            if j not in utr_js:
                if j not in g_info[gene_name]['pcjunctions']:
                    eout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]-1}',
                                        str(exons_before[j]), str(exons_after[j])])+'\n')

        for w in range(len(junc_pass['nuc_rule'])):
            j = junc_pass['nuc_rule'][w]
            if j not in utr_js:
                if j not in g_info[gene_name]['pcjunctions']:
                    nout.write('\t'.join([gene_name, f'{chrom}:{j[0]}-{j[1]-1}',
                                        str(ejc_distances[w])])+'\n')

