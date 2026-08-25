#!/usr/bin/env python

import sys
import bedparse
import argparse
import gzip
import re
import pyfastx
import subprocess
from io import StringIO
from collections import defaultdict
import shutil
import pandas as pd
import csv
import logging
import tempfile
import os
from Bio import bgzf
from Bio import motifs
from Bio.Seq import Seq
import random

def reorder_gtf(gtf_stringio, output_gtf, mode='w'):
    """Reorder GTF lines by gene and transcript, ensuring correct hierarchical sorting,
    and write to output file."""
    # Read the GTF data into a pandas DataFrame
    columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
    # Explicit dtypes to avoid DtypeWarning (score/frame can be '.' or numeric)
    dtype_map = {
        'seqname': 'string',
        'source': 'string',
        'feature': 'string',
        'start': 'int64',
        'end': 'int64',
        'score': 'string',   # keep as string to handle '.'
        'strand': 'string',
        'frame': 'string',   # keep as string to handle '.'
        'attribute': 'string'
    }
    df = pd.read_csv(
        gtf_stringio,
        sep='\t',
        names=columns,
        comment='#',
        dtype=dtype_map
    )
    # Extract gene_id and transcript_id from the attribute column
    df['gene_id'] = df['attribute'].str.extract(r'gene_id "([^"]+)"')
    df['transcript_id'] = df['attribute'].str.extract(r'transcript_id "([^"]+)"')

    # Determine if the line is a gene line
    df['is_gene'] = df['transcript_id'].isnull()

    # Create a column for the start position of the parent gene
    df['start_parent_gene'] = df.groupby('gene_id')['start'].transform('min')

    # Create a column for the start position of the parent transcript
    df['start_parent_transcript'] = df.groupby('transcript_id')['start'].transform('min')

    # Fill NaNs in start_parent_transcript with start_parent_gene for gene lines
    df['start_parent_transcript'] = df['start_parent_transcript'].fillna(df['start_parent_gene'])

    df = df.copy()
    df.loc[:, 'seqname'] = df['seqname'].astype(str)

    # Sort by seqname, start_parent_gene, gene_id, is_gene, start_parent_transcript, transcript_id, start
    df_sorted = df.sort_values(by=['seqname', 'start_parent_gene', 'gene_id','strand', 'is_gene', 'start_parent_transcript', 'transcript_id', 'start'], ascending=[True, True, True, True, False, True, True, True]).reset_index(drop=True)

    # Drop the helper columns before writing to the output
    df_sorted = df_sorted.drop(columns=['gene_id', 'transcript_id', 'is_gene', 'start_parent_gene', 'start_parent_transcript'])

    # Write the sorted DataFrame to the output GTF file
    df_sorted.to_csv(output_gtf, sep='\t', header=False, index=False, mode=mode, quoting=csv.QUOTE_NONE)

def run_bedparse_gtf2bed(gtf_file, *args, n=None):
    """
    Wrapper around command line `bedparse gtf2bed`.
    Ensures input is plain-text GTF (decompresses .gz to a temp file if needed).
    Optionally samples first n lines.
    """
    # Ensure we have a plain-text GTF path (decompress if .gz)
    cleanup_paths = []
    if gtf_file.endswith('.gz'):
        tmp_plain = tempfile.NamedTemporaryFile(delete=False, suffix='.gtf')
        cleanup_paths.append(tmp_plain.name)
        with gzip.open(gtf_file, 'rt') as inf, open(tmp_plain.name, 'w') as outf:
            for line in inf:
                outf.write(line)
        gtf_plain = tmp_plain.name
    else:
        gtf_plain = gtf_file

    # If sampling n lines, create a sampled temp file from plain GTF
    if n is not None:
        tmp_sample = tempfile.NamedTemporaryFile(delete=False, suffix='.gtf')
        cleanup_paths.append(tmp_sample.name)
        with open(gtf_plain, 'r') as f, open(tmp_sample.name, 'w') as out:
            for i, line in enumerate(f):
                if i >= n:
                    break
                out.write(line)
        gtf_input = tmp_sample.name
    else:
        gtf_input = gtf_plain

    # Construct and run bedparse command
    command = ['bedparse', 'gtf2bed', gtf_input]
    command.extend(args)
    try:
        result = subprocess.run(command, capture_output=True, text=True)
    except FileNotFoundError as exc:
        raise RuntimeError(
            "Could not execute `bedparse`. It is a required dependency and must be "
            "on PATH. Install it with `pip install bedparse` or "
            "`conda install -c bioconda bedparse`."
        ) from exc
    finally:
        # Cleanup temp files
        for p in cleanup_paths:
            try:
                os.unlink(p)
            except Exception:
                pass

    if result.returncode != 0:
        raise RuntimeError(
            f"`{' '.join(command)}` failed with exit code {result.returncode}.\n"
            f"bedparse stderr:\n{result.stderr}"
        )
    return result.stdout

def kozak_motif_from_counts(counts_dict=None, pseudocount=0.5):
    """
    Build a Kozak motif from a counts matrix (dict of lists for A,C,G,T).
    Returns (motif, pwm, pssm).
    """
    if counts_dict is None:
        counts_dict = {
            'A': [3620.0, 3166.0, 4450.0, 8952.0, 5277.0, 3288.0, 18870.0, 0.0, 2.0, 3973.0],
            'C': [4423.0, 6320.0, 7562.0, 1691.0, 7700.0, 8999.0, 25.0, 3.0, 0.0, 2652.0],
            'G': [7701.0, 5952.0, 4862.0, 7331.0, 3652.0, 5302.0, 2.0, 0.0, 18895.0, 9712.0],
            'T': [3155.0, 3461.0, 2025.0, 925.0, 2270.0, 1310.0, 2.0, 18896.0, 2.0, 2562.0],
        }

    # Validate
    lengths = {len(v) for v in counts_dict.values()}
    if len(lengths) != 1:
        raise ValueError("All nucleotide lists must have the same length")
    for b in ("A","C","G","T"):
        if b not in counts_dict:
            raise ValueError(f"Missing base '{b}' in counts_dict")

    # Use PositionSpecificCounts (Biopython >= 1.78)
    cnt = matrix.PositionSpecificCounts("ACGT", counts_dict)
    m = motifs.Motif(cnt)
    pwm = cnt.normalize(pseudocounts=pseudocount)
    pssm = pwm.log_odds()
    return m, pwm, pssm

def build_pssm_from_counts_by_sampling(counts_dict, n_instances=500, seed=1):
    """
    Build a PSSM by sampling sequences from a per-position nucleotide counts dict.
    counts_dict: {'A': [...], 'C': [...], 'G': [...], 'T': [...]}, all lists same length (e.g., 10 for Kozak)
    n_instances: number of synthetic sequences to sample
    seed: RNG seed for reproducibility
    Returns (motif, pwm, pssm)
    """
    # Validate
    if not counts_dict or any(b not in counts_dict for b in ('A','C','G','T')):
        raise ValueError("counts_dict must include A,C,G,T")
    lengths = {len(v) for v in counts_dict.values()}
    if len(lengths) != 1:
        raise ValueError("All nucleotide lists must have same length")
    L = lengths.pop()

    # Build per-position population lists for random choice
    # To keep memory modest, normalize counts to probabilities
    probs = []
    for i in range(L):
        a = counts_dict['A'][i]
        c = counts_dict['C'][i]
        g = counts_dict['G'][i]
        t = counts_dict['T'][i]
        total = a + c + g + t
        # Guard against zero total
        if total <= 0:
            probs.append((('A', 0.25), ('C', 0.25), ('G', 0.25), ('T', 0.25)))
        else:
            probs.append((
                ('A', a/total),
                ('C', c/total),
                ('G', g/total),
                ('T', t/total)
            ))

    random.seed(seed)
    seqs = []
    for _ in range(n_instances):
        s = []
        for i in range(L):
            bases, p = zip(*probs[i])
            # random.choices is available in Python 3.6+
            s.append(random.choices(bases, weights=p, k=1)[0])
        seqs.append(Seq(''.join(s)))

    m = motifs.create(seqs)
    pwm = m.counts.normalize(pseudocounts=0.5)
    pssm = pwm.log_odds()
    return m, pwm, pssm

def load_jaspar_pfm_pssm(pfm_path, pseudocount=0.5):
    """Load a JASPAR PFM using Biopython and return (pwm, pssm)."""
    with open(pfm_path) as handle:
        motif_list = list(motifs.parse(handle, "jaspar"))
    if not motif_list:
        raise ValueError("No motifs parsed from JASPAR PFM")
    m = motif_list[0]
    pwm = m.counts.normalize(pseudocounts=pseudocount)
    pssm = pwm.log_odds()
    return pwm, pssm

def extract_kozak_seq(marked_sequence, before=6, after=4):
    """
    Extract sequence surrounding '^' from `before` nt before the marker to `after` nt after the marker.
    Removes '|' and '*' before computing indices. Returns the nucleotide string (no markers).
    If '^' is not present, returns empty string.

    Example:
        marked = "AA|AA|^ATG|GGC*|TT"
        extract_kozak_seq(marked, before=6, after=1) -> "AAATG"  (depending on exact positions)
    """
    if '^' not in marked_sequence:
        return ""

    # Remove splice-junction markers and stop markers but keep the '^' so we can locate it
    cleaned = marked_sequence.replace('|', '').replace('*', '')

    caret_idx = cleaned.find('^')
    if caret_idx == -1:
        return ""

    start = max(0, caret_idx - before)
    end = min(len(cleaned), caret_idx + 1 + after)  # +1 to include the base at '^' position

    window = cleaned[start:end]
    # Remove the caret before returning the nucleotide sequence
    return window.replace('^', '')

def score_kozak_window(marked_sequence, pssm=None, before=6, after=4):
    """
    Extract window around '^' and score with provided PSSM.
    Returns a 3-decimal string; returns 'NA' on any error and logs a helpful message.
    """
    try:
        window = extract_kozak_seq(marked_sequence, before=before, after=after)
        if not window:
            logging.debug("score_kozak_window: no '^' marker found or empty window")
            return 'NA'
        if len(window) != 10:
            logging.debug(f"score_kozak_window: invalid window length {len(window)} (expected 10)")
            return 'NA'
        window = window.upper()
        if any(ch not in 'ACGT' for ch in window):
            logging.debug(f"score_kozak_window: non-ACGT character in window '{window}'")
            return 'NA'
        if pssm is None:
            logging.debug("score_kozak_window: PSSM is None (not initialized)")
            return 'NA'
        s = pssm.calculate(window)
        # Biopython returns scalar when len(window)==motif_len; else list of sliding scores
        if isinstance(s, (list, tuple)):
            if not s:
                logging.debug("score_kozak_window: PSSM.calculate returned empty scores list")
                return 'NA'
            s = s[0]
        return f"{s:.3f}"
    except Exception as e:
        logging.debug(f"score_kozak_window: exception during scoring ({e})")
        return 'NA'

def count_bars_until_n_position(sequence, n):
    bar_count = 0
    non_bar_count = 0
    for i,char in enumerate(sequence):
        if non_bar_count == n:
            break
        if char == "|":
            bar_count += 1
        else:
            non_bar_count += 1
    return bar_count

def AddORF_Marks(sequence, StartMarkBasePosition=None, StopMarkBasePosition=None):
    """
    returns sequence with '^' at StartMarkBasePosition and '*' at StopMarkBasePosition.
    '|' strings are ignored when converting base positions to string indices.
    If StopMarkBasePosition is None, only insert the start mark.
    """
    # Determine insertion indices in the original string accounting for '|' markers
    if StartMarkBasePosition is not None and StartMarkBasePosition >= 0:
        start_offset = count_bars_until_n_position(sequence, StartMarkBasePosition)
        start_idx = StartMarkBasePosition + start_offset
    else:
        start_idx = None

    if StopMarkBasePosition is not None and StopMarkBasePosition >= 0:
        stop_offset = count_bars_until_n_position(sequence, StopMarkBasePosition)
        stop_idx = StopMarkBasePosition + stop_offset
    else:
        stop_idx = None

    # Build the output in a single pass
    if start_idx is None and stop_idx is None:
        return sequence
    if start_idx is not None and (stop_idx is None or stop_idx < start_idx):
        # Insert only start, or start before stop (no overlap)
        return sequence[:start_idx] + '^' + sequence[start_idx:] if stop_idx is None else (
            sequence[:start_idx] + '^' + sequence[start_idx:stop_idx] + '*' + sequence[stop_idx:]
        )
    elif start_idx is None and stop_idx is not None:
        # Only stop mark
        return sequence[:stop_idx] + '*' + sequence[stop_idx:]
    else:
        # start_idx <= stop_idx: insert both
        return sequence[:start_idx] + '^' + sequence[start_idx:stop_idx] + '*' + sequence[stop_idx:]
AddORF_Marks("||ATG|G|ATAGG", 1, 5)

def extract_sequence(self, fasta_obj, AddMarksForORF=False):
    """
    ...method to be used on bedparse.bedline object. will monkey patch it into bedparse.bedline class
    will also add '|' between blocks to mark exon junctions
    """
    # Extract information from the BED entry
    chrom = self.chr
    start = self.start
    end = self.end
    strand = self.strand if hasattr(self, 'strand') else '+'
    block_sizes = list(map(int, self.exLengths.split(','))) if hasattr(self, 'exLengths') else [self.end - self.start]
    block_starts = list(map(int, self.exStarts.split(','))) if hasattr(self, 'exStarts') else [0]

    sequence = ''
    # Retrieve the sequence for each block and concatenate them
    try:
        if block_sizes and block_starts:
            for block_start, block_size in zip(block_starts, block_sizes):
                block_start_genomic = start + block_start
                block_end_genomic = block_start_genomic + block_size
                block_seq = fasta_obj.fetch(chrom, (block_start_genomic+1, block_end_genomic))
                sequence += "|" + block_seq
            sequence = sequence.rstrip('|').lstrip('|')
        else:
            sequence = fasta_obj.fetch(chrom, (start+1, end))
    except KeyError:
        # if contig not present in fasta
        sequence = ""
    # Reverse complement if on the negative strand
    if strand == '-':
        sequence = reverse_complement(sequence)
        # Add ORF markers if requested
    if AddMarksForORF and self.cds():
        if self.utr(which=5):
            cds_relative_start = sum([int(i) for i in self.utr(which=5).exLengths.split(',')])
        else:
            cds_relative_start = 0
        cds_relative_end = cds_relative_start + sum([int(i) for i in self.cds().exLengths.split(',')])
        sequence = AddORF_Marks(sequence, cds_relative_start, cds_relative_end)
    return sequence

def reverse_complement(seq):
    # Helper function to reverse complement a DNA sequence
    complement = str.maketrans('ATCGatcg', 'TAGCtagc')
    return seq.translate(complement)[::-1]

# Monkey-patch the BEDLine class
bedparse.bedline.extract_sequence = extract_sequence
bedparse.bedline.reverse_complement = staticmethod(reverse_complement)

def insert_marks_for_longset_ORF(sequence, require_ATG=True, require_STOP=True, min_ORF_len=0):
    """
    return sequence with "^" to mark start codon, "*" to mark stop for longest ORF that meets minimum length requirement.
    If not require_ATG, can start translation from beginning, which could be useful if transcript starts are not well defined. 
    If not require_STOP, translation does not need stop codon at end. 
    "|" characters (which I use mark splice junctions) are ignored, allowing codons to cross exon junctions.
    min_ORF_len in codons after start codon.
    """
    if require_ATG and require_STOP:
        regex = r"(?=(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)))"
    elif require_ATG and not require_STOP:
        regex = r"(?=(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)*))"
    elif (not require_ATG) and require_STOP:
        regex = r"(?=((?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)))"
    else:
        regex = r"(?=(\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])*(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)*))"
    
    all_orfs = re.findall(regex, sequence, flags=re.IGNORECASE)
    
    # Filter ORFs by minimum length and find the longest
    valid_orfs = []
    for orf_match in all_orfs:
        # Calculate ORF length without splice junction markers
        clean_orf = re.sub(r'[\^\|\*]', '', orf_match)
        orf_length_codons = len(clean_orf) // 3
        
        if orf_length_codons >= min_ORF_len:
            valid_orfs.append((orf_match, len(orf_match)))
    
    if valid_orfs:
        # Get the longest valid ORF
        longest_orf_match = max(valid_orfs, key=lambda x: x[1])[0]
        start_codon_pos = sequence.find(longest_orf_match)
        
        # Handle marking based on whether stop codon is present/required
        if require_STOP or re.search(r'(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)$', longest_orf_match, flags=re.IGNORECASE):
            # Has stop codon or stop required - mark both start and stop
            stop_codon_pos = start_codon_pos + len(longest_orf_match)
            return sequence[0:start_codon_pos] + "^" + longest_orf_match + "*" + sequence[stop_codon_pos:]
        else:
            # No stop codon and stop not required - mark start and let ORF continue to end
            return sequence[0:start_codon_pos] + "^" + sequence[start_codon_pos:]
    else:
        return sequence
insert_marks_for_longset_ORF("ggg|ATGATG|aaaCatgCC|TA|A|gggAAACCCAAACCCAAACCCAAAccc", require_ATG=True, require_STOP=True, min_ORF_len=4)

def insert_marks_for_first_ORF(sequence, require_STOP=True, min_ORF_len=0):
    """
    return sequence with "^" to mark start codon, "*" to mark stop for first ORF. If not require_STOP, translation does not need stop codon at end. "|" characters (which I use mark splice junctions) are ignored, allowing codons to cross exon junctions. min_ORF_len in codons after start codon.
    """
    # Build regex with single capture group for entire ORF
    if require_STOP:
        regex = r"^.*?(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)\|?[ACGTN]\|?[ACGTN]\|?[ACGTN]){" + str(min_ORF_len-2) + r",}?(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)).*$"
    else:
        regex = r"^.*?(\|?A\|?T\|?G\|?(?:(?!\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)(?:\|?[ACGTN]\|?[ACGTN]\|?[ACGTN])(?<!$)){" + str(min_ORF_len-2) + r",}(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A|.{3,5}$)).*$"
    
    # Find first matching ORF directly
    first_orf_match = re.match(regex, sequence, flags=re.IGNORECASE)
    
    if first_orf_match:
        orf_sequence = first_orf_match.group(1)
        start_codon_pos = first_orf_match.start(1)
        
        # Check if ORF ends with stop codon
        if require_STOP or re.search(r'(?:\|?T\|?A\|?A|\|?T\|?A\|?G|\|?T\|?G\|?A)$', orf_sequence, flags=re.IGNORECASE):
            # Has stop codon - mark both start and stop
            orf_end_pos = start_codon_pos + len(orf_sequence)
            return sequence[0:start_codon_pos] + "^" + orf_sequence + "*" + sequence[orf_end_pos:]
        else:
            # No stop codon - mark start and continue to end
            return sequence[0:start_codon_pos] + "^" + sequence[start_codon_pos:]
    
    # No valid ORF found
    return sequence
# # should mark at ORF of at least 4 codons (including start and stop)
# insert_marks_for_first_ORF("ggg|ATG|aaaCCC|TA|A|ggg", require_STOP=True, min_ORF_len=4)
# # should marks an ORF of at least 4 codons including start and stop
# insert_marks_for_first_ORF("ggg|ATG|aaaCCC|TA|A|ggg", require_STOP=False, min_ORF_len=4)
# # should does not mark the ORF of only 3 codons including start and stop
# insert_marks_for_first_ORF("ggg|ATG|aaa|TA|A|ggg", require_STOP=False, min_ORF_len=4)
# # should mark an ORF of at least 4 codons without a stop codon
# insert_marks_for_first_ORF("ggg|ATG|aaaCCC|ggg", require_STOP=False, min_ORF_len=4)
# # should not mark an ORF, when it is only 2 non-start codons.
# insert_marks_for_first_ORF("ggg|ATG|aaaCCC|", require_STOP=False, min_ORF_len=4)


def insert_marks_for_defined_ORF(sequence, start_codon_pos=None):
    """
    return sequence with "^" to mark start codon at provided string position, "*" to mark stop for longest ORF. If no stop codon is found, then no * will be inserted.
    """
    if start_codon_pos != None and start_codon_pos < len(sequence):
        # orf_match = re.search(r"^\|?(A\|?T\|?G\|?(?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])*?)\|?(?:T\|?A\|?A|T\|?A\|?G|T\|?G\|?A)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        # orf_match_nostop = re.search(r"^\|?(A\|?T\|?G\|?(?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])+)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        orf_match = re.search(r"^\|?((?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])*?)\|?(?:T\|?A\|?A|T\|?A\|?G|T\|?G\|?A)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        orf_match_nostop = re.search(r"^\|?((?:\|?[AGCTN]\|?[AGCTN]\|?[ACGTN])+)", sequence[start_codon_pos:], flags=re.IGNORECASE)
        if orf_match:
            # orf match with start and stop codon
            return sequence[:start_codon_pos] + "^" + orf_match.group(0) + "*" + sequence[orf_match.span(0)[1] + start_codon_pos:]
        elif orf_match_nostop:
            # orf match with no stop
            return sequence[:start_codon_pos] + "^" + sequence[start_codon_pos:]
        else:
            return sequence
    else:
        return sequence
# insert_marks_for_defined_ORF("GGGATGAAAGGGAAA|GGG|T|AA|GGGAAA", 3)
# insert_marks_for_defined_ORF("GGGATGAAAGGGAAA|GGG|TT|AA|GGGAAA", 3)

def find_uorfs(seq):
    junk = {"|", "^", "*"}
    stops = {"TAA", "TAG", "TGA"}

    # --- strip junk for matching, but keep map back to original ---
    clean = []
    clean_map = []
    for i, c in enumerate(seq):
        if c not in junk:
            clean.append(c.upper())
            clean_map.append(i)

    clean = "".join(clean)

    # --- main ORF start boundary in clean coordinates ---
    main_start_pos = seq.index("^")
    main_start_clean = sum(1 for x in seq[:main_start_pos] if x not in junk)

    raw_hits = []

    # --- find all start→stop ORFs before ^ ---
    for i in range(0, main_start_clean - 2):
        if clean[i:i+3] == "ATG":

            # find first in-frame stop
            stop_pos = None
            for j in range(i+3, len(clean)-2, 3):
                codon = clean[j:j+3]
                if codon in stops:
                    stop_pos = j + 3  # exclusive
                    break
            if stop_pos is None:
                continue

            start_orig = clean_map[i]
            end_orig   = clean_map[stop_pos-1]

            raw_hits.append((i, stop_pos, start_orig, end_orig))

    # ==========================================================
    #     REMOVE in-frame ORFs that are fully contained inside
    # ==========================================================

    # Group ORFs by frame
    by_frame = {0: [], 1: [], 2: []}
    for cstart, cend, ostart, oend in raw_hits:
        frame = cstart % 3
        by_frame[frame].append((cstart, cend, ostart, oend))

    final_hits = []

    for frame, hits in by_frame.items():
        # sort by cleaned-coordinate length descending
        hits.sort(key=lambda x: (x[1] - x[0]), reverse=True)

        kept = []
        for h in hits:
            cstart, cend, ostart, oend = h

            # reject if fully contained in any kept hit (same frame only)
            contained = False
            for k in kept:
                kcstart, kcend, _, _ = k
                if kcstart <= cstart and cend <= kcend:
                    contained = True
                    break

            if not contained:
                kept.append(h)

        final_hits.extend(kept)

    # return actual sequence segments
    return [seq[o_start:o_end+1] for _, _, o_start, o_end in final_hits]


def Analyze_uORFs(sequence, bedline, pssm=None):
    """
    Analyze upstream ORFs in a sequence marked with '^' for main ORF start and '*' for main ORF stop.
    Returns two genomic position fields (starts and stops) plus classifications, lengths (aa) and relative positions.

    Returns:
        tuple: (uorf_start_genomic_positions, uorf_stop_genomic_positions, classifications, lengths, relative_positions, uorf_kozak_scores)
               each element is a semicolon-separated string or empty string if none.
    """
    try:
        main_orf_index = sequence.index('^')
    except ValueError:
        return "", "", "", "", "", ""
    uorf_seqs = find_uorfs(sequence)
    if not uorf_seqs:
        return "", "", "", "", "", ""
    clean_main_start = len(re.sub(r'[\^\|\*]', '', sequence[:main_orf_index]))
    starts, stops, classes, lengths, relpos, scores = [], [], [], [], [], []
    for uorf in uorf_seqs:
        uorf_start_marked = sequence.find(uorf)
        uorf_end_marked = uorf_start_marked + len(uorf)
        classification = "uORF" if uorf_end_marked <= main_orf_index else "Overlapping_uORF"
        clean_start = len(re.sub(r'[\^\|\*]', '', sequence[:uorf_start_marked]))
        clean_end = len(re.sub(r'[\^\|\*]', '', sequence[:uorf_end_marked])) - 1
        try:
            gstart = get_absolute_pos(bedline, clean_start)
        except Exception:
            gstart = None
        try:
            gstop = get_absolute_pos(bedline, clean_end)
        except Exception:
            gstop = None
        relative_position = clean_start - clean_main_start
        aa_len = len(re.sub(r'[\^\|\*]', '', uorf)) // 3
        starts.append(str(int(gstart)) if gstart is not None else "")
        stops.append(str(int(gstop)) if gstop is not None else "")
        classes.append(classification)
        lengths.append(str(int(aa_len)))
        relpos.append(str(int(relative_position)))
        # build a temporary sequence with '^' at uORF start for scoring
        try:
            temp_seq = AddORF_Marks(sequence.replace('^', ''), StartMarkBasePosition=clean_start, StopMarkBasePosition=None)
            scores.append(score_kozak_window(temp_seq, pssm=pssm))
        except Exception:
            scores.append('NA')
    return (
        ";".join(starts) if any(starts) else "",
        ";".join(stops) if any(stops) else "",
        ";".join(classes) if classes else "",
        ";".join(lengths) if lengths else "",
        ";".join(relpos) if relpos else "",
        ";".join(scores) if scores else ""
    )

def Is_bedline_complete(bedline):
    """
    after run_bedparse_gtf2bed, if a transcript feature is read but only some of the child exons, then transcript feature spans longer than all the child exons and the bedline object is buggy.
    """
    return int(bedline.exStarts.split(',')[-1] ) + int(bedline.exLengths.split(',')[-1] ) + bedline.start == bedline.end


def get_NMD_detective_B_classification(sequence):
    """
    sequence should be marked with '^' for start, '*' for stop, and '|' for splice juncs
    """
    CDS = re.search(r"\^(\w+)\**", sequence.replace("|", ""))
    InternalStopExon = re.search(r"(^|\|)([\^ACGTNacgtn]*\*[ACGTNacgtn]*)\|", sequence)
    if "^" not in sequence or CDS == None:
        return "No CDS"
    elif re.search(r"\^[\w|]+$", sequence):
        return "No stop"
    elif re.search(r"(^|\|)[\^ACGTNacgtn]*\*[ACGTNacgtn]*$", sequence):
        return "Last exon"
    elif len(CDS.group(1)) <= 125:
        return "Start proximal"
    elif len(InternalStopExon.group(2)) >= 407:
        return "Long exon"
    elif re.search(r"\*[ACGTNacgtn]{0,50}\|[ACGTNacgtn]+$", sequence):
        return "50 nt rule"
    else:
        return "Trigger NMD"
# get_NMD_detective_B_classification("ACGTACG|CACGT")
# get_NMD_detective_B_classification("A^ATGACG|CACGT")
# get_NMD_detective_B_classification("A^CGTACG|CAC*GT")
# get_NMD_detective_B_classification("A^CGTA*CG|CACGT")
# get_NMD_detective_B_classification("ACG^TACG|" + "A" * 407 +"A*A|CACGT")
# get_NMD_detective_B_classification("ACG^TA" + "A" * 100 + "CG|" + "A" * 40 +"A*A|CACGT")
# get_NMD_detective_B_classification("ACG^TA" + "A" * 150 + "CG|A*" + "A"*60 + "|CACGT")

def get_NMD_detective_B_classification_number(NMD_detective_B_classification):
    OrdinalDict = {"Last exon":1,"Start proximal":2,"50 nt rule":3, "Long exon":4,"Trigger NMD":5, "No stop":6, "No CDS":7}
    try: 
        return OrdinalDict[NMD_detective_B_classification]
    except KeyError:
        return 8

def get_NMD_detective_B_classification_color(NMD_detective_B_classification):
    OrdinalDict = {"Last exon":"#08519c","Start proximal":"#6baed6","50 nt rule":"#c6dbef", "Long exon":"#fcbba1","Trigger NMD":"#de2d26", "No stop":"#a50f15", "No CDS":"#969696"}
    try: 
        return OrdinalDict[NMD_detective_B_classification]
    except KeyError:
        return "#252525"

def calculate_frames(bedline):
    """Calculate the frame for each CDS block."""
    CDS_bed12 = bedline.cds()
    if CDS_bed12:
        frames = []
        current_frame = 0
        bed6_blocks = list(CDS_bed12.bed12tobed6())
        if bedline.strand == '+':
            for cds in bed6_blocks:
                frames.append(current_frame)
                block_length = cds.end - cds.start
                current_frame = (current_frame + block_length) % 3
        elif bedline.strand == '-':
            for cds in reversed(bed6_blocks):
                frames.insert(0, current_frame)
                block_length = cds.end - cds.start
                current_frame = (current_frame + block_length) % 3
        # Adjusting the frames gets the right answer below. haven't figured out why. but i'm sure it makes the function work.
        frames = [2 if frame == 1 else 1 if frame == 2 else frame for frame in frames]
        return frames
    
def extract_codon(bedline, codon_type='start'):
    """
    Extract the start or stop codon from a BED12 object.
    
    Args:
        bedline (bedparse.bedline): A BED12 object.
        codon_type (str): Either 'start' or 'stop' to specify which codon to extract.
    
    Returns:
        bedparse.bedline: A BED12 object representing the codon.
    """
    if bedline.cds() is None:
        raise ValueError("The BED line does not contain any CDS region.")

    CDS_bed12 = bedline.cds()
    bed6_blocks = list(CDS_bed12.bed12tobed6())
    
    if codon_type == 'start':
        codon_blocks = []
        codon_length = 3
        if bedline.strand == '+':
            for block in bed6_blocks:
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.start, block.start + codon_length))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length
        elif bedline.strand == '-':
            for block in reversed(bed6_blocks):
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.end - codon_length, block.end))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length
            codon_blocks = list(reversed(codon_blocks))

    elif codon_type == 'stop':
        codon_blocks = []
        codon_length = 3
        if bedline.strand == '+':
            for block in reversed(bed6_blocks):
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.end - codon_length, block.end))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length
            codon_blocks = list(reversed(codon_blocks))
        elif bedline.strand == '-':
            for block in bed6_blocks:
                block_length = block.end - block.start
                if block_length >= codon_length:
                    codon_blocks.append((block.start, block.start + codon_length))
                    break
                else:
                    codon_blocks.append((block.start, block.end))
                    codon_length -= block_length

    else:
        raise ValueError("codon_type must be either 'start' or 'stop'.")

    block_sizes = ",".join([str(e - s) for s, e in codon_blocks])
    block_starts = ",".join([str(s - codon_blocks[0][0]) for s, _ in codon_blocks])

    return bedparse.bedline([
        bedline.chr,
        codon_blocks[0][0],
        codon_blocks[-1][1],
        bedline.name,
        bedline.score,
        bedline.strand,
        codon_blocks[0][0],
        codon_blocks[-1][1],
        bedline.color,
        len(codon_blocks),
        block_sizes,
        block_starts
    ])

def gtf_formatted_bedline_tx(bedline, source='.', attributes_str=''):
    """
    input bedparse.bedline transcript object, return 9 field gtf line string, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object 
    """
    gtf_fields = [bedline.chr, source, 'transcript', str(bedline.start + 1), str(bedline.end), str(bedline.score), bedline.strand, '.', attributes_str]
    return '\t'.join(gtf_fields) + '\n'

def gtf_formatted_bedline_exons(bedline, source='.', attributes_str=''):
    """
    input bedparse.bedline transcript object, return 9 field gtf line strings, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object 
    """
    string_to_return = ""
    for exon in bedline.bed12tobed6(appendExN=True):
        gtf_fields = [exon.chr, source, 'exon', str(exon.start + 1), str(exon.end), str(exon.score), exon.strand, '.', attributes_str + f' exon_id "{exon.name}";']
        string_to_return += '\t'.join(gtf_fields) + '\n'
    return string_to_return

def gtf_formatted_bedline_cds(bedline, source='.', attributes_str=''):
    """
    input bedparse.bedline transcript object, return 9 field gtf line strings, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object 
    """
    string_to_return = ""
    if bedline.cds():
        for cds, frame in zip(bedline.cds().bed12tobed6(), calculate_frames(bedline)):
            gtf_fields = [cds.chr, source, 'CDS', str(cds.start + 1), str(cds.end), str(cds.score), cds.strand, str(frame), attributes_str]
            string_to_return += '\t'.join(gtf_fields) + attributes_str + '\n'
    return string_to_return



def gtf_formatted_bedline_utr_start_stop(bedline, source='.', attributes_str='', filter_NF_in_attributes_str=True):
    """
    input bedparse.bedline transcript object, return 9 field gtf line strings, where the source and attributes fields are provided as optional argument, since those information cannot be inferred from the bedline object. Note that start and stop have frame attribute that may be non-zero if ATG cross exon boundary. In Gencode (and maybe ensembl... have to check), when a start codon traverses an exon junction, exists as two lines, one for each piece of the codon. If the attributes string contains 'NF' (as in 'mRNA_end_NF'), don't write out start or stop codon.
    """
    string_to_return = ""
    if filter_NF_in_attributes_str and '_NF' in attributes_str:
        return string_to_return
    bedline_cds = bedline.cds()
    if bedline_cds:
        start = extract_codon(bedline_cds, codon_type='start')
        for cds, frame in zip(start.bed12tobed6(), calculate_frames(start)):
            gtf_fields = [cds.chr, source, 'start_codon', str(cds.start + 1), str(cds.end), str(cds.score), cds.strand, str(frame), attributes_str]
            string_to_return += '\t'.join(gtf_fields) + '\n'
        stop = extract_codon(bedline_cds, codon_type='stop')
        for cds, frame in zip(stop.bed12tobed6(), calculate_frames(stop)):
            gtf_fields = [cds.chr, source, 'stop_codon', str(cds.start + 1), str(cds.end), str(cds.score), cds.strand, str(frame), attributes_str]
            string_to_return += '\t'.join(gtf_fields) + '\n'
        UTR_fiveprime = bedline.utr(which=5)
        if UTR_fiveprime:
            for exon in UTR_fiveprime.bed12tobed6():
                gtf_fields = [exon.chr, source, 'UTR', str(exon.start + 1), str(exon.end), str(exon.score), exon.strand, '.', attributes_str]
                string_to_return += '\t'.join(gtf_fields) + '\n'
        UTR_threeprime = bedline.utr(which=3)
        if UTR_threeprime:
            for exon in UTR_threeprime.bed12tobed6():
                gtf_fields = [exon.chr, source, 'UTR', str(exon.start + 1), str(exon.end), str(exon.score), exon.strand, '.', attributes_str]
                string_to_return += '\t'.join(gtf_fields) + '\n'
    return string_to_return

def bed12_formatted_bedline(bedline, attributes_str='', color=''):
    bed12_fields = [bedline.chr, bedline.start, bedline.end, bedline.name, bedline.score, bedline.strand, bedline.cdsStart, bedline.cdsEnd, color, bedline.nEx, bedline.exLengths, bedline.exStarts]
    return '\t'.join([str(i) for i in bed12_fields]) + attributes_str + '\n'
    
def get_transcript_length(bedline):
    return sum([int(i) for i in bedline.exLengths.split(',')])

def get_tx_stats(bedline, fasta_obj):
    """
    Return dict of potential useful stats/attributes about the transcript
    """
    dict_out = {"FiveUTR_nEx":0, "FiveUTR_len":0, "ThreeUTR_nEx":0, "ThreeUTR_len":0, "ThreeUTR_nEx_AfterFirst50":0, "CDSLen":0, "Introns":""}
    FiveUTR = bedline.utr(which=5)
    ThreeUTR = bedline.utr(which=3)
    CDS = bedline.cds()
    Introns = bedline.introns()
    if CDS:
        dict_out["CDSLen"] = get_transcript_length(CDS)
    if FiveUTR:
        dict_out["FiveUTR_nEx"] = FiveUTR.nEx
        dict_out["FiveUTR_len"] = get_transcript_length(FiveUTR)
    if ThreeUTR:
        dict_out["ThreeUTR_nEx"] = ThreeUTR.nEx
        dict_out["ThreeUTR_len"] = get_transcript_length(ThreeUTR)
        dict_out["ThreeUTR_nEx_AfterFirst50"] = ThreeUTR.nEx - ThreeUTR.extract_sequence(fasta_obj)[0:50].count("|")
    if Introns:
        dict_out["Introns"]=','.join([f'{i.chr}:{i.start}-{i.end}:{i.strand}' for i in Introns.bed12tobed6()])
    return dict_out

def get_thickStart_thickStop_from_marked_seq(bedline, ORF_marked_sequence):
    CDS = re.search(r"\^\w+\*", ORF_marked_sequence.replace("|", ""))
    if CDS:
        TranscriptLength = sum([int(i) for i in bedline.exLengths.split(',')])
        relativeStart = CDS.span()[0]
        relativeStop = CDS.span()[1] - 2
        FivePrimeEdge = get_absolute_pos(bedline, relativeStart)
        ThreePrimeEdge = get_absolute_pos(bedline, relativeStop)
        return tuple(sorted([FivePrimeEdge, ThreePrimeEdge]))
    return bedline.cdsStart, bedline.cdsEnd

def get_absolute_pos(bedline, coord):
    """
    Similar to bedline.tx2genome(coord, stranded=True), but I think that function is buggy and doesn't work for getting last base in tx
    """
    cumulative_block_size = 0
    exLengths = [int(i) for i in bedline.exLengths.split(',')]
    exStarts = [int(i) for i in bedline.exStarts.split(',')]
    if bedline.strand == '+':
        for i, block_size in enumerate(exLengths):
            remainder = coord - cumulative_block_size
            cumulative_block_size += block_size
            if cumulative_block_size >= coord:
                absolute_pos = bedline.start + remainder + exStarts[i]
                break
    elif bedline.strand == '-':
        for i, block_size in enumerate(list(reversed(exLengths))):
            remainder = coord - cumulative_block_size
            cumulative_block_size += block_size
            if cumulative_block_size >= coord:
                absolute_pos = bedline.start + list(reversed(exStarts))[i] + block_size - remainder
                break
    return absolute_pos

def add_gene_type_to_gtf(gtf_io, gene_type_dict):
    """
    Add gene_type attribute to a GTF file stored in a StringIO object.
    Parameters:
    gtf_io (StringIO): StringIO object containing the GTF file data.
    gene_type_dict (dict): Dictionary mapping gene_name to gene_type.
    Returns:
    StringIO: A new StringIO object with the updated GTF data.
    """
    # Create a new StringIO object to store the modified GTF content
    updated_gtf_io = StringIO()
    # Seek to the beginning of the input StringIO object
    gtf_io.seek(0)
    # Read and process each line
    for line in gtf_io:
        # Skip comment lines
        if line.startswith("#"):
            updated_gtf_io.write(line)
            continue
        # Parse the gene_id attribute using a regex
        match = re.search(r'gene_id\s+"([^"]+)"', line)
        if match:
            gene_id = match.group(1)
            gene_type = gene_type_dict.get(gene_id, "unknown_gene_type")  # Default to "unknown_gene_type" if not found
            # Add the gene_type attribute to the line
            if 'gene_type' not in line:
                line = line.rstrip() + f' gene_type "{gene_type}";\n'
        # Write the modified (or unmodified) line to the new StringIO object
        updated_gtf_io.write(line)
    # Reset the cursor of the new StringIO object to the beginning
    updated_gtf_io.seek(0)
    return updated_gtf_io



def parse_args(args=None):
    parser = argparse.ArgumentParser(description="Helper script to reformat GTF file for compatibility with SpliceJunctionClassifier.py. More specifically, for a GTF to be compatible with that script, each gene feature and transcript feature must have a 'gene_type' and 'gene_id' attributes where 'gene_type' attribute value is 'protein_coding' for protein coding genes, and 'gene_id' attribute value is unique for each gene. Each child transcript feature must also have 'transcript_type' and 'transcript_id' attributes, where 'transcript_type' is 'protein_coding' for protein_coding (productive) transcript isoforms. Also, each protein_coding isoform must have a child 'exon' 'CDS', 'start_codon' and 'stop_codon' features (which also have the gene_type, gene_id, transcript_type, and transcript_id attributes of their parent features). As in Gencode v43 GTFs, transcripts that are not 'protein_coding' may also have 'CDS', 'start_codon' and 'start_codon' child features but if the transcript_type attribute is not 'protein_coding', these will not determine inclusion into the set of productive start/stop codons by SpliceJunctionClassifier.py. Furthermore, this script adds NMDetectiveB attributes to the transcript (see options for this feature). This is useful in cases when a GTF may be missing 'transcript_type' attribute, and if it is even missing 'CDS' features, this script can add the 'CDS', 'start_codon', and 'stop_codon' features after attempting to translate the transcript sequence and assign values to the 'transcript_type' attribute (eg, 'protein_coding' if NMDetectiveB result is 'Last exon' and 'noncoding if NMDetectiveB result is 'Trigger NMD'). This script also optionally can output bed12 file for each transcript. I have tested this on some GTF files from Gencode, Ensembl, and UCSC. Depending on the exact nature (eg attribute names) in the input GTF, you may have to alter some options to make the output suitable for SpliceJunctionClassifier.py.")
    parser.add_argument('-i', dest='transcripts_in', required=True, help='Input file containing transcript structures. Can be GTF format or BED12 format (see -input_type). Gzipped files (.gz) are automatically detected and supported.')
    parser.add_argument('-o', dest='gtf_out', help='Output GTF file (optional if -bed12_out is specified)')
    parser.add_argument('-input_type', dest='input_type', choices=['gtf', 'bed12'], default='gtf', help='"-i" input file that contains transcript structures can be either gtf format, or bed12 format. If bed12 format, use --bed12_column_indexes to specify column mapping. default: "%(default)s"')
    parser.add_argument('-fa', dest='fasta_in', required=True, help='Input FASTA file')
    parser.add_argument('-bed12_out', dest='bed12_out', help='Optional bed12 out. One line per transcript. Attributes as extra columns')
    parser.add_argument('-bed12_column_indexes', dest='bed12_column_indexes', nargs=4, type=int, metavar=('TRANSCRIPT_ID', 'GENE_ID', 'TRANSCRIPT_TYPE', 'GENE_TYPE'), help='When input_type is bed12, specify 1-based column indexes for transcript_id, gene_id, transcript_type, and gene_type (in that order). Required when input_type=bed12. Example: --bed12_column_indexes 13 14 15 16 or --bed12_column_indexes 4 4 4 4 to use BED name column for all.')
    parser.add_argument('-n', dest='n_lines', help='Number of lines to read in gtf. Useful for quick debugging, but reading only n lines may cause buggy behavior if a transcript feature line is processed but not all of its child (eg exon, or CDS) features', type=int)
    parser.add_argument('-infer_transcript_type_approach', dest='infer_transcript_type_approach', choices=['A', 'B', 'C'], default='B', help='Approach to determine the transcript_type attribute value in output. (A) Use existing value. (B) Value is "protein_coding" if a the transcript contains child features of type "CDS". Value is "noncoding" otherwise. (C) Translate the transcript and use NMDetectiveB on the ORF. See -translation_approach and -NMDetectiveB_coding_threshold options to specify how this is done. default: "%(default)s"')
    parser.add_argument('-infer_gene_type_approach', dest='infer_gene_type_approach', choices=['A', 'B'], default='B', help='Approach to determine the gene_type attribute value in output. (A) Use existing value. (B) Value is "protein_coding" if gene has any child transcripts have attribute transcript_type "protein_coding". Value is "noncoding" otherwise. default: "%(default)s"')
    parser.add_argument('-transcript_id_attribute_name', dest='transcript_id_attribute_name', default='transcript_id', help='Name of the attribute in the input gtf file that defines the transcript_id attribute in the output file. Only relevant if -infer_transcript_type_approach==A. default: "%(default)s"')
    parser.add_argument('-gene_id_attribute_name', dest='gene_id_attribute_name', default='gene_id', help='Name of the attribute in the input gtf file that defines the gene_id attribute in the output file. Only relevant if -infer_gene_type_approach==A. default: "%(default)s"')
    parser.add_argument('-transcript_type_attribute_name', dest='transcript_type_attribute_name', default='transcript_type', help='Name of the attribute in the input gtf file that defines the transcript_type attribute in the output file. Only relevant if -infer_transcript_type_approach==A. default: "%(default)s"')
    parser.add_argument('-gene_type_attribute_name', dest='gene_type_attribute_name', default='gene_type', help='Name of the attribute in the input gtf file that defines the gene_type attribute in the output file. Only relevant if -infer_gene_type_approach==A. default: "%(default)s"')
    parser.add_argument('-extra_attributes', dest='extra_attributes', default='transcript_name,gene_name,transcript_support_level,tag,ccds_id',
                        help='Extra transcript attributes (comma delimited quoted string). default: "%(default)s"')
    parser.add_argument('-NMDetectiveB_coding_threshold', dest='NMDetectiveB_coding_threshold', type=int, choices=[1,2,3,4,5,6,7], default=5, help='NMDetectiveB classifies each transcript into 7 ordinal categories, from the most coding potential to the least coding potential: (1) Last exon (2) Start proximal (3) 50nt rule (4) Long exon (5) Trigger NMD (6) No stop (7) No CDS. Transcripts with classified as this NMDetective value or greater will be assigned transcript_type attribute value of "noncoding", while others will be value of "protein_coding". default: "%(default)s"')
    parser.add_argument('-translation_approach', dest='translation_approach', choices=['A', 'B', 'C', 'D', 'E', 'F'], default='B', help='Approach to use NMDetective to annotate CDS in output for genes where gene_biotype/gene_type == "protein_coding", some of which may not have annotated CDS in input (eg transcript_type=="processed_transcript"). Possible approaches: (A) using annotated ORF if ORF is present in input with 5UTR and 3UTR. If 3UTR is present and 5UTR is absent (suggesting stop codon is annotated but start codon may be outside of the transcript bounds), search for longest ORF within transcript bounds where a start codon is not required at the beginning of ORF. Similarly, if 5UTR is present but 3UTR is absent, search for longest ORF without requiring stop codon. If neither UTR is present, or if no CDS is annotated in input, search for longest ORF, not requiring start or stop within transcript bounds. I think this approach might be useful to correctly identify the ORF, even if transcript bounds are not accurate, but it has the downside that true "processed_transcripts" with a internal TSS that eliminates the correct start codon, may be erronesously be classified as "Last exon" (ie productive) transcripts by NMDFinderB. (B) Use only annotated CDS. In effect, output gtf is the same except start_codon and stop_codon features are added even if not present in input. This would be useful for dealing with "processed_transcripts" properly by NMDFinder, but I havent checked whether CDS annotations in poorly annotated species (eg lamprey, chicken, etc) are reasonable, which could be a problem for Yangs script. (C) use annotated CDS if present, and use first ATG if no CDS present (minimum ORF length controlled by --min_new_ORF_length). (D) Use first ORF regardless of annotation (minimum ORF length controlled by --min_new_ORF_length). (E) Use longest ORF with both start and stop codons required (minimum ORF length controlled by --min_new_ORF_length). (F) Use longest ORF with start codon required but no stop codon required (minimum ORF length controlled by --min_new_ORF_length) - useful for detecting no-stop decay in full-length reads.')
    parser.add_argument('--min_new_ORF_length', dest='min_new_ORF_length', type=int, default=50, help='Minimum ORF length (in codons) for newly predicted ORFs in translation approaches C, D, E, and F. Only relevant when using translation approach C, D, E, or F. default: %(default)s')
    parser.add_argument('--include_uorf_analysis', action='store_true', help='Include upstream ORF (uORF) analysis in output. Adds uORF_start_genomic_positions, uORF_stop_genomic_positions, uORF_classifications, uORF_lengths_aa, and uORF_relative_positions attributes.', default=False)
    parser.add_argument('-v', '--verbose', action='store_true', help='Increase verbosity')
    parser.add_argument('--sort_output', action='store_true', default=False,
                        help='Sort outputs before writing: BED12 will be coordinate-sorted; GTF will be sorted by gene, then transcript, then features within transcript. If not set, BED12 preserves read order and GTF preserves construction order.')
    parser.add_argument('--kozak_jaspar_pfm', dest='kozak_jaspar_pfm', default=None,
                        help='Optional path to a JASPAR-like PFM file to override the default Kozak motif. If provided, its PWM/PSSM will be used for scoring.')
    return parser.parse_args(args)


def setup_logging(verbose):
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format='%(asctime)s - %(levelname)s - %(message)s')

def _open_bgzip_aware(path, text_mode=True):
    """Open path for writing, using BGZF for .gz/.bgz outputs.
    text_mode=True uses text mode; otherwise binary."""
    mode = 'wt' if text_mode else 'wb'
    if path.endswith('.gz') or path.endswith('.bgz'):
        return bgzf.open(path, mode)
    # fallback to regular open
    return open(path, 'w' if text_mode else 'wb')

def main(args=None):
    args = parse_args(args)
    setup_logging(args.verbose)

    # Check that at least one output is specified
    if not args.gtf_out and not args.bed12_out:
        raise ValueError("At least one output must be specified: -o (GTF output) or -bed12_out (BED output)")

    # Validation for BED12 input
    if args.input_type == "bed12" and not args.bed12_column_indexes:
        raise ValueError("--bed12_column_indexes is required when input_type is bed12")

    # Example usage of the parsed arguments
    logging.info("Input transcript file: %s", args.transcripts_in)
    if args.gtf_out:
        logging.info("Output GTF file: %s", args.gtf_out)
    else:
        logging.info("GTF output: disabled")
    if args.bed12_out:
        logging.info("Output BED file: %s", args.bed12_out)
    else:
        logging.info("BED output: disabled")
    logging.info("Input FASTA file: %s", args.fasta_in)
    # Handle extra attributes based on input type
    if args.input_type == "gtf" and args.extra_attributes:
        extra_attributes = args.extra_attributes.split(',')
        logging.info("Extra attributes: %s", extra_attributes)
    elif args.input_type == "bed12":
        extra_attributes = []  # BED12 input doesn't have GTF extra attributes
        logging.info("BED12 input: no extra attributes")
    else:
        extra_attributes = []
        logging.info("No extra attributes provided.")

    # Open the FASTA file using pyfastx
    fasta_obj = pyfastx.Fasta(args.fasta_in)

    # Handle attribute extraction based on input type and inference approaches
    if args.input_type == "gtf":
        # Build attributes to extract - handle missing transcript_type/gene_type gracefully
        required_attributes = [args.transcript_id_attribute_name, args.gene_id_attribute_name]
        
        # Always try to extract type attributes (user's responsibility to use inference approaches)
        required_attributes.extend([args.transcript_type_attribute_name, args.gene_type_attribute_name])
        
        attributes_to_extract = ','.join(required_attributes + extra_attributes)
        
        # gtf2bed
        logging.info(f'Reading in transcripts and converting transcripts to bedparse.bedline objects...')
        bed12 = run_bedparse_gtf2bed(args.transcripts_in, '--extraFields', attributes_to_extract, n=args.n_lines)
        
    elif args.input_type == "bed12":
        logging.info(f'Reading BED12 input with column mapping: {args.bed12_column_indexes}')
        
        # Handle -n argument for BED12 input with gzip auto-detection
        if args.n_lines is not None:
            bed12_lines = []
            open_func = gzip.open if args.transcripts_in.endswith('.gz') else open
            mode = 'rt' if args.transcripts_in.endswith('.gz') else 'r'
            with open_func(args.transcripts_in, mode) as f:
                for i, line in enumerate(f):
                    if i >= args.n_lines:
                        break
                    bed12_lines.append(line.rstrip('\n'))
            bed12 = '\n'.join(bed12_lines)
            logging.info(f'Reading first {args.n_lines} lines from BED12 file')
        else:
            open_func = gzip.open if args.transcripts_in.endswith('.gz') else open
            mode = 'rt' if args.transcripts_in.endswith('.gz') else 'r'
            with open_func(args.transcripts_in, mode) as f:
                bed12 = f.read()
        
        # Validate column indexes against actual BED12 data
        first_line = bed12.split('\n')[0].split('\t')
        max_col_needed = max(args.bed12_column_indexes)
        if len(first_line) < max_col_needed:
            raise ValueError(f"BED12 file has {len(first_line)} columns but column index {max_col_needed} was specified")
        
        logging.info(f"BED12 column mapping: transcript_id=col{args.bed12_column_indexes[0]}, gene_id=col{args.bed12_column_indexes[1]}, transcript_type=col{args.bed12_column_indexes[2]}, gene_type=col{args.bed12_column_indexes[3]}")
        
    NumTrancsripts = len(bed12.splitlines())
    logging.info(f'Read in {NumTrancsripts} from input.')

    logging.info(f'Processing transcripts.')
    
    # Initialize GTF-specific variables only if GTF output is requested
    if args.gtf_out:
        gtf_stringio = StringIO()
        # Initialize an empty dictionary to keep coordinates of each transcript for each gene. Will need to write out gene coordinates as most extensive span of child transcripts on each chromosome:strand pair, because some genes appear on more than one chromosome, (eg X and Y in Gencode gtf), or even on the same chromosome but different strands (eg some TRNAA gene appears on chr6 + strand and chr6 - strand in human RefSeq annotations from UCSC)
        gene_coords_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda : defaultdict(set))))
        # gene_coords_dict[gene_id][transcript.chr][transcript.strand]['end'] = SetOfTranscriptEnds
    else:
        gtf_stringio = None
        gene_coords_dict = None
    
    # Always initialize gene_types_dict - needed for both GTF and BED output
    gene_types_dict = defaultdict(lambda: defaultdict((set)))

    # BED output setup (sorting controlled by --sort_output)
    bed_out_fh = None
    bed_header_line = None
    sort_bed = bool(args.bed12_out) and args.sort_output
    if args.bed12_out:
        # Build BED header
        bed_header_fields = [
            "chr", "start", "end", "name", "score", "strand", 
            "thickStart", "thickEnd", "color", "blockCount", "blockSizes", "blockStarts",
            "gene_id", "transcript_id", "gene_type", "transcript_type", "NMDFinderB"
        ]
        extra_calc_attrs = ["FiveUTR_nEx", "FiveUTR_len", "ThreeUTR_nEx", "ThreeUTR_len", 
                           "ThreeUTR_nEx_AfterFirst50", "CDSLen", "Introns"]
        bed_header_fields.extend(extra_calc_attrs)
        if args.include_uorf_analysis:
            # Include separate start and stop genomic position columns for uORFs
            bed_header_fields.extend(["uORF_start_genomic_positions", "uORF_stop_genomic_positions", "uORF_classifications", "uORF_lengths_aa", "uORF_relative_positions", "uORF_kozak_log_odds", "start_codon_kozak_log_odds"])
        else:
            bed_header_fields.extend(["start_codon_kozak_log_odds"])
        bed_header_line = "#" + "\t".join(bed_header_fields) + "\n"
        if sort_bed:
            bed_lines = []
        else:
            bed_out_fh = _open_bgzip_aware(args.bed12_out, text_mode=True)
            bed_out_fh.write(bed_header_line)

    # Initialize Kozak PSSM once (override with JASPAR if provided)
    pssm = None
    try:
        if args.kozak_jaspar_pfm:
            _, pssm = load_jaspar_pfm_pssm(args.kozak_jaspar_pfm)
            logging.info(f"Loaded Kozak PFM from {args.kozak_jaspar_pfm}")
        else:
            # Use your empirical counts to synthesize instances and build PSSM
            counts_dict = {
                'A': [3620.0, 3166.0, 4450.0, 8952.0, 5277.0, 3288.0, 18870.0, 0.0, 2.0, 3973.0],
                'C': [4423.0, 6320.0, 7562.0, 1691.0, 7700.0, 8999.0, 25.0, 3.0, 0.0, 2652.0],
                'G': [7701.0, 5952.0, 4862.0, 7331.0, 3652.0, 5302.0, 2.0, 0.0, 18895.0, 9712.0],
                'T': [3155.0, 3461.0, 2025.0, 925.0, 2270.0, 1310.0, 2.0, 18896.0, 2.0, 2562.0],
            }
            _, _, pssm = build_pssm_from_counts_by_sampling(counts_dict, n_instances=800, seed=1)
    except Exception as e:
        logging.error(f"Failed to initialize Kozak PSSM: {e}")
        pssm = None

    for i,l in enumerate(bed12.splitlines()):
        if i >= 0:
            if i % 1000 == 0: logging.debug(f'processed {i} trancsripts for output...')
            lsplit = l.split('\t')
            
            # Handle attribute extraction based on input type
            if args.input_type == "gtf":
                # GTF input: attributes come from bedparse gtf2bed extraction (columns 12-15 + extras)
                transcript_id, gene_id, transcript_type, gene_type = lsplit[12:16]
                extra_attribute_values = lsplit[16:]
            elif args.input_type == "bed12":
                # BED12 input: use column mapping (convert from 1-based to 0-based indexing)
                transcript_id = lsplit[args.bed12_column_indexes[0] - 1]  # transcript_id
                gene_id = lsplit[args.bed12_column_indexes[1] - 1]        # gene_id 
                transcript_type = lsplit[args.bed12_column_indexes[2] - 1]  # transcript_type
                gene_type = lsplit[args.bed12_column_indexes[3] - 1]        # gene_type
                extra_attribute_values = []  # BED12 doesn't have extra attributes by default
            transcript = bedparse.bedline(lsplit[0:12])
            if not Is_bedline_complete(transcript):
                continue
            source = "input_gtf"
            transcript_attributes = f'transcript_id "{transcript_id}"; gene_id "{gene_id}";'
            transcript_out = transcript
            NMDFinderB = "NA"
            #Determine NMDetectiveB classification
            try:
                if args.translation_approach == 'A':
                    if transcript.cds() and transcript.utr(which=5) and transcript.utr(which=3):
                        source = "input_gtf"
                        sequence = transcript.extract_sequence(fasta_obj, AddMarksForORF=True)
                    elif transcript.cds() and transcript.utr(which=5) and not transcript.utr(which=3):
                        source = "LongestORF_NoStopRequired"
                        sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    elif transcript.cds() and not transcript.utr(which=5) and transcript.utr(which=3):
                        source = "LongestORF_NoStartRequired"
                        sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    else:
                        source = "LongestORF_NeitherRequired"
                        sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_ATG=False, require_STOP = False)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence)
                elif args.translation_approach == 'B':
                    sequence = transcript.extract_sequence(fasta_obj, AddMarksForORF=True)
                    NMDFinderB = get_NMD_detective_B_classification(sequence)
                elif args.translation_approach == 'C':
                    if transcript.cds():
                        sequence = transcript.extract_sequence(fasta_obj, AddMarksForORF=True)
                        source = "input_gtf"
                    else:
                        source = "FirstORF_NoStopRequired"
                        sequence = insert_marks_for_first_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False, min_ORF_len = args.min_new_ORF_length)
                        thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                        transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence)
                elif args.translation_approach == 'D':
                    source = "FirstORF_NoStopRequired"
                    sequence = insert_marks_for_first_ORF(transcript.extract_sequence(fasta_obj), require_STOP = False, min_ORF_len = args.min_new_ORF_length)
                    thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                    transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence)
                elif args.translation_approach == 'E':
                    source = "LongestORF_WithMinLength"
                    sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_ATG=True, require_STOP=True, min_ORF_len=args.min_new_ORF_length)
                    thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                    transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence)
                elif args.translation_approach == 'F':
                    source = "LongestORF_NoStopRequired"
                    sequence = insert_marks_for_longset_ORF(transcript.extract_sequence(fasta_obj), require_ATG=True, require_STOP=False, min_ORF_len=args.min_new_ORF_length)
                    thickStart, thickStop = get_thickStart_thickStop_from_marked_seq(transcript, sequence)
                    transcript_out = bedparse.bedline([transcript.chr, transcript.start, transcript.end, transcript.name, transcript.score, transcript.strand, thickStart, thickStop, transcript.color, transcript.nEx, transcript.exLengths, transcript.exStarts])
                    NMDFinderB = get_NMD_detective_B_classification(sequence)
            except NameError as e:
                logging.error(f"NameError encountered while processing transcript {transcript.name}: {e}")
                continue
            NMDFinderB_NoWhitespace = NMDFinderB.replace(' ', '_')
            NMDFinderB_number = get_NMD_detective_B_classification_number(NMDFinderB)
            # Determine transcript_type
            if args.infer_transcript_type_approach == 'A':
                transcript_type_out = transcript_type
            elif args.infer_transcript_type_approach == 'B':
                transcript_type_out = "protein_coding" if transcript.cds() else "noncoding"
            elif  args.infer_transcript_type_approach == 'C':
                transcript_type_out = "protein_coding" if NMDFinderB_number < args.NMDetectiveB_coding_threshold else "noncoding"
            extra_calculated_transcript_attributes = get_tx_stats(transcript_out, fasta_obj)
            
            # Build all transcript attributes in one place
            transcript_attributes = f'transcript_id "{transcript_id}"; gene_id "{gene_id}"; transcript_type "{transcript_type_out}"; tag "NMDFinderB:{NMDFinderB_NoWhitespace}";'

            # Add extra calculated attributes
            extra_calculated_transcript_attributes = get_tx_stats(transcript_out, fasta_obj)
            for k, v in extra_calculated_transcript_attributes.items():
                transcript_attributes += f' tag "{k}":"{v}";'

            # Conditionally add upstream ORF analysis
            if args.include_uorf_analysis:
                uorf_start_genomic_pos, uorf_stop_genomic_pos, uorf_classes, uorf_lengths, uorf_relative_pos, uorf_kozak_scores = Analyze_uORFs(sequence, transcript_out, pssm=pssm)
                transcript_attributes += f' tag "uORF_start_genomic_positions":"{uorf_start_genomic_pos}";'
                transcript_attributes += f' tag "uORF_stop_genomic_positions":"{uorf_stop_genomic_pos}";'
                transcript_attributes += f' tag "uORF_classifications":"{uorf_classes}";'
                transcript_attributes += f' tag "uORF_lengths_aa":"{uorf_lengths}";'
                transcript_attributes += f' tag "uORF_relative_positions":"{uorf_relative_pos}";'
                transcript_attributes += f' tag "uORF_kozak_log_odds":"{uorf_kozak_scores}";'
            # Always compute main start_codon Kozak score
            start_codon_kozak = score_kozak_window(sequence, pssm=pssm)
            transcript_attributes += f' tag "start_codon_kozak_log_odds":"{start_codon_kozak}";'

            # After transcript_out is finalized, track per-gene coordinate span for GTF gene features
            if args.gtf_out and gene_coords_dict is not None:
                chrom = transcript_out.chr
                strand = transcript_out.strand
                # initialize nested sets
                _ = gene_coords_dict[gene_id][chrom][strand]['start']
                _ = gene_coords_dict[gene_id][chrom][strand]['end']
                gene_coords_dict[gene_id][chrom][strand]['start'].add(transcript_out.start)
                gene_coords_dict[gene_id][chrom][strand]['end'].add(transcript_out.end)

            # Track gene types and transcript types seen (for later gene_type inference)
            gtinfo = gene_types_dict[gene_id]
            if 'gene_types_in_input' not in gtinfo:
                gtinfo['gene_types_in_input'] = set()
            if 'trancscript_types' not in gtinfo:
                gtinfo['trancscript_types'] = set()
            if gene_type:
                gtinfo['gene_types_in_input'].add(gene_type)
            if transcript_type:
                gtinfo['trancscript_types'].add(transcript_type)

            # Write to GTF only if GTF output is enabled
            if args.gtf_out:
                _ = gtf_stringio.write(gtf_formatted_bedline_tx(transcript_out, source=source, attributes_str=transcript_attributes))
                _ = gtf_stringio.write(gtf_formatted_bedline_exons(transcript_out, source=source, attributes_str=transcript_attributes))
                _ = gtf_stringio.write(gtf_formatted_bedline_cds(transcript_out, source=source, attributes_str=transcript_attributes))
                _ = gtf_stringio.write(gtf_formatted_bedline_utr_start_stop(transcript_out, source=source, attributes_str=transcript_attributes))

            # For BED output, create a clean list of values
            if args.bed12_out:
                bed_extra_fields = [gene_id, transcript_id, gene_type, transcript_type_out, NMDFinderB_NoWhitespace]
                bed_extra_fields.extend([str(v) for v in extra_calculated_transcript_attributes.values()])
                if args.include_uorf_analysis:
                    bed_extra_fields.extend([uorf_start_genomic_pos, uorf_stop_genomic_pos, uorf_classes, uorf_lengths, uorf_relative_pos, uorf_kozak_scores, start_codon_kozak])
                else:
                    bed_extra_fields.extend([start_codon_kozak])
                bed_line = bed12_formatted_bedline(transcript_out, 
                                                   color=get_NMD_detective_B_classification_color(NMDFinderB), 
                                                   attributes_str='\t' + '\t'.join(bed_extra_fields))
                if sort_bed:
                    bed_lines.append(bed_line)
                else:
                    _ = bed_out_fh.write(bed_line)
    if args.bed12_out:
        if sort_bed:
            with _open_bgzip_aware(args.bed12_out, text_mode=True) as fh:
                fh.write(bed_header_line)
                for line in sorted(bed_lines, key=lambda s: (s.split('\t')[0], int(s.split('\t')[1]), int(s.split('\t')[2]))):
                    fh.write(line)
        else:
            if bed_out_fh is not None:
                bed_out_fh.close()

    # Only process GTF output if it was requested
    if args.gtf_out:
        logging.info('Writing gene (parent) feature coordiantes and based on child (transcript) feature coordinates')
        for gene, chrom_dict in gene_coords_dict.items():
            for chrom, strand_dict in chrom_dict.items():
                for StrandIteration, (strand, info_dict) in enumerate(strand_dict.items()):
                    min_start = min(info_dict['start'])
                    max_stop = max(info_dict['end'])
                    _ = gtf_stringio.write(f'{chrom}\tinput_gtf\tgene\t{min_start+1}\t{max_stop}\t.\t{strand}\t.\tgene_id "{gene}";\n')
                if StrandIteration > 0:
                    logging.warning(f"Transcripts for gene {gene} on {chrom} are on different strands. Writing {gene} gene feature on {chrom} for more than one strand")

        logging.info('Writing gene_type attributes')
        gene_types_dict_final = dict()
        for gene, info_dict in gene_types_dict.items():
            if args.infer_gene_type_approach == "A":
                if len(info_dict['gene_types_in_input']) != 1: 
                    raise ValueError(f"{gene} does not have a single gene_type in input")
                gene_types_dict_final[gene] = list(info_dict['gene_types_in_input'])[0]
            elif args.infer_gene_type_approach == "B":
                gene_types_dict_final[gene] = "protein_coding" if "protein_coding" in info_dict['trancscript_types'] else "noncoding"
        # Write out gene_type attributes
        gtf_stringio_updated = add_gene_type_to_gtf(gtf_stringio, gene_types_dict_final)

        logging.info('writing out gtf')
        with _open_bgzip_aware(args.gtf_out, text_mode=True) as output_fh:
            if args.input_type == "gtf":
                input_open_func = gzip.open if args.transcripts_in.endswith('.gz') else open
                input_mode = 'rt' if args.transcripts_in.endswith('.gz') else 'r'
                with input_open_func(args.transcripts_in, input_mode) as input_fh:
                    for l in input_fh:
                        if l.startswith('#'):
                            _ = output_fh.write(l)
                        else:
                            break
            _ = output_fh.write(f"#! args: {args}\n")
            # If not sorting by coordinates, reorder by gene→transcript→features
            if not args.sort_output:
                reorder_gtf(gtf_stringio_updated, output_fh, mode='a')
            else:
                # Coordinate sort by seqname, start, end
                columns = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
                dtype_map = {
                    'seqname': 'string', 'source': 'string', 'feature': 'string',
                    'start': 'int64', 'end': 'int64', 'score': 'string',
                    'strand': 'string', 'frame': 'string', 'attribute': 'string'
                }
                df = pd.read_csv(gtf_stringio_updated, sep='\t', names=columns, comment='#', dtype=dtype_map)
                df_sorted = df.sort_values(by=['seqname', 'start', 'end']).reset_index(drop=True)
                df_sorted.to_csv(output_fh, sep='\t', header=False, index=False, mode='a', quoting=csv.QUOTE_NONE)
    else:
        logging.info('GTF output disabled - skipping GTF processing and sorting')

if __name__ == "__main__":
    main()