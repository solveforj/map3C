import re
import pysam
import pandas as pd
import numpy as np
import subprocess
import shlex
import gzip
import bisect
import os
from collections import OrderedDict

rng = np.random.default_rng(1)

def get_mate_from_tag(read):
    if read.is_read1:
        return "1"
    elif read.is_read2:
        return "2"
    else:
        return "1"
        #raise Exception(f"Mate not defined for read {read.query_name}")


def process_chrom_sizes(chrom_sizes_file):
    chrom_sizes = {}

    with open(chrom_sizes_file) as csf:
        for line in csf:
            line = line.strip().split()
            chrom = line[0]
            size = int(line[1])
            chrom_sizes[chrom] = size
    return chrom_sizes

def process_chrom_orders(chrom_sizes_file):
    chrom_orders = {}
    order = 0
        
    with open(chrom_sizes_file) as csf:
        for line in csf:
            line = line.strip().split()
            chrom = line[0]
            order += 1
            chrom_orders[chrom] = order
    return chrom_orders

def process_bed(bed):
    bed_dict = {}

    #if not os.path.isfile(bed):
    #    raise Exception("Provided BED file path does not exist")
            
    if bed.endswith(".gz"):
        f = gzip.open(bed, 'rt')
    else:
        f = open(bed)
        
    for line in f:
        line = line.strip().split()[:3]
        chrom = line[0]
        start = int(line[1])
        end = int(line[2])
        if chrom not in bed_dict:
            bed_dict[chrom] = []
        bed_dict[chrom].append({"start":start, "end":end, "idx":None})
    for chrom in bed_dict:
        bed_dict[chrom] = sorted(bed_dict[chrom], key=lambda x:x["start"])
        for i in range(len(bed_dict[chrom])):
            bed_dict[chrom][i]["idx"] = i
    return bed_dict

def process_variants(snps):
    
    snp_dict = {}

    #if not os.path.isfile(snps):
    #    raise Exception("Provided variants file path does not exist")

    #if not os.path.isfile(snps):
    #    raise Exception("Provided variants file path does not exist")

    if snps.endswith(".gz"):
        f = gzip.open(snps, 'rt')
    else:
        f = open(snps)
        
    for line in f:
        line = line.strip().split()
        chrom = line[0]
        pos = int(line[1]) - 1
        a1 = line[2]
        a2 = line[3]
        if chrom not in snp_dict:
            snp_dict[chrom] = []
        snp_dict[chrom].append({"pos":pos, "a1":a1, "a2":a2})
    for chrom in snp_dict:
        snp_dict[chrom] = sorted(snp_dict[chrom], key=lambda x : x["pos"])
    
    return snp_dict
            
    
def process_restriction_sites(restriction_sites, restriction_enzymes, chrom_sizes):

    restriction_sites_dict = {}
    restriction_sites = sum(restriction_sites, [])
    restriction_enzymes = sum(restriction_enzymes, [])

    if len(restriction_sites) != len(restriction_enzymes):
        raise Exception("You didn't specify the same number of restriction enzyme names and site position files!")
        
    if len(restriction_sites) == 0:
        return restriction_sites_dict
            
    for i in range(len(restriction_sites)):
        file = restriction_sites[i]
        enzyme = restriction_enzymes[i]
        #enzyme = os.path.basename(file).split("_")[-1].split(".txt")[0]
        restriction_sites_dict[enzyme] = {}
        with open(file) as f:
            for line in f:
                line = line.strip().split()
                chrom = line[0]
                sites = [int(i) for i in line[1:]]
                if chrom not in restriction_sites_dict[enzyme]:
                    restriction_sites_dict[enzyme][chrom] = sites
                else:
                    restriction_sites_dict[enzyme][chrom] += sites

    return restriction_sites_dict

class PairsStats:

    def parse_base(self, line, contact_type):
        chrom1 = line[1]
        chrom2 = line[3]
        pos1 = line[2]
        pos2 = line[4]

        if chrom1 == chrom2:
            dist = np.abs(int(pos2) - int(pos1))
            if dist >= 1000:
                self.pairs_stats[f"{contact_type}_intra1kb"] += 1
            if dist >= 10000:
                self.pairs_stats[f"{contact_type}_intra10kb"] += 1
            if dist >= 20000:
                self.pairs_stats[f"{contact_type}_intra20kb"] += 1
        else:
            self.pairs_stats[f"{contact_type}_inter"] += 1

        self.pairs_stats[f"{contact_type}_total"] += 1


    def parse_phase(self, line, contact_type):
        phase1 = line[self.phase1_idx]
        phase2 = line[self.phase2_idx]

        phase_id = "".join(sorted([phase1, phase2]))

        self.pairs_stats[f"{contact_type}_phased_{phase_id}"] += 1

    def parse_metadata(self, line, contact_type):
        pair_type = line[self.pair_type_idx]
        if pair_type == "RU":
            pair_type = "UR"
        rule = line[self.rule_idx]
        reads = line[self.reads_idx]

        pair_type_rule = f"{pair_type}_{rule}"

        self.pairs_stats[f"{contact_type}_{reads}"] += 1
        self.pairs_stats[f"{contact_type}_{pair_type_rule}"] += 1
        
    def parse_file(self):
        self.parsers = [self.parse_base]
        
        with gzip.open(self.pair_file, "rt") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#columns"):
                    columns = line[10:].split()

                    self.pair_type_idx = columns.index("pair_type")

                    if "rule" in columns:
                        self.rule_idx = columns.index("rule")
                        self.reads_idx = columns.index("reads")
                        self.contact_class_idx = columns.index("contact_class")
                        #multimap_overlap_idx = columns.index("multimap_overlap")
                        #cs_locs_idx = columns.index("cut_site_locs")

                        self.types = ["enzyme", "enzymeless"]
                        
                        for i in self.types:
                            self.pairs_stats.update({
                                f"pairs_{i}_R1" : 0,
                                f"pairs_{i}_R2" : 0,
                                f"pairs_{i}_R1_rescue" : 0,
                                f"pairs_{i}_R2_rescue" : 0,
                                f"pairs_{i}_R1-2" : 0,
                                f"pairs_{i}_R1&2" : 0,
                                f"pairs_{i}_comb" : 0,
                                f"pairs_{i}_UU_all" : 0,
                                f"pairs_{i}_UU_mask" : 0,
                                f"pairs_{i}_UR_mask" : 0
                            })
                        
                        self.parsers.append(self.parse_metadata)
                    else:
                        self.types = ["pairs"]
                        
                    if "phase0" in columns and "phase1" in columns:
                        self.phase1_idx = columns.index("phase0")
                        self.phase2_idx = columns.index("phase1")
                        for i in self.types:
                            self.pairs_stats.update({
                                f"pairs_{i}_phased_.0" : 0,
                                f"pairs_{i}_phased_.1" : 0,
                                f"pairs_{i}_phased_.." : 0,
                                f"pairs_{i}_phased_11" : 0,
                                f"pairs_{i}_phased_00" : 0,
                                f"pairs_{i}_phased_01" : 0,
                            })
                        
                        self.parsers.append(self.parse_phase)

                    for i in self.types:
                        self.pairs_stats.update({
                            f"pairs_{i}_total" : 0,
                            f"pairs_{i}_intra1kb" : 0,
                            f"pairs_{i}_intra10kb" : 0,
                            f"pairs_{i}_intra20kb" : 0,
                            f"pairs_{i}_inter" : 0})
                                    
                        
                elif line.startswith("#"):
                    continue

                else:
                    line = line.split()
                    if len(self.types) == 1:
                        contact_type = self.types[0]
                    else:
                        contact_type = "pairs_enzymeless" if "enzymeless" in line[self.contact_class_idx] else "pairs_enzyme"

                    for parser in self.parsers:
                        parser(line, contact_type)
                    
    def parse_dedup_file(self):
        pairs_dup_rate = np.nan
        with open(self.pairtools_dedup_file) as f:
            for line in f:
                if "summary/frac_dups" in line:
                    try:
                        pairs_dup_rate = float(line.strip().split()[1])
                    except ValueError:
                        pairs_dup_rate = np.nan
                    break

        self.pairs_stats.update({f"pairs_dup_rate":pairs_dup_rate})
                        
    def __init__(self, pair_file, pairtools_dedup_file=None):
        
        self.pair_file = pair_file
        
        self.pairs_stats = {}
        
        self.parse_file()

        if pairtools_dedup_file:
            self.pairtools_dedup_file = pairtools_dedup_file
            self.parse_dedup_file()


def pair_stats(out_prefix,
                    input_pairs,
                    pairs_dedup_stats=None,
                    pairs_filterbycov_stats=None
                   ):

    stats_path = f"{out_prefix}_pairs_stats.txt"
    
    contacts_stats = PairsStats(input_pairs, pairs_dedup_stats).pairs_stats

    if pairs_filterbycov_stats:
        if os.path.isfile(pairs_filterbycov_stats):
            with open(pairs_filterbycov_stats) as f:
                highcov = 0.0
                for line in f:
                    if "pair_types/FF" in line:
                        try:
                            highcov = float(line.strip().split()[1])
                        except ValueError:
                            highcov = 0.0
                        break
            contacts_stats["high_coverage_pairs"] = highcov
                
    full_stats = {}

    full_stats.update(contacts_stats)
    
    stats_df = pd.DataFrame.from_dict(full_stats, orient="index").T
    stats_df.to_csv(stats_path, index=False, sep="\t")

def compute_genome_coverage(bam, min_mapq, min_base_quality, keep_dup=False):
    # Compute genome coverage

    ff_val = "QCFAIL,UNMAP"
    if not keep_dup:
        ff_val = "DUP," + ff_val

    command1 = f"samtools coverage -q {min_mapq} -Q {min_base_quality} -H --ff {ff_val} {bam}"
    command2 = "awk -v OFS='\t' -v nn=0 '{{nn += $5}} END {print nn } '"
    p1 = subprocess.Popen(shlex.split(command1), stdout=subprocess.PIPE)
    p2 = subprocess.run(shlex.split(command2), stdin=p1.stdout, capture_output=True, text=True)
    out = p2.stdout.strip()
    reference_coverage = int(out)
    return reference_coverage

def compute_mapped_nucleotides(bam, min_mapq, min_base_quality, keep_dup=False):
    # Compute mapped nucleotides

    add_flags = "SECONDARY"
    if keep_dup:
        add_flags += ",DUP"
    
    command1 = f"samtools depth -Q {min_mapq} -q {min_base_quality} -g {add_flags} {bam}"
    command2 = "awk -v nn=0 '{nn+=$3} END { print nn}'"
    p1 = subprocess.Popen(shlex.split(command1), stdout=subprocess.PIPE)
    p2 = subprocess.run(shlex.split(command2), stdin=p1.stdout, capture_output=True, text=True)
    out = p2.stdout.strip()
    mapped_read_bases = int(out)
    return mapped_read_bases

def aggregate_qc_stats(job,
                       out_prefix, 
                       mode):

    if mode == "bsdna":
        
        txt_paths = [f"{out_prefix}_trim_stats.txt",
                     f"{out_prefix}_contam_stats.txt",
                     f"{out_prefix}_dupsifter_stats.txt", 
                     f"{out_prefix}_alignment_stats.txt",
                     f"{out_prefix}_pairs_stats.txt",
                     f"{out_prefix}_methylation_stats.txt"
                    ]

    elif mode == "dna":
        
        txt_paths = [f"{out_prefix}_trim_stats.txt",
                     f"{out_prefix}_alignment_stats.txt",
                     f"{out_prefix}_pairs_stats.txt",
                    ]

    elif mode == "snmCTseq":

        txt_paths = [f"{out_prefix}_trim_stats.txt",

                     # DNA
                     f"{out_prefix}_biscuit_mod_contam_stats.txt",
                     f"{out_prefix}_biscuit_dupsifter_stats.txt",
                     f"{out_prefix}_biscuit_alignment_stats.txt",
                     f"{out_prefix}_methylation_stats.txt",

                     # RNA
                     f"{out_prefix}_STAR_contam_stats.txt",
                     f"{out_prefix}_picard_rna_stats.txt",
                     f"{out_prefix}_STAR_stats.txt",
                     f"{out_prefix}_featureCount_stats.txt",
                    ]
    
    stat_dfs = [pd.DataFrame([job], columns=["job"])]
    for path in txt_paths:
        if os.path.exists(path):
            stat_dfs.append(pd.read_table(path).astype(float))                 
    
    all_stats = pd.concat(stat_dfs, axis=1)
    all_stats = all_stats.T

    out_file = f"{out_prefix}_qc_stats.txt"

    all_stats.to_csv(out_file, sep="\t", header=False)    
    