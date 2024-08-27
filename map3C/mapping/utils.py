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
            
    
def process_restriction_sites(restriction_sites, restriction_enzymes):

    restriction_sites_dict = {}
    restriction_sites = sum(restriction_sites, [])
    restriction_enzymes = sum(restriction_enzymes, [])
    print(restriction_sites)
    print(restriction_enzymes)
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

    def parse_base(self, line):
        chrom1 = line[1]
        chrom2 = line[3]
        pos1 = line[2]
        pos2 = line[4]

        if chrom1 == chrom2:
            dist = np.abs(int(pos2) - int(pos1))
            if dist >= 1000:
                self.pairs_stats[f"{self.prefix}_intra1kb"] += 1
            if dist >= 10000:
                self.pairs_stats[f"{self.prefix}_intra10kb"] += 1
            if dist >= 20000:
                self.pairs_stats[f"{self.prefix}_intra20kb"] += 1
        else:
            self.pairs_stats[f"{self.prefix}_inter"] += 1

        self.pairs_stats[f"{self.prefix}_total"] += 1


    def parse_phase(self, line):
        phase1 = line[self.phase1_idx]
        phase2 = line[self.phase2_idx]

        phase_id = "".join(sorted([phase1, phase2]))

        self.pairs_stats[f"{self.prefix}_phased_{phase_id}"] += 1

    def parse_metadata(self, line):
        pair_type = line[self.pair_type_idx]
        if pair_type == "RU":
            pair_type = "UR"
        rule = line[self.rule_idx]
        reads = line[self.reads_idx]

        pair_type_rule = f"{pair_type}_{rule}"

        self.pairs_stats[f"{self.prefix}_{reads}"] += 1
        self.pairs_stats[f"{self.prefix}_{pair_type_rule}"] += 1
        
    def parse_file(self):
        self.parsers = [self.parse_base]
        
        with gzip.open(self.pair_file, "rt") as f:
            for line in f:
                line = line.strip()
                if line.startswith("#columns"):
                    columns = line[10:].split()

                    self.pair_type_idx = columns.index("pair_type")
                    
                    if "phase0" in columns and "phase1" in columns:
                        self.phase1_idx = columns.index("phase0")
                        self.phase2_idx = columns.index("phase1")
                        self.pairs_stats.update({
                            f"{self.prefix}_phased_.0" : 0,
                            f"{self.prefix}_phased_.1" : 0,
                            f"{self.prefix}_phased_.." : 0,
                            f"{self.prefix}_phased_11" : 0,
                            f"{self.prefix}_phased_00" : 0,
                            f"{self.prefix}_phased_01" : 0,
                        })
                        
                        self.parsers.append(self.parse_phase)
                        
                    if "rule" in columns:
                        self.rule_idx = columns.index("rule")
                        self.reads_idx = columns.index("reads")
                        #contact_class_idx = columns.index("contact_class")
                        #multimap_overlap_idx = columns.index("multimap_overlap")
                        #cs_locs_idx = columns.index("cut_site_locs")
                        
                        self.pairs_stats.update({
                            f"{self.prefix}_R1" : 0,
                            f"{self.prefix}_R2" : 0,
                            f"{self.prefix}_R1_rescue" : 0,
                            f"{self.prefix}_R2_rescue" : 0,
                            f"{self.prefix}_R1-2" : 0,
                            f"{self.prefix}_R1&2" : 0,
                            f"{self.prefix}_UU_all" : 0,
                            f"{self.prefix}_UU_mask" : 0,
                            f"{self.prefix}_UR_mask" : 0
                        })
                        
                        self.parsers.append(self.parse_metadata)
                        
                elif line.startswith("#"):
                    continue

                else:
                    line = line.split()
                    for parser in self.parsers:
                        parser(line)
                    
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

        self.pairs_stats.update({f"{self.prefix}_dup_rate":pairs_dup_rate})
                        
    def __init__(self, pair_file, prefix, pairtools_dedup_file=None):
        
        self.pair_file = pair_file
        self.prefix = prefix
        
        self.pairs_stats = {
            f"{self.prefix}_total" : 0,
            f"{self.prefix}_intra1kb" : 0,
            f"{self.prefix}_intra10kb" : 0,
            f"{self.prefix}_intra20kb" : 0,
            f"{self.prefix}_inter" : 0
        }
        
        self.parse_file()

        if pairtools_dedup_file:
            self.pairtools_dedup_file = pairtools_dedup_file
            self.parse_dedup_file()


def pairtools_stats(out_prefix,
                    contacts,
                    artefacts,
                    contacts_dedup_stats=None,
                    artefacts_dedup_stats=None,
                    filterbycov_stats=None
                   ):

    stats_path = f"{out_prefix}_pairs_stats.txt"
    
    contacts_stats = PairsStats(contacts, "contacts", contacts_dedup_stats).pairs_stats

    artefacts_stats = PairsStats(artefacts, "artefacts", artefacts_dedup_stats).pairs_stats

    if filterbycov_stats:
        if os.path.isfile(filterbycov_stats):
            with open(filterbycov_stats) as f:
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
    full_stats.update(artefacts_stats)
    
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
                       mode,
                       min_mapq = 30, 
                       min_base_quality = 20):

    if mode == "bsdna":
        
        txt_paths = [f"{out_prefix}_trim_stats.txt",
                     f"{out_prefix}_contam_stats.txt",
                     f"{out_prefix}_dupsifter_stats.txt", 
                     f"{out_prefix}_alignment_stats.txt",
                     f"{out_prefix}_pairs_stats.txt",
                     f"{out_prefix}.allc.tsv.gz_methylation_stats.txt"
                    ]

    elif mode == "dna":
        
        txt_paths = [f"{out_prefix}_trim_stats.txt",
                     f"{out_prefix}_alignment_stats.txt",
                     f"{out_prefix}_pairs_stats.txt",
                    ]
    
    stat_dfs = [pd.DataFrame([job], columns=["job"])]
    for path in txt_paths:
        if os.path.exists(path):
            stat_dfs.append(pd.read_table(path).astype(float))                 
    
    cov_map_stats = []
    cov_map_columns = []
    trimmed_bam = f"{out_prefix}_trimmed_sorted.bam"
    if os.path.exists(trimmed_bam):
        genome_cov_dedup_trim = compute_genome_coverage(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dedup_trim = compute_mapped_nucleotides(trimmed_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dedup_trim, mapped_nuc_dedup_trim]
        cov_map_columns += ["reference_coverage_dedup_trim", "mapped_nucleotides_dedup_trim"]

    mkdup_bam = f"{out_prefix}_mkdup_sorted.bam"
    if os.path.exists(mkdup_bam):
        genome_cov_dup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
        genome_cov_dedup = compute_genome_coverage(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=True)
        mapped_nuc_dedup = compute_mapped_nucleotides(mkdup_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dup, genome_cov_dedup, mapped_nuc_dup, mapped_nuc_dedup]
        cov_map_columns += ["reference_coverage_dup", "reference_coverage_dedup", "mapped_nucleotides_dup", "mapped_nucleotides_dedup"]

    masked_bam = f"{out_prefix}_masked_sorted.bam"
    if os.path.exists(masked_bam):
        genome_cov_dedup_mask = compute_genome_coverage(masked_bam, min_mapq, min_base_quality, keep_dup=False)
        mapped_nuc_dedup_mask = compute_mapped_nucleotides(masked_bam, min_mapq, min_base_quality, keep_dup=False)
        cov_map_stats += [genome_cov_dedup_mask, mapped_nuc_dedup_mask]
        cov_map_columns += ["reference_coverage_dedup_mask", "mapped_nucleotides_dedup_mask"]

    cov_map_df = pd.DataFrame([cov_map_stats], columns=cov_map_columns)
    stat_dfs.append(cov_map_df)

    all_stats = pd.concat(stat_dfs, axis=1)
    all_stats = all_stats.T

    out_file = f"{out_prefix}_qc_stats.txt"

    all_stats.to_csv(out_file, sep="\t", header=False)    
    