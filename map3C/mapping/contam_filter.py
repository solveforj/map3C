import pysam 
import pandas as pd
from .utils import *

class ContaminationFilter:

    def compute_non_cg_methylation_star(self, read):

        m = 0
        um = 0
        
        pairs_info = read.get_aligned_pairs(with_seq=True)

        read_seq = read.query_sequence

        for i in range(len(pairs_info) - 1):
            pair = pairs_info[i]
            next_pair = pairs_info[i+1]

            ref_nuc = pair[-1]
            read_nuc_pos = pair[0]
            
            next_ref_nuc = next_pair[-1]
            
            if ref_nuc == None: 
                continue
            if next_ref_nuc == None:
                continue
            if read_nuc_pos == None:
                continue
            
            ref_nuc = ref_nuc.upper()
            next_ref_nuc = next_ref_nuc.upper()

            read_nuc = read_seq[read_nuc_pos].upper()

            if ref_nuc != "C":
                continue
            if next_ref_nuc == "G":
                continue

            if read_nuc == "C":
                m += 1
            elif read_nuc == "T":
                um += 1
        
        return m, um
        
    def compute_non_cg_methylation_biscuit(self, read):
        # Note: Only compatible with Biscuit
        
        m = 0
        um = 0

        zn = read.get_tag("ZN").split(",")
        zn_dict = {}
        for context in zn:
            context = context.split("_")
            context_type = context[0]
            if context_type != "CG":
                counts = context[1][1:].split("C")
                m += int(counts[0]) # retained
                um += int(counts[1]) # converted
    
        return m, um

    def evaluate_ch_v1(self, um, m):
        """
        Original map3C contamination removal, with some modifications relative to the original snm3c_bam_filter.py:
        1. There is a bug, where 1 is subtracted from methylated site count.
        2. mCH fraction should be <0.7, not <=0.7
        This was used in map3C versions up to and including 0.0.4.

        Rules:
        Discard all reads with same query name (i.e. primary/secondary alignments)
        if (>= 3 CH sites) OR (> 70% of CH sites are methylated). Originally performed at the mate, instead of pair, level
        """        
        
        m = m - 1

        if (um + m) < 3 or (m / (um+m)) <= 0.7:
            return True
        else:
            return False

    def evaluate_ch_v2(self, um, m):

        total_ch = um + m
        if total_ch == 0:
            m = m + 1e-10 # Add a pseudocount to prevent divide by zero error
            total_ch = um + m

        mCH_fraction = m / total_ch
        mCH_valid = self.min_mch_fraction <= mCH_fraction <= self.max_mch_fraction
        CH_valid = total_ch >= self.min_ch_sites
        CH_low = total_ch < self.min_ch_sites

        # Reads that fail mCH levels but have low number of CH sites are rescued
        if not mCH_valid and self.rescue_low_ch_sites and CH_low:
            return True

        if mCH_valid and CH_valid:
            return True
        
        return False

    def filter_reads_mate(self):

        for mate in ["1", "2"]:
            read_group_reads = self.read_group_reads[mate]
    
            if len(read_group_reads) == 0:
                return 
    
            um = self.counts[mate]["um"]
            m = self.counts[mate]["m"]
            
            if self.evaluate_ch(um, m):
                for r in read_group_reads:
                    self.bam_out.write(r)
                self.stats_dict[f"R{mate}_contam_pass"] += 1
            else:
                self.stats_dict[f"R{mate}_contam_fail"] += 1

    def filter_reads_pair(self):

        um = self.counts["1"]["um"] + self.counts["2"]["um"]
        m = self.counts["1"]["m"] + self.counts["2"]["m"]
        
        read_group_reads = self.read_group_reads["1"] + self.read_group_reads["2"]

        if len(read_group_reads) == 0:
            return 
        
        if self.evaluate_ch(um, m):
            for r in read_group_reads:
                self.bam_out.write(r)
            self.stats_dict[f"pairs_contam_pass"] += 1
        else:
            self.stats_dict[f"pairs_contam_fail"] += 1
        
    
    def filter_bam(self):

        """
        Mates are analyzed individually, since single cell prep creates chimeras,
        but mates could be analyzed jointly.
        """
        
        with pysam.AlignmentFile(self.bam) as bam_in, \
            pysam.AlignmentFile(self.out, "wb", header=bam_in.header) as self.bam_out:
            
            
            self.header = bam_in.header.to_dict()
            for i in self.header["PG"]:
                cl = i["CL"]
                if "biscuit" in cl:
                   self.compute_non_cg_methylation = self.compute_non_cg_methylation_biscuit
                   break
                if "STAR" in cl:
                    self.compute_non_cg_methylation = self.compute_non_cg_methylation_star
                    break
    
            read_group = None
            read_group_reads = None
            iter_count = 0

            self.counts = {"1" : {"um" : 0, "m" : 0},
                           "2" : {"um" : 0, "m" : 0}}

            
            
            for read in bam_in:
                
                read_name = self.get_read_id(read)
                mate = self.get_read_mate(read)
                
                
                if read.mapping_quality < self.min_mapq:
                      
                    # Do not count unmethylated/methylated sites for read with low MAPQ
                    r_m, r_um = 0, 0
                    
                else:
                    r_m, r_um = self.compute_non_cg_methylation(read)
        
                if iter_count == 0:
                    read_group = read_name
                    self.read_group_reads = {"1" : [], "2" : []}
                    iter_count = 1
                if read_name == read_group:
                    self.read_group_reads[mate].append(read)
                    self.counts[mate]["um"] += r_um
                    self.counts[mate]["m"] += r_m
                else:
                    
                    self.filter_reads()
                    
                    read_group = read_name
                    self.read_group_reads = {"1" : [], "2" : []}
                    self.read_group_reads[mate].append(read)
                    
                    self.counts = {"1" : {"um" : 0, "m" : 0},
                                   "2" : {"um" : 0, "m" : 0}}
                    
                    self.counts[mate]["um"] += r_um
                    self.counts[mate]["m"] += r_m

            if read_group == None:
                return
                
            self.filter_reads()

    
    def __init__(self, 
                 bam, 
                 out_prefix,
                 mate_annotation,
                 filter_type="mate",
                 min_mapq=30, 
                 min_mch_fraction=0.0,
                 max_mch_fraction=0.7,
                 min_ch_sites=3,
                 rescue_low_ch_sites=False,
                 old_contam_filter=False
                ):

        self.bam = bam
        self.min_mapq = min_mapq
        self.min_mch_fraction = min_mch_fraction
        self.max_mch_fraction = max_mch_fraction
        self.min_ch_sites = min_ch_sites
        self.rescue_low_ch_sites = rescue_low_ch_sites
        self.old_contam_filter = old_contam_filter
        self.out_prefix = out_prefix
        self.filter_type = filter_type

        self.out = f"{out_prefix}_contam_filtered.bam"
        self.stats = f"{out_prefix}_contam_stats.txt"
    
        if mate_annotation == "qname":
            self.get_read_id = lambda x: x.query_name.split("_")[0] 
            self.get_read_mate = lambda x: x.query_name.split("_")[1] 
        else:
            self.get_read_id = lambda x: x.query_name
            self.get_read_mate = get_mate_from_tag

        if self.old_contam_filter:
            self.evaluate_ch = self.evaluate_ch_v1
            self.filter_type = "mate"
        else:
            self.evaluate_ch = self.evaluate_ch_v2

        
        if self.filter_type == "mate":
            self.filter_reads = self.filter_reads_mate
            self.stats_dict = {"R1_contam_pass" : 0,
                               "R1_contam_fail" : 0,
                               "R2_contam_pass" : 0,
                               "R2_contam_fail": 0
                              }
        elif self.filter_type == "pair":
            self.filter_reads = self.filter_reads_pair
            self.stats_dict = {"pairs_contam_pass" : 0,
                               "pairs_contam_fail" : 0
                              }

            
        self.filter_bam()

        stats_df = pd.DataFrame.from_dict(self.stats_dict, orient="index").T
        stats_df.to_csv(self.stats, index=False, sep="\t")
