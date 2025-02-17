import pysam 
import pandas as pd
from .utils import *

class ContaminationFilter:
    
    def compute_non_cg_methylation(self, read):
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

    def filter_reads(self, mate):
        
        read_group_reads = self.read_group_reads[mate]

        if len(read_group_reads) == 0:
            return 

        um = self.counts[mate]["um"]
        m = self.counts[mate]["m"]
        #print(mate, m, um, (m / (um + m)))
        
        if (um + m) < self.max_ch_sites or (m / (um + m)) <= self.max_mc_ch:
            for r in read_group_reads:
                self.bam_out.write(r)
            self.stats_dict[f"R{mate}_contam_pass"] += 1
        else:
            self.stats_dict[f"R{mate}_contam_fail"] += 1
        
    
    def filter_bam(self):

        """
        Mates are analyzed individually, since single cell prep creates chimeras,
        but mates could be analyzed jointly.
        """
        
        with pysam.AlignmentFile(self.bam) as bam_in, \
            pysam.AlignmentFile(self.out, "wb", header=bam_in.header) as self.bam_out:
            
            read_group = None
            read_group_reads = None
            iter_count = 0

            self.counts = {"1" : {"um" : 0, "m" : -1},
                           "2" : {"um" : 0, "m" : -1}}
            
            for read in bam_in:
                
                #read_name = read.query_name
                #mate = read_name.split("_")[-1]

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
                    """
                    Discard all reads with same query name (i.e. primary/secondary alignments)
                    if (>= 3 CH sites) OR (> 70% of CH sites are methylated)
                    """
                    self.filter_reads("1")
                    self.filter_reads("2")
                    
                    read_group = read_name
                    self.read_group_reads = {"1" : [], "2" : []}
                    self.read_group_reads[mate].append(read)
                    
                    self.counts = {"1" : {"um" : 0, "m" : -1},
                                   "2" : {"um" : 0, "m" : -1}}
                    
                    self.counts[mate]["um"] += r_um
                    self.counts[mate]["m"] += r_m

            if read_group == None:
                return
                
            self.filter_reads("1")
            self.filter_reads("2")        

    
    def __init__(self, 
                 bam, 
                 out_prefix,
                 mate_annotation,
                 min_mapq, 
                 max_mc_ch,
                 max_ch_sites):

        self.bam = bam
        self.min_mapq = min_mapq
        self.max_mc_ch = max_mc_ch
        self.max_ch_sites = max_ch_sites
        self.out_prefix = out_prefix

        self.out = f"{out_prefix}_contam_filtered.bam"
        self.stats = f"{out_prefix}_contam_stats.txt"

        self.stats_dict = {"R1_contam_pass" : 0,
                           "R1_contam_fail" : 0,
                           "R2_contam_pass" : 0,
                           "R2_contam_fail": 0
                          }
    
        if mate_annotation == "qname":
            self.get_read_id = lambda x: x.query_name.split("_")[0] 
            self.get_read_mate = lambda x: x.query_name.split("_")[1] 
        else:
            self.get_read_id = lambda x: x.query_name
            self.get_read_mate = get_mate_from_tag
            
        self.filter_bam()

        stats_df = pd.DataFrame.from_dict(self.stats_dict, orient="index").T
        stats_df.to_csv(self.stats, index=False, sep="\t")
