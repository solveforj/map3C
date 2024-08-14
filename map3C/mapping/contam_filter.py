import pysam 

class ContaminationFilter:
    
    def compute_non_cg_methylation(self, read):
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
    
    def filter_bam(self):

        """
        Mates are analyzed individually, since single cell prep creates chimeras,
        but mates could be analyzed jointly.
        """
        
        self.r1_kept = 0
        self.r1_removed = 0
        self.r2_kept = 0
        self.r2_removed = 0
        
        with pysam.AlignmentFile(self.bam) as bam_in, \
            pysam.AlignmentFile(self.out, "wb", header=bam_in.header) as bam_out:
            
            read_group = None
            read_group_reads = None
            iter_count = 0

            um = 0
            m = -1
            
            for read in bam_in:
                
                read_name = read.query_name
                mate = read_name.split("_")[-1]
                
                if read.mapping_quality < self.min_mapq:
                      
                    # Do not count unmethylated/methylated sites for read with low MAPQ
                    r_m, r_um = 0, 0
                    
                else:
                    r_m, r_um = self.compute_non_cg_methylation(read)
        
                if iter_count == 0:
                    read_group = read_name
                    read_group_reads = []
                    iter_count = 1
                if read_name == read_group:
                    read_group_reads.append(read)
                    um += r_um
                    m += r_m
                else:
                    """
                    Discard all reads with same query name (i.e. primary/secondary alignments)
                    if (>= 3 CH sites) OR (> 70% of CH sites are methylated)
                    """
                    if (um + m) < self.max_ch_sites or (m / (um + m)) <= self.max_mc_ch:
                        for r in read_group_reads:
                            bam_out.write(r)
                        if mate == "1":
                            self.r1_kept += 1
                        else:
                            self.r2_kept += 1
                    else:
                        if mate == "1":
                            self.r1_removed += 1
                        else:
                            self.r2_removed += 1
                    read_group = read_name
                    read_group_reads = [read]
                    um = 0
                    m = -1
                    um += r_um
                    m += r_m
            if read_group == None:
                return
            if (um + m) < self.max_ch_sites or (m / (um + m)) <= self.max_mc_ch:
                for r in read_group_reads:
                    bam_out.write(r)
                if mate == "1":
                    self.r1_kept += 1
                else:
                    self.r2_kept += 1
            else:
                if mate == "1":
                    self.r1_removed += 1
                else:
                    self.r2_removed += 1          

    
    def __init__(self, 
                 bam, 
                 min_mapq, 
                 max_mc_ch,
                 max_ch_sites,
                 out_prefix):

        self.bam = bam
        self.min_mapq = min_mapq
        self.max_mc_ch = max_mc_ch
        self.max_ch_sites = max_ch_sites
        self.out_prefix = out_prefix

        self.out = f"{out_prefix}_contam_filtered.bam"
        self.stats = f"{out_prefix}_contam_stats.txt"

        self.filter_bam()

        with open(self.stats, "w") as f:
            f.write("\t".join(["R1_contam_pass", "R1_contam_fail", 
                               "R2_contam_pass", "R2_contam_fail"]) + "\n")
            f.write("\t".join([str(i) for i in [
                self.r1_kept, self.r1_removed, 
                self.r2_kept, self.r2_removed]]) + "\n")  
