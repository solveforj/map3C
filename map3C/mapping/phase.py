from bisect import bisect_left, bisect_right
from .utils import *

class ReadPhaser:

    def phase_read(self, read):
    
        phase = "."
        status = "N"
        eval_snps = 0
        a1 = 0
        a2 = 0

        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        
        if chrom not in self.variants:
            return phase, status, eval_snps, a1, a2 
    
        chrom_variants = self.variants[chrom]

        snp_s = bisect_left(chrom_variants, start, key=lambda x : x["pos"])
        snp_e = bisect_right(chrom_variants, end, key=lambda x : x["pos"])

        overlap_variants = chrom_variants[snp_s:snp_e]

        if len(overlap_variants) == 0:
            return phase, status, eval_snps, a1, a2 

        cigartuples = read.cigartuples
        seq = read.query_sequence
        qual = read.query_qualities
        qual_seq_agree = len(seq) == len(qual)
        
        x = start
        y = 0


        for cig in cigartuples:
            op = cig[0]
            length = cig[1]
            if op == 0: # M
                for i in range(len(overlap_variants)):
                    v = overlap_variants[i]
                    #print(v)
                    pos = v["pos"]
                    if x <= pos and pos < x + length:
                        pos = pos -x + y
                        nt = seq[pos]
                        #print(nt)
                        q = qual[pos] if qual_seq_agree else self.min_base_quality
                        if q >= self.min_base_quality:
                            if nt == v["a1"]:
                                a1 += 1
                            elif nt == v["a2"]:
                                a2 += 1
                            eval_snps += 1
                x += length
                y += length
            elif op == 1: # I
                y += length
            elif op == 2: # D
                x += length
            elif op == 4: # S
                y += length

        if a1 > 0 and a2 == 0:
            phase = "0"
            status = "P"
        elif a1 == 0 and a2 > 0:
            phase = "1"
            status = "P"
        elif a1 > 0 and a2 > 0:
            status = "C"
        
        return phase, status, eval_snps, a1, a2

    def __init__(self, 
                 variants,
                 min_base_quality=20       
                ):
        
        self.variants = process_variants(variants)
        self.conflicting_phase_count = 0
        self.successful_phase_count = 0
        self.na_phase_count = 0
            
        self.min_base_quality = min_base_quality

class ReadPhaserBisulfite(ReadPhaser):
        
    def phase_read(self, read):
        
        phase = "."
        status = "N"
        eval_snps = 0
        a1 = 0
        a2 = 0
        
        chrom = read.reference_name
        start = read.reference_start
        end = read.reference_end
        
        if chrom not in self.variants:
            return phase, status, eval_snps, a1, a2    
    
        chrom_variants = self.variants[chrom]

        snp_s = bisect_left(chrom_variants, start, key=lambda x : x["pos"])
        snp_e = bisect_right(chrom_variants, end, key=lambda x : x["pos"])

        overlap_variants = chrom_variants[snp_s:snp_e]

        if len(overlap_variants) == 0:
           return phase, status, eval_snps, a1, a2

        cigartuples = read.cigartuples
        seq = read.query_sequence
        qual = read.query_qualities
        qual_seq_agree = len(seq) == len(qual)
        ref_seq = read.get_reference_sequence().upper()
        direction = read.get_tag("YD")
        
        x = start
        y = 0

        rx = start
        ry = 0

        for cig in cigartuples:
            op = cig[0]
            length = cig[1]
            if op == 0: # M
                
                for i in range(len(overlap_variants)):
                    v = overlap_variants[i]
                    pos = v["pos"]
                    
                    if (x <= pos < x + length) and (rx <= pos < rx + length):
                        qpos = pos - x + y
                        qnt = seq[qpos]
                        
                        rpos = pos - rx + ry
                        rnt = ref_seq[rpos]

                        if direction == "f":
                            if rnt == "C" and qnt == "T":
                                continue
                        elif direction == "r":
                            if rnt == "G" and qnt == "A":
                                continue

                        q = qual[qpos] if qual_seq_agree else self.min_base_quality
                        if q >= self.min_base_quality:
                            if qnt == v["a1"]:
                                a1 += 1
                            elif qnt == v["a2"]:
                                a2 += 1
                            eval_snps += 1
                                
                x += length
                y += length
                rx += length
                ry += length
                
            elif op == 1: # I
                y += length
            elif op == 2: # D
                x += length
            elif op == 4: # S
                y += length

        if a1 > 0 and a2 == 0:
            phase = "0"
            status = "P"
        elif a1 == 0 and a2 > 0:
            phase = "1"
            status = "P"
        elif a1 > 0 and a2 > 0:
            status = "C"
        
        return phase, status, eval_snps, a1, a2
