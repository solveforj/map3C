import pysam
import numpy as np
from collections import OrderedDict

rng = np.random.default_rng(1)

def mask_qualities(read, start, end):

    original_qualities = pysam.qualities_to_qualitystring(read.query_qualities)

    original_qualities = [*original_qualities]

    full_cigar = []
    for i in read.cigartuples:
        if i[0] != 5: # Do not include hard clipping (on secondary alignments)
            full_cigar += [i[0]] * i[1]

    qual_idx = 0
    algn_idx = read.reference_start
    mask = False
    for i in range(len(full_cigar)):
        cigar_val = full_cigar[i]
        #print(cigar_val, qual_idx, algn_idx)
        if cigar_val in [0, 1, 4]:
            if cigar_val == 0:
                if not mask:
                    mask = True
                if mask:
                    if algn_idx >= start:
                        original_qualities[qual_idx] = "!"
                    algn_idx += 1
            qual_idx += 1
        elif cigar_val == 2:
            algn_idx += 1
        if algn_idx == end:
            break

    original_qualities = "".join(original_qualities)
    
    read.query_qualities = pysam.qualitystring_to_array(original_qualities)
    return read

def get_mate_from_tag(read):
    if read.is_read1:
        return "1"
    elif read.is_read2:
        return "2"
    else:
        raise Exception(f"Mate not defined for read {read.query_name}")

# Masked R1 and R2 as input
def mask_overlaps(R1, R2):
    all_reads = []
    for read in R1:
        all_reads.append(R1[read])
    for read in R2:
        all_reads.append(R2[read])

    rng.shuffle(all_reads)
    all_reads = sorted(all_reads, key = lambda x : (x['trimmed_read'].mapping_quality, 
                                                    len(x['trimmed_read'].query_alignment_sequence)), reverse=True)

    for i in range(len(all_reads) - 1):
        read1 = all_reads[i]["trimmed_read"]
        mate1 = all_reads[i]["mate"]
        for j in range(i + 1, len(all_reads)):
            read2 = all_reads[j]["trimmed_read"]
            mate2 = all_reads[j]["mate"]
            
            if mate1 == mate2:
                continue
            if read1.reference_name != read2.reference_name:
                continue

            overlap_start = max(read1.reference_start, read2.reference_start)
            overlap_end = min(read1.reference_end, read2.reference_end)

            if overlap_end <= overlap_start:
                continue
            
            read2 = mask_qualities(
              read2, 
              overlap_start, 
              overlap_end
            )

            all_reads[j]["trimmed_read"] = read2

    new_R1 = []
    new_R2 = []

    for read in all_reads:
        mate = read["mate"]
        if mate == "1":
            new_R1.append(read["trimmed_read"])
        elif mate == "2":
            new_R2.append(read["trimmed_read"])
    return new_R1, new_R2

class OverlapMask:
    
    def process_mate(self, all_alignments, mate, read_group_name):

        ordered_reads = OrderedDict()

            
        for i in range(len(all_alignments)):
            read = all_alignments[i]
            ordered_reads[i] = {"trimmed_read" : read,
                                "mate" : mate}
                        
        return ordered_reads
    
    def process_read_group(self, read_group, read_group_name, bam_out):

        if read_group == None:
            return
        r1 = []
        
        r2 = []

        # Divide reads into R1 and R2
        for read in read_group:
            # Handles low-MAPQ primary alignments for chimeric reads
            if read.mapping_quality < self.min_mapq:
                continue
            if self.get_read_mate(read) == "1":
                r1.append(read)
            elif self.get_read_mate(read) == "2":
                r2.append(read)

        # Order reads from 5' to 3'
        R1 = self.process_mate(r1, "1", read_group_name)
        R2 = self.process_mate(r2, "2", read_group_name)
        
        masked_R1, masked_R2 = mask_overlaps(R1, R2)
        
        for read in masked_R1:
            bam_out.write(read)

        for read in masked_R2:
            bam_out.write(read)

    def process_bam(self):

        iter_count = 0
        read_group = None
        read_group_name = None
        with pysam.AlignmentFile(self.trimmed_bam, index_filename=None) as bam_in:
            with pysam.AlignmentFile(self.masked_bam, 'wb', template=bam_in) as bam_out:
                                
                for read in bam_in:
                    read_id = self.get_read_id(read)
                    #print(self.get_read_mate(read))
                    #print(get_mate_from_tag(read))
                    if iter_count == 0:
                        read_group_name = read_id
                        read_group = []
                        iter_count = 1
    
                    if read_group_name == read_id:
                        read_group.append(read)
                    else:
                        self.process_read_group(read_group, read_group_name, bam_out)
                        
                        read_group_name = read_id
                        read_group = [read]
                self.process_read_group(read_group, read_group_name, bam_out)
        
    def __init__(self, bam, out_prefix, min_mapq=30, manual_mate_annotation=False):

        self.min_mapq = min_mapq
        self.trimmed_bam = bam
        self.masked_bam = f'{out_prefix}_masked.bam'
        
        if manual_mate_annotation:
            self.get_read_id = lambda x: x.query_name.split("_")[0] 
            self.get_read_mate = lambda x: x.query_name.split("_")[1] 
        else:
            self.get_read_id = lambda x: x.query_name
            self.get_read_mate = get_mate_from_tag

        self.process_bam()
        