import re
import pysam
import numpy as np
import bisect
from collections import OrderedDict

rng = np.random.default_rng(1)

def get_biscuit_tags(pairs, direction):

    ZC = 0
    ZR = 0
    NM = 0
    possible_conv = 0
    
    for pair in pairs:
        if pair["cigar"] == 0:
            if pair["ref_nuc"].islower():
                NM += 1
                if pair["ref_nuc"] == "c" and pair["read_nuc"] == "T":
                    possible_conv += 1
                elif pair["ref_nuc"] == "g" and pair["read_nuc"] == "A":
                    possible_conv += 1
            if direction == "f":
                if pair["ref_nuc"] == "c" and pair["read_nuc"] == "T":
                    ZC += 1
                elif pair["ref_nuc"] == "C" and pair["read_nuc"] == "C":
                    ZR += 1
            elif direction == "r":
                if pair["ref_nuc"] == "g" and pair["read_nuc"] == "A":
                    ZC += 1
                elif pair["ref_nuc"] == "G" and pair["read_nuc"] == "G":
                    ZR += 1
        elif pair["cigar"] in [1, 2]:
            NM += 1
    
    NM = NM - possible_conv

    return ZC, ZR, NM
    
def get_bwa_tags(pairs):
    NM = 0
    
    for pair in pairs:
        if pair["cigar"] == 0:
            if pair["ref_nuc"].islower():
                NM += 1
        elif pair["cigar"] in [1, 2]:
            NM += 1

    return NM

def get_md_tag(pairs):
    MD = ""
    match_count = 0
    deletion = False

    for i in range(len(pairs)):
        pair = pairs[i]
        cigar_id = pair["cigar"]
        if cigar_id == 0:
            deletion = False
            if pair["ref_nuc"].isupper():
                match_count += 1
            else:
                MD += str(match_count)
                MD += pair["ref_nuc"].upper()
                match_count = 0
        elif cigar_id == 2:
            if deletion == False:
                deletion = True
                MD += str(match_count)
                MD += "^"
                match_count = 0
            MD += pair["ref_nuc"].upper()
    MD += str(match_count)

    return MD

def softclip_pairs(pairs):
    new_pairs = []
    for i in range(len(pairs)):
        if pairs[i]["cigar"] in [0,1]:
            new_algn = {
                "cigar" : 4,
                "read_pos" : None,
                "ref_pos" : None,
                "ref_nuc" : None,
                "read_nuc" : None,
                "read_qual" : None,
            }
            new_pairs.append(new_algn)
        elif pairs[i]["cigar"] in [2]:
            continue
        else:
            new_pairs.append(pairs[i])
        
    return new_pairs            

def reverse_complement(seq):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return "".join(nn[n] for n in reversed(seq))

def complement(nuc):
    nn = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'}
    return nn[nuc.upper()]

def alignment_info(read, seg, original_sequence):
    pairs = read.get_aligned_pairs(with_seq = True)
    cigar = read.cigartuples
    read_seq = read.query_sequence
    read_qualities = pysam.qualities_to_qualitystring(read.query_qualities)
    
    if read.is_forward:
        read_start = seg[0]
        increment = 1
    else:
        read_start = seg[1] - 1 # Due to indexing
        increment = -1
    
    full_cigar = []
    for i in cigar:
        if i[0] != 5: # Do not include hard clipping (on secondary alignments)
            full_cigar += [i[0]] * i[1]

    new_pairs = []

    seq_idx = 0
    for i in range(len(pairs)):
        pair = pairs[i]
        cigar_val = full_cigar[i]
        algn_pos = pair[0]
        ref_pos = pair[1]
        ref_nuc = pair[2]
        read_nuc = None
        qual_val = None
        read_pos = None
        original_read_nuc = None
        
        if cigar_val in [0,1, 4]:
            read_nuc = read_seq[seq_idx]
            qual_val = read_qualities[seq_idx]
            if cigar_val != 4:
                read_pos = read_start 
                #original_read_nuc = original_sequence[read_pos]
                read_start += increment
            seq_idx += 1            

        pair_dict = {
            "cigar":cigar_val,
            "algn_pos": algn_pos,
            "read_pos": read_pos,
            "ref_pos": ref_pos,
            "ref_nuc": ref_nuc,
            "read_nuc" : read_nuc,
            "read_qual" : qual_val
        }
        new_pairs.append(pair_dict)

    return new_pairs

def overlap_reference_alignment(aligned_pairs, overlap_start, overlap_end, is_forward):
    seq = ""
    collect = False
    if is_forward:
        iterator = range(len(aligned_pairs))
    else:
        iterator = reversed(range(len(aligned_pairs)))

    for i in iterator:
        pair = aligned_pairs[i]
        if pair["read_pos"] == overlap_start:
            collect = True
        if collect:
            nuc = pair["ref_nuc"]
            if pair['cigar'] == 1:
                nuc = "-" # Add gap for insertion
            elif not is_forward:
                nuc = complement(nuc)
            seq += nuc.upper()
            
        if pair["read_pos"] == overlap_end-1:
            break

    return seq

def get_cut_site_spans(site, seq):
    spans = [i.span() for i in re.finditer(site, seq)]
    return spans

def get_5_adj_point(pairs, original_point, is_forward):
    # For a given aligned segment of a read, adjust 5' end of segment in case given end is not mapped
    if is_forward:
        index_5 = 0
        iterator = range(len(pairs))
    else:
        index_5 = 1 # Index will be at end of range
        iterator = reversed(range(len(pairs)))
        
    for i in iterator:
        pair = pairs[i]
        if pair["cigar"] == 0:
            # Adjusted pos should be equal to the original point or closer to the 3' end
            if pair["read_pos"] >= original_point: 
                pos_5 = pair["ref_pos"]
                read_5 = pair["read_pos"]
                index_5 += i
                break
    return pos_5, index_5, read_5

def get_3_adj_point(pairs, original_point, is_forward):
    # For a given aligned segment of a read, adjust 3' end of segment in case given end is not mapped
    if is_forward:
        index_3 = 1 # Index will be at end of range
        iterator = reversed(range(len(pairs)))
    else:
        index_3 = 0
        iterator = range(len(pairs))
        
    for i in iterator:
        pair = pairs[i]
        if pair["cigar"] == 0:
            # Adjusted pos should be equal to the original point or closer to the 5' end
            if pair["read_pos"] <= original_point:
                pos_3 = pair["ref_pos"]
                read_3 = pair["read_pos"]
                index_3 += i
                break
    return pos_3, index_3, read_3

def build_cigar(read, pairs, start_index, end_index):
    full_cigar = []

    if not read.is_secondary:
      mask_val = 4
    else:
      mask_val = 5
        
    mask = True
    for i in range(len(pairs)):
        pair = pairs[i]
        cigar_val = pair["cigar"]
        if i < start_index or i >= end_index:
            if cigar_val in [0, 1]:
                cigar_val = mask_val
            elif cigar_val in [2]:
                cigar_val = "*"
        full_cigar.append(cigar_val)

    if read.is_secondary:
        expanded_left = []
        left_original = read.cigartuples[0]
        if left_original[0] == 5:
            expanded_left = [5] * left_original[1]

        expanded_right = []
        right_original = read.cigartuples[-1]
        if right_original[0] == 5:
            expanded_right = [5] * right_original[1]

        full_cigar = expanded_left + full_cigar + expanded_right

    
    updated_cigartuples = []

    full_cigar = [i for i in full_cigar if i != "*"]

    counter = 1
    for i in range(1, len(full_cigar)):
        val = full_cigar[i]
        previous_val = full_cigar[i-1]
        if val != previous_val:
            updated_cigartuples.append((previous_val, counter))
            counter = 0
        counter += 1
    updated_cigartuples.append((val, counter))

    return updated_cigartuples

def get_ref_pos(pairs, new_pos, is_forward):
    # Determines position of pairs where specific reference position is located
    
    # Handles issues where new position is not mapped
    if is_forward:
        index_5 = 0
        iterator = range(len(pairs))
        for i in iterator:
            pair = pairs[i]
            if pair["cigar"] == 0:
                # Adjusted pos should be equal to the new position or closer to the 3' end
                if pair["ref_pos"] >= new_pos: 
                    pos_5 = pair["ref_pos"]
                    read_5 = pair["read_pos"]
                    index_5 += i
                    break
    else:
        index_5 = 1 # Index will be at end of range
        iterator = reversed(range(len(pairs)))
        for i in iterator:
            pair = pairs[i]
            if pair["cigar"] == 0:
                # Adjusted pos should be equal to the new position or closer to the 3' end
                if pair["ref_pos"] <= new_pos: 
                    pos_5 = pair["ref_pos"]
                    read_5 = pair["read_pos"]
                    index_5 += i
                    break

    return pos_5, read_5

def adjust_read5_cut(pairs5, read5, r5_closest, seg5):
    # For the 5' alignment in sequential alignments on the same read, trim alignment range relative to  
    # original read to remove cut site
    if read5.is_forward:
        if r5_closest < read5.reference_end:
            trim_stop = r5_closest - 1
            _, adj_trim_stop = get_ref_pos(pairs5, trim_stop, read5.is_forward)
            seg5 = (seg5[0], adj_trim_stop + 1, seg5[2])
    else:
        if r5_closest > read5.reference_start - 4:
            trim_stop = r5_closest + 4
            _, adj_trim_stop = get_ref_pos(pairs5, trim_stop, read5.is_forward)
            seg5 = (seg5[0], adj_trim_stop + 1, seg5[2])

    return seg5

def adjust_read3_cut(pairs3, read3, r3_closest, seg3):
    # For the 3' alignment in sequential alignments on the same read, trim alignment range relative to 
    # original read to remove cut site
    if read3.is_forward:
        if r3_closest >= read3.reference_start:
            trim_stop = r3_closest + 4
            _, adj_trim_stop = get_ref_pos(pairs3, trim_stop, read3.is_forward)
            seg3 = (adj_trim_stop, seg3[1], seg3[2])
    else:
        if r3_closest <= read3.reference_end - 4:
            trim_stop = r3_closest - 1
            _, adj_trim_stop = get_ref_pos(pairs3, trim_stop, read3.is_forward)
            seg3 = (adj_trim_stop, seg3[1], seg3[2])
            
    return seg3

def get_cut_site_spans(site, seq):
    spans = [i.span() for i in re.finditer(site, seq)]
    return spans


def closest_restriction_site(chrom, pos, restriction_sites_dict, rule="closest"):

    results = {}
    
    for enzyme in restriction_sites_dict:
        results[enzyme] = {}
        
        restriction_sites_chrom = restriction_sites_dict[enzyme][chrom]

        rs_chrom_len = len(restriction_sites_chrom)
        
        insert = bisect.bisect_left(restriction_sites_chrom, pos)

        if insert <= 0:
            fragment = f"{chrom}_{insert}_{insert}"
            
            upstream_pos = restriction_sites_chrom[insert] - 1
            downstream_pos = restriction_sites_chrom[insert] - 1

        elif insert >= rs_chrom_len:
            fragment = f"{chrom}_{insert-1}_{insert-1}"
            
            upstream_pos = restriction_sites_chrom[insert-1] - 1
            downstream_pos = restriction_sites_chrom[insert-1] - 1

        else:
            fragment = f"{chrom}_{insert-1}_{insert}"
    
            upstream_pos = restriction_sites_chrom[insert-1] - 1
            downstream_pos = restriction_sites_chrom[insert] - 1
    
        upstream_dist = np.abs(upstream_pos - pos)
        downstream_dist = np.abs(downstream_pos - pos)

        if rule == "closest":
            if upstream_dist < downstream_dist:
                dist = upstream_dist
                site_pos = upstream_pos
            else:
                dist = downstream_dist
                site_pos = downstream_pos
        elif rule == "upstream":
            dist = upstream_dist
            site_pos = upstream_pos
        elif rule == "downstream":
            dist = downstream_dist
            site_pos = downstream_pos

        results[enzyme] = {"site_pos":site_pos, "dist":dist, "fragment": fragment}

    return results

def classify_breakpoint(r5_site_info, r3_site_info, max_cut_site_distance):

    contact_classes = []
    total_trim = {}
    contact_enzyme = None

    for enzyme in r5_site_info:
        r5_dist = r5_site_info[enzyme]["dist"]
        r3_dist = r3_site_info[enzyme]["dist"]
        r5_less = r5_dist <= max_cut_site_distance
        r3_less = r3_dist <= max_cut_site_distance

        if r5_less and r3_less:
            contact_classes.append(f"{enzyme}_both")
            total_trim[enzyme] = r5_site_info[enzyme]["dist"] + r5_site_info[enzyme]["dist"] 
        elif r5_less:
            contact_classes.append(f"{enzyme}_5")
        elif r3_less:
            contact_classes.append(f"{enzyme}_3")

    # Assign enzyme to valid contacts
    # If there are multiple enzymes assigned to a valid contact, one will be chosen which will require the least trimming
    if len(total_trim) == 0:
        contact_enzyme = "artefact"
    elif len(total_trim) == 1:
        contact_enzyme = list(total_trim.keys())[0]
    elif len(total_trim) > 1:
        min_trim = min(total_trim.values())
        possible_enzymes = [k for k, v in total_trim.items() if v==min_trim]
        if len(possible_enzymes) == 1:
            contact_enzyme = possible_enzymes[0]
        elif len(possible_enzymes) > 1:
            contact_enzyme = rng.choice(possible_enzymes)

    if len(contact_classes) == 0:
        contact_classes = ["na"]

    return contact_classes, contact_enzyme

def split_pair_to_restriction_site(read5, read3, restriction_sites, max_cut_site_distance):

    if read5.is_forward:
        r5_rs = closest_restriction_site(read5.reference_name, 
                                read5.reference_end,
                                restriction_sites)
    else:
        r5_rs = closest_restriction_site(read5.reference_name, 
                        read5.reference_start,
                        restriction_sites) 

    if read3.is_forward:
        r3_rs = closest_restriction_site(read3.reference_name, 
                                        read3.reference_start,
                                        restriction_sites)
    else:
        r3_rs = closest_restriction_site(read3.reference_name, 
                                        read3.reference_end,
                                        restriction_sites)

    bp_class, bp_enzyme = classify_breakpoint(r5_rs, r3_rs, max_cut_site_distance)

    return bp_class, bp_enzyme, r5_rs, r3_rs

def gap_pair_to_restriction_site(read5, read3, restriction_sites, max_cut_site_distance):

    if read5.is_forward:
        r5_rule = "downstream"
        r5_pos = read5.reference_end
    else:
        r5_rule = "upstream"
        r5_pos = read5.reference_start

    r5_rs = closest_restriction_site(read5.reference_name, 
                                     r5_pos, restriction_sites, rule=r5_rule) 

    if read3.is_forward:        
        r3_rule = "downstream"
        r3_pos = read3.reference_end
    else:        
        r3_rule = "upstream"
        r3_pos = read3.reference_start
        
    r3_rs = closest_restriction_site(read3.reference_name, 
                                     r3_pos, restriction_sites, rule=r3_rule) 
    
    bp_class, bp_enzyme = classify_breakpoint(r5_rs, r3_rs, max_cut_site_distance)

    return bp_class, bp_enzyme, r5_rs, r3_rs

class ReadTrimmer:

    def adjust_split_enzyme(self):
        
        if self.overlap > 0:
            self.first_trim += 1
        
        if self.bp_enzyme != "artefact":
            r5_closest = self.r5_rs[self.bp_enzyme]["site_pos"]
            r3_closest = self.r3_rs[self.bp_enzyme]["site_pos"]
            seg5 = adjust_read5_cut(self.pairs5, self.read5, r5_closest, self.seg5)
            seg3 = adjust_read3_cut(self.pairs3, self.read3, r3_closest, self.seg3)
            self.overlap = seg5[1] - seg3[0]
            if self.overlap > 0:
                self.second_trim += 1
            self.cut_site_pairs += 1
            return seg5, seg3
        else:
            return self.seg5, self.seg3

    def adjust_split_win(self):

        seg5, seg3 = self.adjust_split_enzyme()
        
        if self.overlap > 0:
            # Handle overlap even if both reads had cut site trimming or neither reads had cut site
            if self.read5.mapping_quality > self.read3.mapping_quality:
                r5_cover = True
            elif self.read5.mapping_quality == self.read3.mapping_quality:
                r5_len = len(self.read5.query_alignment_sequence)
                r3_len = len(self.read3.query_alignment_sequence)
                if r5_len > r3_len:
                    r5_cover = True
                elif r5_len == r3_len:
                    r5_cover = rng.choice([True, False])
                else:
                    r5_cover = False
            else:
                r5_cover = False
        
            if r5_cover:
                # Trim 3' read
                seg3 = (seg5[1], seg3[1], seg3[2])
            else:
                seg5 = (seg5[0], seg3[0], seg5[2])

        return seg5, seg3

    def adjust_split_complete(self):

        seg5, seg3 = self.adjust_split_enzyme()
        
        if self.overlap > 0:
            # Handle overlap even if both reads had cut site trimming or neither reads had cut site
            min_idx = min(seg5[1], seg3[0])
            max_idx = max(seg5[1], seg3[0])
        
            seg5 = (seg5[0], min_idx, seg5[2])
            seg3 = (max_idx, seg3[1], seg3[2])

        return seg5, seg3

    def adjust_split_none(self):

        return self.seg5, self.seg3
  
    def create_trimmed_mate(self, read, original_span, adjusted_span, pairs):

        #mate = self.mate
        header = self.header
        full_bam = self.full_bam
        
        pos_5, index_5, read_5 = get_5_adj_point(pairs, adjusted_span[0], read.is_forward)
        # Index is at end of range; subtract 1 to get actual final index
        pos_3, index_3, read_3 = get_3_adj_point(pairs, adjusted_span[1] - 1, read.is_forward)
    
        if read.is_forward:
            start_index = index_5
            end_index = index_3
            pos = pos_5
        else:
            start_index = index_3
            end_index = index_5
            pos = pos_3
    
        if start_index >= end_index:
            return None, None
        
        new_pairs = pairs[start_index:end_index]
        if not read.is_secondary:
            left_pairs = softclip_pairs(pairs[:start_index])
            right_pairs = softclip_pairs(pairs[end_index:])
            new_pairs = left_pairs + new_pairs + right_pairs
    
        if read.is_secondary:
            read_seq = ""
            read_qual = ""
            for pair in new_pairs:
                if pair["cigar"] in [0, 1]:
                    read_seq += pair["read_nuc"]
                    read_qual += pair["read_qual"]
        else:
            read_seq = read.query_sequence
            read_qual = pysam.qualities_to_qualitystring(read.query_qualities)
    
        cigartuples = build_cigar(read, pairs, start_index, end_index)
        new_MD = get_md_tag(new_pairs)

        if self.bisulfite:
                        
            direction = read.get_tag("YD")
            new_ZC, new_ZR, new_NM = get_biscuit_tags(new_pairs, direction)
            
            new_tags = {
                "ZC" : new_ZC,
                "ZR" : new_ZR,
                "NM" : new_NM,
                "MD" : new_MD
            }
            
        else:
            new_NM = get_bwa_tags(new_pairs)

            new_tags = {
                "NM" : new_NM
            }
    
        tags = read.get_tags()
    
        for i in range(len(tags)):
            tag = tags[i]
            if tag[0] in new_tags:
                tags[i] = (tag[0], new_tags[tag[0]])
    
        real_trim_5 = read_5 - original_span[0]
        real_trim_3 = original_span[1] - read_3 - 1 
    
        if full_bam:
            tags.append(("ZU", real_trim_5))
            tags.append(("ZD", real_trim_3))
        
        new_read = pysam.AlignedSegment(header=pysam.AlignmentHeader().from_dict(header))
        new_read.query_name = read.query_name
        new_read.query_sequence = read_seq
        new_read.flag = read.flag
        new_read.reference_id = read.reference_id
        new_read.reference_start = pos
        new_read.mapping_quality = read.mapping_quality
        new_read.cigar = cigartuples
        new_read.query_qualities = pysam.qualitystring_to_array(read_qual)
        new_read.tags = tags
    
        return new_read, new_pairs
    
    def trim_mate(self):

        ordered_reads = self.ordered_reads
        seg_keys = self.seg_keys
        original_sequence = self.original_sequence
        
        cut_site_keys = {}
        overlap_keys = {}
        cut_site_labels = {}
        cut_site_classes = {}
    
        adjusted_seg_keys = seg_keys.copy()
        
        for i in range(len(adjusted_seg_keys)-1):
            self.seg5 = adjusted_seg_keys[i]
            self.pairs5 = ordered_reads[i]["pairs"]
            self.read5 = ordered_reads[i]["read"]
            
            self.seg3 = adjusted_seg_keys[i+1]
            self.pairs3 = ordered_reads[i+1]["pairs"]
            self.read3 = ordered_reads[i+1]["read"]

            self.overlap = self.seg5[1] - self.seg3[0]
            original_overlap = self.overlap
            
            self.total_pairs += 1

            
            bp_class, bp_enzyme, r5_rs, r3_rs = \
                split_pair_to_restriction_site(self.read5, self.read3, self.restriction_sites, self.max_cut_site_split_algn_dist)

            self.bp_enzyme = bp_enzyme
            self.r5_rs = r5_rs
            self.r3_rs = r3_rs
            
            aseg5, aseg3 = self.adjust_split()
    
            adjusted_seg_keys[i] = aseg5
            adjusted_seg_keys[i+1] = aseg3
            cut_site_keys[(i, i+1)] = bp_class
            cut_site_classes[(i, i+1)] = bp_enzyme
            overlap_keys[(i, i+1)] = original_overlap
    
            if i not in cut_site_labels:
                cut_site_labels[i] = []
            if i+1 not in cut_site_labels:
                cut_site_labels[i+1] = []
    
            if bp_enzyme != "artefact":
                cut_site_labels[i].append("D")
                cut_site_labels[i+1].append("U")
    
        all_keys = list(ordered_reads.keys())
        
        for key in all_keys:
            ordered_reads[key]["adjusted_span"] = adjusted_seg_keys[key]

            if ordered_reads[key]["span"] == ordered_reads[key]["adjusted_span"]:
                trimmed_read = ordered_reads[key]["read"]
                new_pairs = ordered_reads[key]["pairs"]
                read_is_trimmed = False
            else:
                trimmed_read, new_pairs = self.create_trimmed_mate(
                    ordered_reads[key]["read"], 
                    ordered_reads[key]["span"], 
                    ordered_reads[key]["adjusted_span"], 
                    ordered_reads[key]["pairs"]
                )
                read_is_trimmed = True
    
            if trimmed_read == None: # Read is eliminated/erroneous by trimming
                del ordered_reads[key] # Delete read from dictionary 
                self.illegal_post_trim += 1
            else:
                ordered_reads[key]["trimmed_read"] = trimmed_read
                ordered_reads[key]["new_pairs"] = new_pairs
                ordered_reads[key]["is_trimmed"] = read_is_trimmed

        self.cut_site_keys = cut_site_keys
        self.cut_site_labels = cut_site_labels
        self.cut_site_classes = cut_site_classes
        self.overlap_keys = overlap_keys
        
        return
        
    def process_mate(self):

        seg_keys = self.seg_keys
        has_split = self.has_split
        read_parts = self.read_parts
        original_sequence = self.original_sequence
        
        ordered_reads = OrderedDict()

        if len(seg_keys) == 1:
            read = read_parts[seg_keys[0]]
            span = seg_keys[0]
            ordered_reads[0] = {"read" : read,
                                "span" : span,
                                "pairs" : None,
                                "adjusted_span" : span,
                                "trimmed_read" : read,
                                "new_pairs" : None,
                                "is_trimmed" : False
                               }

            self.ordered_reads = ordered_reads
            self.cut_site_keys = {}
            self.cut_site_labels = {}
            self.cut_site_classes = {}
            self.overlap_keys = {}

            return
            
        for i in range(len(seg_keys)):
            read = read_parts[seg_keys[i]]
            span = seg_keys[i]
            ordered_reads[i] = {"read" : read,
                                "span" : span,
                                "pairs" : alignment_info(read, span, original_sequence),
                               }

        self.ordered_reads = ordered_reads

        self.trim_mate()

        return
    
    def __init__(self,
                 seg_keys, 
                 original_sequence = "",
                 read_parts = {},
                 restriction_sites = {},
                 header = None,
                 full_bam = False,
                 trim_method = "winner",
                 bisulfite=False,
                 max_cut_site_split_algn_dist=20,
                 max_cut_site_whole_algn_dist=500
                ):

        self.illegal_post_trim = 0
        self.cut_site_pairs = 0
        self.first_trim = 0
        self.second_trim = 0
        self.total_pairs = 0
        
        if len(seg_keys) == 0:
            self.ordered_reads = OrderedDict()
            self.has_split = False
            self.cut_site_keys = {}
            self.cut_site_labels = {}
            self.overlap_keys = {}
            self.cut_site_classes = {}
            self.original_sequence = None
            return
            
        self.has_split = read_parts[seg_keys[0]].has_tag("SA")
            
        self.seg_keys = seg_keys
        self.original_sequence = original_sequence
        self.read_parts = read_parts
        self.restriction_sites = restriction_sites
        self.header = header
        self.full_bam = full_bam
        
        if trim_method == "winner":
            self.adjust_split = self.adjust_split_win
        elif trim_method == "complete":
            self.adjust_split = self.adjust_split_complete
        elif trim_method == "none":
            self.adjust_split = self.adjust_split_none

        self.bisulfite = bisulfite
        self.max_cut_site_split_algn_dist = max_cut_site_split_algn_dist
        self.max_cut_site_whole_algn_dist = max_cut_site_whole_algn_dist

        self.process_mate()
        