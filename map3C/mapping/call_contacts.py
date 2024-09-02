import pysam
import subprocess
import numpy as np
from collections import OrderedDict

from .utils import *
from .pairtools import *
from .pairs_generator import *
from .read_trimmer import *
from .phase import *

def divide_reads_default(read_group):
    
    r1 = []
    r1_primary = None
    
    r2 = []
    r2_primary = None

    for read in read_group:
        if read.is_read1:
            if not read.is_secondary:
                r1_primary = read
            r1.append(read)
        elif read.is_read2:
            if not read.is_secondary:
                r2_primary = read
            r2.append(read)
        else:
            if not read.is_secondary:
                r1_primary = read
            r1.append(read)

    return r1, r1_primary, r2, r2_primary


def divide_reads_manual_annotation(read_group):
    
    r1 = []
    r1_primary = None
    
    r2 = []
    r2_primary = None

    for read in read_group:
        if read.query_name.split("_")[1] == "1":
            if not read.is_secondary:
                r1_primary = read
            r1.append(read)
        elif read.query_name.split("_")[1] == "2":
            if not read.is_secondary:
                r2_primary = read
            r2.append(read)

    return r1, r1_primary, r2, r2_primary

def get_loc(read, original_sequence):
    cigar = read.cigartuples
    if not cigar:
        return (0, 0, read.mapping_quality)

    if read.is_forward:
        clip_5 = cigar[0][1] if cigar[0][0] in [4, 5] else 0
        clip_3 = cigar[-1][1] if cigar[-1][0] in [4, 5] else 0
    else:
        clip_5 = cigar[-1][1] if cigar[-1][0] in [4, 5] else 0
        clip_3 = cigar[0][1] if cigar[0][0] in [4, 5] else 0

    start = clip_5
    end = len(original_sequence) - clip_3

    return (start, end, read.mapping_quality)

def illegal_overlap(seg_keys):
    illegal_overlap = False
    # If split, determine if sequential overlaps condition is violated
    if len(seg_keys) > 2:
        for i in range(len(seg_keys) - 1):
            key = seg_keys[i]
            next_key = seg_keys[i+1]
            if i == 0:
                previous_key = (-1, -1, 0) # Dummy key
            else:
                previous_key = seg_keys[i-1]
        
            # Current key is contained within previous key
            if key[0] >= previous_key[0] and key[1] <= previous_key[1]:
                illegal_overlap = True
                break
            # Current key contains previous key
            if key[0] <= previous_key[0] and key[1] >= previous_key[1]:
                illegal_overlap = True
                break
            # Current key is contained within next key
            elif key[0] >= next_key[0] and key[1] <= next_key[1]:
                illegal_overlap = True
                break
            # Current key contains next key
            elif key[0] <= next_key[0] and key[1] >= next_key[1]:
                illegal_overlap = True
                break
        
            # Check for overlaps of current key with all other keys that are 
            # not the previous or next key
            for j in range(len(seg_keys)):
                test_seg_key = seg_keys[j]
                if test_seg_key == previous_key:
                    continue
                elif test_seg_key == next_key:
                    continue
                elif test_seg_key == key:
                    continue
                elif test_seg_key[0] < key[0] < test_seg_key[1]:
                    illegal_overlap = True
                    break
                elif test_seg_key[0] < key[1] < test_seg_key[1]:
                    illegal_overlap = True
                    break

    elif len(seg_keys) == 2:
        key = seg_keys[0]
        next_key = seg_keys[1]

        # Key is contained within next key
        if key[0] >= next_key[0] and key[1] <= next_key[1]:
            illegal_overlap = True
        # Key contains next key
        elif key[0] <= next_key[0] and key[1] >= next_key[1]:
            illegal_overlap = True

    return illegal_overlap

def remove_illegal_overlap(seg_keys):
    seg_keys_copy = seg_keys.copy()
    rng.shuffle(seg_keys_copy)
    sorted_seg_keys = sorted(seg_keys_copy, key = lambda x : (x[2], x[1] - x[0]), reverse=True) 

    # This enables skipping already deleted segments
    to_analyze = {}
    for i in sorted_seg_keys:
        to_analyze[i] = True

    for i in range(len(sorted_seg_keys) - 1):
        seg0 = sorted_seg_keys[i]
        if not to_analyze[seg0]:
            continue
            
        for j in range(i + 1, len(sorted_seg_keys)):
            seg1 = sorted_seg_keys[j]

            if seg0[0] >= seg1[0] and seg0[1] <= seg1[1]:
                to_analyze[seg1] = False
            elif seg0[0] <= seg1[0] and seg0[1] >= seg1[1]:
                to_analyze[seg1] = False

    filtered_seg_keys = []
    removed_keys = 0
    for i in seg_keys:
        if to_analyze[i]:
            filtered_seg_keys.append(i)
        else:
            removed_keys += 1

    return filtered_seg_keys, removed_keys
    
class ContactGenerator:
    
    def add_cs_tags(self, **kwargs):
        
        read = kwargs['read_data']["trimmed_read"]
        key = kwargs['key']
        cs_labels = kwargs['cs_labels']
        
        if key in cs_labels:
            cs_label = cs_labels[key]
            if len(cs_label) == 0:
                read.set_tag("ZL", "N")
            elif len(cs_label) == 1:
                read.set_tag("ZL", cs_label[0])
            elif cs_label == ["D", "U"] or cs_label == ["U", "D"]:
                read.set_tag("ZL", "B")

    def add_phase_tags(self, **kwargs):

        trim = kwargs["read_data"]["trimmed_read"]
        pre_trim = kwargs["read_data"]["read"]
        is_trimmed = kwargs["read_data"]["is_trimmed"]
        
        phase_t, status_t, eval_snps_t, a1_t, a2_t = self.read_phaser.phase_read(trim)

        trim.set_tag("ZP", f"{status_t},{eval_snps_t},{a1_t},{a2_t},{phase_t}", value_type="Z")
        
        if is_trimmed:
            phase_pt, status_pt, eval_snps_pt, a1_pt, a2_pt = self.read_phaser.phase_read(pre_trim)
            trim.set_tag("ZQ", f"{status_pt},{eval_snps_pt},{a1_pt},{a2_pt},{phase_pt}", value_type="Z")
    
    def tag_and_write_ordered_reads(self, ordered_reads, cs_labels, primary, has_split, primary_unique, mate):

        # dupsifter marks all secondary alignments as duplicate if primary alignment is duplicate
        if not primary.is_duplicate:
            if has_split:
                self.stats_dict[f"{mate}_split_aligned_mates_dedup"] += 1
            else:
                self.stats_dict[f"{mate}_whole_aligned_mates_dedup"] += 1

            # Read is split, but primary alignment will not be included in processed mates
            if has_split and not primary_unique:
                self.bam_out.write(primary)

        if has_split:
            self.stats_dict[f"{mate}_split_aligned_mates_dup"] += 1
        else:
            self.stats_dict[f"{mate}_whole_aligned_mates_dup"] += 1

        for i in ordered_reads:
            read = ordered_reads[i]["trimmed_read"]
            if not read.is_duplicate:
                
                self.stats_dict[f"{mate}_total_alignments_dedup"] += 1
                
                for func in self.tag_funcs:
                    func(read_data = ordered_reads[i], key=i, cs_labels=cs_labels)
                
                self.bam_out.write(read)
                
            self.stats_dict[f"{mate}_total_alignments_dup"] += 1
                
        
    def process_mate(self, all_alignments, primary_alignment, mate, read_group_name):

        if len(all_alignments) == 0:
            read_trimmer = ReadTrimmer(seg_keys = [])
            return read_trimmer
    
        read_parts = {}
            
        # Get original whole read sequence from primary alignment
        original_sequence = primary_alignment.get_forward_sequence()
        
        for read in all_alignments:
            if read.is_secondary:
                if "S" in read.cigarstring:
                    read.flag = read.flag - 256
            read_parts[get_loc(read, original_sequence)] = read
    
        # Sort segments in order from 5' to 3' by starting position
        # Only choose segments with high MAPQ
        seg_keys = sorted([i for i in list(read_parts.keys()) if i[2] >= self.min_mapq], 
                          key = lambda x : x[0])

        if len(seg_keys) == 0:
            read_trimmer = ReadTrimmer(seg_keys = [])
            return read_trimmer

        # Throw out these reads
        if illegal_overlap(seg_keys):
            seg_keys, removed_keys = remove_illegal_overlap(seg_keys)
            self.illegal_overlap_alignments += removed_keys
            self.illegal_overlap_reads += 1

        #print(mate)
        #print(seg_keys)
        
        read_trimmer = ReadTrimmer(seg_keys, 
                                   original_sequence,
                                   read_parts,
                                   self.restriction_sites,
                                   self.header,
                                   self.full_bam,
                                   self.trim_method,
                                   self.bisulfite,
                                   self.max_cut_site_split_algn_dist,
                                   self.max_cut_site_whole_algn_dist
                                  )
        
        self.illegal_post_trim_count += read_trimmer.illegal_post_trim
        
        return read_trimmer
    
    def process_read_group(self, read_group, read_group_name):

        bam_out = self.bam_out
        
        if read_group == None:
            return

        r1, r1_primary, r2, r2_primary = self.divide_reads(read_group)
        if r1_primary != None:
            r1_primary_unique = r1_primary.mapping_quality >= self.min_mapq
        if r2_primary != None:
            r2_primary_unique = r2_primary.mapping_quality >= self.min_mapq
        
        # Order reads from 5' to 3'
        R1_trimmer = self.process_mate(r1, r1_primary, "1", read_group_name)
        R1 = R1_trimmer.ordered_reads
        R1_overlap_keys = R1_trimmer.overlap_keys
        R1_cs_labels = R1_trimmer.cut_site_labels
        R1_has_split = R1_trimmer.has_split
        
        R2_trimmer = self.process_mate(r2, r2_primary, "2", read_group_name)
        R2 = R2_trimmer.ordered_reads
        R2_overlap_keys = R2_trimmer.overlap_keys
        R2_cs_labels = R2_trimmer.cut_site_labels
        R2_has_split = R2_trimmer.has_split

        self.total_chimera_pairs += R1_trimmer.total_pairs
        self.total_chimera_pairs += R2_trimmer.total_pairs
        
        self.cut_site_chimera_pairs += R1_trimmer.cut_site_pairs
        self.cut_site_chimera_pairs += R2_trimmer.cut_site_pairs

        self.first_trim += R1_trimmer.first_trim
        self.first_trim += R2_trimmer.first_trim
        
        self.second_trim += R1_trimmer.second_trim
        self.second_trim += R2_trimmer.second_trim
        
        if len(R1) + len(R2) >= 2:
            self.at_least_two_alignments += 1

        if len(R1) > 0:
            self.tag_and_write_ordered_reads(R1, R1_cs_labels, r1_primary, R1_has_split, r1_primary_unique, "R1")
        if len(R2) > 0:
            self.tag_and_write_ordered_reads(R2, R2_cs_labels, r2_primary, R2_has_split, r2_primary_unique, "R2")

        contacts, rule = contact_iter(R1, R2, min_mapq=self.min_mapq, 
                                      max_molecule_size=self.max_molecule_size, 
                                      max_inter_align_gap=self.max_inter_align_gap)

        for (hic_algn1, hic_algn2, pair_index) in contacts:

            self.pairs_gen.write_pairs(hic_algn1, hic_algn2, read_group_name, pair_index, R1_trimmer, R2_trimmer, rule)
   
    def process_bam(self):

        self.illegal_overlap_alignments = 0
        self.illegal_overlap_reads = 0
        self.illegal_post_trim_count = 0
        self.at_least_two_alignments = 0
        
        self.total_chimera_pairs = 0
        self.cut_site_chimera_pairs = 0
        self.first_trim = 0
        self.second_trim = 0
        
        self.stats_dict = {
            "R1_total_alignments_dup" : 0,
            "R1_whole_aligned_mates_dup" : 0,
            "R1_split_aligned_mates_dup" : 0,
            "R1_total_aligned_mates_dup" : 0,

            "R1_total_alignments_dedup" : 0,
            "R1_whole_aligned_mates_dedup" : 0,
            "R1_split_aligned_mates_dedup" : 0,
            "R1_total_aligned_mates_dedup" : 0,

            "R2_total_alignments_dup" : 0,
            "R2_whole_aligned_mates_dup" : 0,
            "R2_split_aligned_mates_dup" : 0,
            "R2_total_aligned_mates_dup" : 0,

            "R2_total_alignments_dedup" : 0,
            "R2_whole_aligned_mates_dedup" : 0,
            "R2_split_aligned_mates_dedup" : 0,
            "R2_total_aligned_mates_dedup" : 0,
        }
        
        iter_count = 0
        read_group_name = None
        read_group = None
        
        count = 0
        with pysam.AlignmentFile(self.bam, index_filename=None) as bam_in, \
            pysam.AlignmentFile(self.trimmed_bam, 'wb', template=bam_in) as bam_out, \
            PairsGenerator(self.contacts, 
                           self.chrom_sizes,
                           self.chrom_orders,
                           self.restriction_sites,
                           self.artefacts, 
                           self.blacklist,
                           self.min_blacklist_overlap_length,
                           self.min_blacklist_overlap_ratio,
                           self.read_phaser,
                           self.remove_all,
                           self.full_pairs, 
                           self.flip_reads,
                           self.genome,
                           self.chrom_regex, 
                           self.min_inward_dist_contacts,
                           self.min_outward_dist_contacts,
                           self.min_same_strand_dist_contacts,
                           self.min_inward_dist_artefacts,
                           self.min_outward_dist_artefacts,
                           self.min_same_strand_dist_artefacts,
                           self.include_artefacts,
                           self.max_cut_site_whole_algn_dist) as self.pairs_gen:
                
                self.header = bam_in.header.to_dict()
                self.bam_in = bam_in
                self.bam_out = bam_out

                for read in bam_in:
                    read_id = read.query_name.split("_")[0]
                    if iter_count == 0:
                        read_group_name = read_id
                        read_group = []
                        iter_count = 1
                    count += 1
                    
                    if read_group_name == read_id:
                        read_group.append(read)
                    else:
                        self.process_read_group(read_group, read_group_name)
                        
                        read_group_name = read_id
                        read_group = [read]

                self.process_read_group(read_group, read_group_name)

                self.phase_stats = self.pairs_gen.phase_stats

    def generate_stats(self):

        self.stats_dict.update({
            "R1_total_aligned_mates_dup" : self.stats_dict["R1_whole_aligned_mates_dup"] + self.stats_dict["R1_split_aligned_mates_dup"],
            "R1_total_aligned_mates_dedup" : self.stats_dict["R1_whole_aligned_mates_dedup"] + self.stats_dict["R1_split_aligned_mates_dedup"],

            "R2_total_aligned_mates_dup" : self.stats_dict["R2_whole_aligned_mates_dup"] + self.stats_dict["R2_split_aligned_mates_dup"],
            "R2_total_aligned_mates_dedup" : self.stats_dict["R2_whole_aligned_mates_dup"] + self.stats_dict["R2_split_aligned_mates_dup"],

            "discarded_alignments_illegal_overlap" : self.illegal_overlap_alignments,
            "reads_with_illegal_overlap" : self.illegal_overlap_reads,
            "discarded_alignments_post_trim" : self.illegal_post_trim_count,
            "read_pairs_with_multiple_valid_alignments" : self.at_least_two_alignments,

            "adjacent_split_alignment_total" : self.total_chimera_pairs,
            "adjacent_split_alignment_with_cut_site" : self.cut_site_chimera_pairs,
            "adjacent_split_alignment_first_trim" : self.first_trim,
            "adjacent_split_alignment_second_trim" : self.second_trim,

        })

        if self.variants:

            self.stats_dict.update(self.phase_stats)
            
        stats_df = pd.DataFrame.from_dict(self.stats_dict, orient="index").T
        stats_df.to_csv(self.stats_path, index=False, sep="\t")

    def __init__(self, 
                 bam, 
                 out_prefix, 
                 chrom_sizes,
                 reference_name,
                 restriction_sites,
                 restriction_enzymes,
                 min_mapq=30, 
                 max_molecule_size=750, 
                 max_inter_align_gap=20,
                 trim_reporting="minimal",
                 trim_method="winner",
                 pairs_reporting="minimal",
                 variants=None,
                 phase_alignments=False,
                 min_base_quality=20,
                 chrom_regex=None,
                 blacklist=None,
                 min_blacklist_overlap_length=1,
                 min_blacklist_overlap_ratio=0.5,
                 no_flip=False,
                 remove_all=False,
                 min_inward_dist_contacts=1000,
                 min_outward_dist_contacts=1000,
                 min_same_strand_dist_contacts=0,
                 min_inward_dist_artefacts=1000,
                 min_outward_dist_artefacts=1000,
                 min_same_strand_dist_artefacts=0,
                 include_artefacts=False,
                 read_type="bisulfite",
                 manual_mate_annotation=False,
                 max_cut_site_split_algn_dist = 10,
                 max_cut_site_whole_algn_dist = 500
                ):

        self.min_mapq = min_mapq
        self.max_molecule_size = max_molecule_size
        self.max_inter_align_gap = max_inter_align_gap

        self.min_inward_dist_contacts = min_inward_dist_contacts
        self.min_outward_dist_contacts = min_outward_dist_contacts
        self.min_same_strand_dist_contacts = min_same_strand_dist_contacts

        self.min_inward_dist_artefacts = min_inward_dist_artefacts
        self.min_outward_dist_artefacts = min_outward_dist_artefacts
        self.min_same_strand_dist_artefacts = min_same_strand_dist_artefacts

        self.include_artefacts = include_artefacts

        self.trim_method = trim_method

        self.tag_funcs = []
        
        self.full_bam = trim_reporting == "full"
        if self.full_bam:
            self.tag_funcs.append(self.add_cs_tags)
            
        self.full_pairs = pairs_reporting == "full"
        
        self.out_prefix = out_prefix
        self.bam = bam
        self.bisulfite = read_type == "bisulfite"
        self.divide_reads = divide_reads_manual_annotation if manual_mate_annotation else divide_reads_default
        self.max_cut_site_split_algn_dist = max_cut_site_split_algn_dist
        self.max_cut_site_whole_algn_dist = max_cut_site_whole_algn_dist

        self.restriction_sites = process_restriction_sites(restriction_sites, restriction_enzymes)
        
        self.chrom_sizes_file = chrom_sizes
        self.chrom_sizes = process_chrom_sizes(chrom_sizes)
        self.genome = reference_name
        self.chrom_orders = process_chrom_orders(chrom_sizes)

        self.variants = variants 
        self.phase_alignments = phase_alignments
        if self.phase_alignments:
            self.tag_funcs.append(self.add_phase_tags)
        
        if variants:
            if self.bisulfite:
                self.read_phaser = ReadPhaserBisulfite(variants, min_base_quality)
            else:
                self.read_phaser = ReadPhaser(variants, min_base_quality)
        else:
            self.read_phaser = None

        self.flip_reads = not no_flip
        
        self.chrom_regex=chrom_regex
        if blacklist:
            self.blacklist=process_bed(blacklist) 
        else:
            self.blacklist=None
            
        self.min_blacklist_overlap_length = min_blacklist_overlap_length
        self.min_blacklist_overlap_ratio = min_blacklist_overlap_ratio

        self.remove_all=remove_all
        
        self.contacts = f'{out_prefix}_contacts.pairs.gz'
        self.artefacts = f'{out_prefix}_artefacts.pairs.gz'
        self.stats_path = f"{out_prefix}_alignment_stats.txt" 
        self.trimmed_bam = f'{out_prefix}_trimmed.bam'

        self.process_bam()
        
        self.generate_stats()
