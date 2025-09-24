import pysam
import subprocess
import numpy as np
from collections import OrderedDict

from .utils import *
from .pairtools import *
from .pairs_generator import *
from .cut_analysis import *
from .phase import *

def divide_reads_default(read_group):
    
    r1 = []
    r1_primary = None
    
    r2 = []
    r2_primary = None

    for read in read_group:
        if not read.cigarstring:
            continue
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
        if not read.cigarstring:
            continue
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
    

class Alignment:

    def assign_restriction_sites(self, nosplit_cs_direction, restriction_sites):
        
        # Want the closest cut site to 5' and 3' end of read
        self.split_cut_site3 = closest_restriction_site(self.chrom, self.pos3, restriction_sites, rule = "closest")
        self.split_cut_site5 = closest_restriction_site(self.chrom, self.pos5, restriction_sites, rule = "closest")
        
        # Only care about cut site at 3' end of read
        self.gap_cut_site3 = closest_restriction_site(self.chrom, self.pos3, restriction_sites, rule=nosplit_cs_direction)
        self.gap_cut_site5 = None

    def __init__(self, read, span, restriction_sites):

        if read.mapping_quality == 0:
            return 
        
        self.has_split = read.has_tag("SA")
        
        self.read = read
        self.trimmed_read = read

        self.span = span
        self.adjusted_span = span

        self.pairs = None
        self.new_pairs = None

        self.is_trimmed = False

        self.chrom = read.reference_name
        
        if read.is_forward:
            self.pos5 = read.reference_start
            self.pos3 = read.reference_end - 1
            nosplit_cs_direction = "downstream"
        else:
            self.pos5 = read.reference_end - 1
            self.pos3 = read.reference_start
            nosplit_cs_direction = "upstream"

        self.assign_restriction_sites(nosplit_cs_direction, restriction_sites)


class AlignmentNoEnzyme(Alignment):

    def assign_restriction_sites(self, nosplit_cs_direction, restriction_sites):
        
        # Want the closest cut site to 5' and 3' end of read
        self.split_cut_site3 = {}
        self.split_cut_site5 = {}
        
        # Only care about cut site at 3' end of read
        self.gap_cut_site3 = {}
        self.gap_cut_site5 = None

class ContactGenerator:
    
    def add_cs_tags(self, **kwargs):
        
        read = kwargs['read_data'].trimmed_read
        key = kwargs['key']
        cs_labels = kwargs['cs_tag_info']
        
        if key in cs_labels:
            cs_label = cs_labels[key]
            if len(cs_label) == 0:
                read.set_tag("ZL", "N")
            elif len(cs_label) == 1:
                read.set_tag("ZL", cs_label[0])
            elif cs_label == ["D", "U"] or cs_label == ["U", "D"]:
                read.set_tag("ZL", "B")

    def add_phase_tags(self, **kwargs):

        trim = kwargs["read_data"].trimmed_read
        pre_trim = kwargs["read_data"].read
        is_trimmed = kwargs["read_data"].is_trimmed
        
        phase_t, status_t, eval_snps_t, a1_t, a2_t = self.read_phaser.phase_read(trim)

        trim.set_tag("ZP", f"{status_t},{eval_snps_t},{a1_t},{a2_t},{phase_t}", value_type="Z")
        
        if is_trimmed:
            phase_pt, status_pt, eval_snps_pt, a1_pt, a2_pt = self.read_phaser.phase_read(pre_trim)
            trim.set_tag("ZQ", f"{status_pt},{eval_snps_pt},{a1_pt},{a2_pt},{phase_pt}", value_type="Z")

    def tag_and_write_ordered_reads_no_write(self, readcuts, primary, primary_unique, mate):

        has_split = readcuts.has_split

        if has_split:
            self.stats_dict[f"{mate}_split_aligned_mates"] += 1
        else:
            self.stats_dict[f"{mate}_whole_aligned_mates"] += 1

        ordered_reads = readcuts.ordered_reads
        for i in ordered_reads:
            read = ordered_reads[i].trimmed_read
                
            self.stats_dict[f"{mate}_total_alignments"] += 1
            
    def tag_and_write_ordered_reads_dup(self, readcuts, primary, primary_unique, mate):

        has_split = readcuts.has_split

        # Read is split, but primary alignment will not be included in processed mates
        if has_split:
            self.stats_dict[f"{mate}_split_aligned_mates"] += 1
            if not primary_unique:
                self.bam_out.write(primary)
        else:
            self.stats_dict[f"{mate}_whole_aligned_mates"] += 1
                
        ordered_reads = readcuts.ordered_reads
        for i in ordered_reads:
            read = ordered_reads[i].trimmed_read

            self.stats_dict[f"{mate}_total_alignments"] += 1
                        
            for func in self.tag_funcs:
                func(read_data = ordered_reads[i], key=i, cs_tag_info=readcuts.cut_site_tag_info)
            
            self.bam_out.write(read)

    def tag_and_write_ordered_reads_dedup(self, readcuts, primary, primary_unique, mate):

        has_split = readcuts.has_split

        if has_split:
            self.stats_dict[f"{mate}_split_aligned_mates"] += 1
        else:
            self.stats_dict[f"{mate}_whole_aligned_mates"] += 1
            
        # Read is split, but primary alignment will not be included in processed mates
        if not primary.is_duplicate:
            if has_split:
                self.stats_dict[f"{mate}_split_aligned_mates_dedup"] += 1
                if not primary_unique:
                    self.bam_out.write(primary)                
            else:
                self.stats_dict[f"{mate}_whole_aligned_mates_dedup"] += 1

        ordered_reads = readcuts.ordered_reads
        for i in ordered_reads:
            read = ordered_reads[i].trimmed_read

            self.stats_dict[f"{mate}_total_alignments"] += 1
            
            if read.is_duplicate:
                continue
                
            self.stats_dict[f"{mate}_total_alignments_dedup"] += 1
            
            for func in self.tag_funcs:
                func(read_data = ordered_reads[i], key=i, cs_tag_info=readcuts.cut_site_tag_info)
            
            self.bam_out.write(read)
               
        
    def process_mate(self, all_alignments, primary_alignment):

        if len(all_alignments) == 0:
            readcuts = self.cutanalysis(seg_keys = [])
            return readcuts
    
        read_parts = {}

        # Get original whole read sequence from primary alignment
        original_sequence = primary_alignment.get_forward_sequence()

        primary_count = 0
        
        for read in all_alignments:
            if read.is_secondary:
                if "S" in read.cigarstring:
                    read.flag = read.flag - 256
            else:
                primary_count += 1
            loc = get_loc(read, original_sequence)
            read_parts[loc] = self.alignment(read, loc, self.restriction_sites)

        if primary_count > 1:
            readcuts = self.cutanalysis(seg_keys = [])
            self.multiple_primary_reads += 1
            return readcuts
            
        # Sort segments in order from 5' to 3' by starting position
        # Only choose segments with high MAPQ
        seg_keys = sorted([i for i in list(read_parts.keys()) if i[2] >= self.min_mapq], 
                          key = lambda x : x[0])

        if len(seg_keys) == 0:
            readcuts = self.cutanalysis(seg_keys = [])
            return readcuts

        # Throw out these reads
        if illegal_overlap(seg_keys):
            seg_keys, removed_keys = remove_illegal_overlap(seg_keys)
            self.illegal_overlap_alignments += removed_keys
            self.illegal_overlap_reads += 1
        
        readcuts = self.cutanalysis(seg_keys, 
                                    original_sequence,
                                    read_parts,
                                    self.header,
                                    self.full_bam,
                                    self.aligner,
                                    self.max_cut_site_split_algn_dist,
                                    self.max_cut_site_whole_algn_dist
                                  )
        
        self.illegal_post_trim_count += readcuts.illegal_post_trim
        
        return readcuts
    
    def process_read_group(self, read_group, read_group_name):
        
        if read_group == None:
            return

        r1, r1_primary, r2, r2_primary = self.divide_reads(read_group)
        
        if r1_primary != None:
            r1_primary_unique = r1_primary.mapping_quality >= self.min_mapq

        if r2_primary != None:
            r2_primary_unique = r2_primary.mapping_quality >= self.min_mapq
        
        # Order reads from 5' to 3'

        R1_readcuts = self.process_mate(r1, r1_primary)
        R2_readcuts = self.process_mate(r2, r2_primary)

        #print("Done")
        
        # if read_group_name == "x":
        #     print("R1")
        #     for i in R1_readcuts.ordered_reads:
        #         item = R1_readcuts.ordered_reads[i]
        #         print(item.span, item.adjusted_span, "+" if item.read.is_forward else "-")
        #         print(item.read.reference_start, item.read.reference_end)
        #         print(item.trimmed_read.reference_start, item.trimmed_read.reference_end)
        #     print("R2")
        #     for i in R2_readcuts.ordered_reads:
        #         item = R2_readcuts.ordered_reads[i]
        #         print(item.span, item.adjusted_span, "+" if item.read.is_forward else "-")
        #         print(item.read.reference_start, item.read.reference_end)
        #         print(item.trimmed_read.reference_start, item.trimmed_read.reference_end)
            
        
        self.total_chimera_pairs += R1_readcuts.total_pairs
        self.total_chimera_pairs += R2_readcuts.total_pairs
        
        self.cut_site_chimera_pairs += R1_readcuts.cut_site_pairs
        self.cut_site_chimera_pairs += R2_readcuts.cut_site_pairs

        self.first_trim += R1_readcuts.first_trim
        self.first_trim += R2_readcuts.first_trim
        
        self.second_trim += R1_readcuts.second_trim
        self.second_trim += R2_readcuts.second_trim
        
        if len(R1_readcuts.ordered_reads) + len(R2_readcuts.ordered_reads) >= 2:
            self.at_least_two_alignments += 1

        if len(R1_readcuts.ordered_reads) > 0:
            self.tag_and_write_ordered_reads(R1_readcuts, r1_primary, r1_primary_unique, "R1")
        if len(R2_readcuts.ordered_reads) > 0:
            self.tag_and_write_ordered_reads(R2_readcuts, r2_primary, r2_primary_unique, "R2")

        contacts, rule = contact_iter(R1_readcuts.ordered_reads, 
                                      R2_readcuts.ordered_reads, 
                                      min_mapq=self.min_mapq, 
                                      max_molecule_size=self.max_molecule_size, 
                                      max_inter_align_gap=self.max_inter_align_gap,
                                      comb=self.comb
                                     )

        for c in contacts:

            # if read_group_name == "x":
            #     print(c[0])
            #     print(c[1])
            #     print()

            self.pairs_gen.write_pairs(c, read_group_name, R1_readcuts, R2_readcuts, rule)
   
    def process_bam(self):

        self.illegal_overlap_alignments = 0
        self.illegal_overlap_reads = 0
        self.illegal_post_trim_count = 0
        self.at_least_two_alignments = 0
        self.multiple_primary_reads = 0
        
        self.total_chimera_pairs = 0
        self.cut_site_chimera_pairs = 0
        self.first_trim = 0
        self.second_trim = 0
        
        self.stats_dict = {
            "R1_total_alignments" : 0,
            "R1_whole_aligned_mates" : 0,
            "R1_split_aligned_mates" : 0,

            "R2_total_alignments" : 0,
            "R2_whole_aligned_mates" : 0,
            "R2_split_aligned_mates" : 0,

        }

        if not self.keep_duplicates:
            self.stats_dict.update({
    
                "R1_total_alignments_dedup" : 0,
                "R1_whole_aligned_mates_dedup" : 0,
                "R1_split_aligned_mates_dedup" : 0,
                
                "R2_total_alignments_dedup" : 0,
                "R2_whole_aligned_mates_dedup" : 0,
                "R2_split_aligned_mates_dedup" : 0,
                
            })
        
        iter_count = 0
        read_group_name = None
        read_group = None
        
        count = 0
        
        with pysam.AlignmentFile(self.bam, index_filename=None) as bam_in, \
            PairsGenerator(self.contacts, 
                           self.chrom_sizes,
                           self.chrom_orders,
                           self.blacklist,
                           self.min_blacklist_overlap_length,
                           self.min_blacklist_overlap_ratio,
                           self.read_phaser,
                           self.remove_all,
                           self.full_pairs, 
                           self.flip_reads,
                           self.genome,
                           self.chrom_regex, 
                           self.min_inward_dist_enzyme,
                           self.min_outward_dist_enzyme,
                           self.min_same_strand_dist_enzyme,
                           self.min_inward_dist_enzymeless,
                           self.min_outward_dist_enzymeless,
                           self.min_same_strand_dist_enzymeless,
                           self.max_cut_site_whole_algn_dist) as self.pairs_gen:

                if not self.no_output_bam:
                    bam_out = pysam.AlignmentFile(self.trimmed_bam, 'wb', template=bam_in) 
                    self.bam_out = bam_out
                    
                self.bam_in = bam_in
                
                for read in bam_in:
                    read_id = read.query_name.replace("/", " ").replace("_", " ").split()[0]
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
                self.pair_stats = self.pairs_gen.pair_stats

                if not self.no_output_bam:
                    bam_out.close()

    def generate_stats(self):

        self.stats_dict.update({

            "discarded_alignments_illegal_overlap" : self.illegal_overlap_alignments,
            "reads_with_illegal_overlap" : self.illegal_overlap_reads,
            "discarded_reads_multiple_primary" : self.multiple_primary_reads,
            "discarded_alignments_post_trim" : self.illegal_post_trim_count,
            "read_pairs_with_multiple_valid_alignments" : self.at_least_two_alignments,

            "adjacent_split_alignment_total" : self.total_chimera_pairs,
            "adjacent_split_alignment_with_cut_site" : self.cut_site_chimera_pairs,
            "adjacent_split_alignment_first_trim" : self.first_trim,
            "adjacent_split_alignment_second_trim" : self.second_trim,

        })

        self.stats_dict.update(self.pair_stats)
        
        if self.variants:

            self.stats_dict.update(self.phase_stats)

        stats_df = pd.DataFrame.from_dict(self.stats_dict, orient="index").T
        stats_df.to_csv(self.stats_path, index=False, sep="\t")

    def get_aligner(self):
        # Stores header and finds aligner from input BAM file
        
        with pysam.AlignmentFile(self.bam, index_filename=None) as bam_in:
            self.header = bam_in.header.to_dict()
            for i in self.header["PG"]:
                cl = i["CL"]
                if "bsbolt" in cl:
                   self.aligner = "bsbolt"
                   break
                if "biscuit" in cl:
                    self.aligner = "biscuit"
                    break
                else:
                    self.aligner = "bwa mem"
                    
            print(self.aligner)

            if self.variants:
                if self.aligner in ["biscuit", "bsbolt"]:
                    self.read_phaser = ReadPhaserBisulfite(self.variants, self.min_base_quality)
                else:
                    self.read_phaser = ReadPhaser(self.variants, self.min_base_quality)
            else:
                self.read_phaser = None
    
    def __init__(self, 
                 bam, 
                 out_prefix, 
                 chrom_sizes,
                 reference_name,
                 restriction_sites,
                 restriction_enzymes,
                 mate_annotation,
                 keep_duplicates=False,
                 no_output_bam=False,
                 min_mapq=30, 
                 max_molecule_size=750, 
                 max_inter_align_gap=20,
                 trim_reporting="minimal",
                 trim_reads=False,
                 pairs_reporting="minimal",
                 variants=None,
                 phase_bam=False,
                 min_base_quality=20,
                 chrom_regex=None,
                 blacklist=None,
                 min_blacklist_overlap_length=1,
                 min_blacklist_overlap_ratio=0.5,
                 no_flip=False,
                 remove_all=False,
                 min_inward_dist_enzyme=1000,
                 min_outward_dist_enzyme=1000,
                 min_same_strand_dist_enzyme=0,
                 min_inward_dist_enzymeless=1000,
                 min_outward_dist_enzymeless=1000,
                 min_same_strand_dist_enzymeless=0,
                 max_cut_site_split_algn_dist = 10,
                 max_cut_site_whole_algn_dist = 500,
                 pair_combinations=False
                ):

        self.min_mapq = min_mapq
        self.max_molecule_size = max_molecule_size
        self.max_inter_align_gap = max_inter_align_gap

        self.min_inward_dist_enzyme = min_inward_dist_enzyme
        self.min_outward_dist_enzyme = min_outward_dist_enzyme
        self.min_same_strand_dist_enzyme = min_same_strand_dist_enzyme

        self.min_inward_dist_enzymeless = min_inward_dist_enzymeless
        self.min_outward_dist_enzymeless = min_outward_dist_enzymeless
        self.min_same_strand_dist_enzymeless = min_same_strand_dist_enzymeless

        self.tag_funcs = []
        
        if trim_reads:
            self.cutanalysis = CutAnalysis
            if trim_reporting == "full":
                self.full_bam = True
                self.tag_funcs.append(self.add_cs_tags)
            elif trim_reporting == "minimal":
                self.full_bam = False
        else:
            self.full_bam = False
            self.cutanalysis = CutAnalysisNoTrim            
            
        self.full_pairs = pairs_reporting == "full"
        
        self.out_prefix = out_prefix
        self.bam = bam

        # Figure out aligner from SAM header

        if mate_annotation == "qname":
            self.divide_reads = divide_reads_manual_annotation
        else:
            self.divide_reads = divide_reads_default

        self.max_cut_site_split_algn_dist = max_cut_site_split_algn_dist
        self.max_cut_site_whole_algn_dist = max_cut_site_whole_algn_dist

        self.keep_duplicates = keep_duplicates
        
        if self.keep_duplicates:
            self.tag_and_write_ordered_reads = self.tag_and_write_ordered_reads_dup
        else:
            self.tag_and_write_ordered_reads = self.tag_and_write_ordered_reads_dedup

        self.no_output_bam = no_output_bam
        if no_output_bam:
            self.tag_and_write_ordered_reads = self.tag_and_write_ordered_reads_no_write
        
        self.chrom_sizes_file = chrom_sizes
        self.chrom_sizes = process_chrom_sizes(chrom_sizes)
        self.genome = reference_name
        self.chrom_orders = process_chrom_orders(chrom_sizes)

        if len(restriction_sites) > 0 and len(restriction_enzymes) > 0:
            self.restriction_sites = process_restriction_sites(restriction_sites, restriction_enzymes, self.chrom_sizes)
            self.alignment = Alignment
        else:
            self.restriction_sites = {}
            self.alignment = AlignmentNoEnzyme
            

        self.variants = variants 
        self.phase_bam = phase_bam
        if self.phase_bam:
            self.tag_funcs.append(self.add_phase_tags)

        self.flip_reads = not no_flip
        
        self.chrom_regex=chrom_regex
        if blacklist:
            self.blacklist=process_bed(blacklist) 
        else:
            self.blacklist=None

        self.comb = pair_combinations

        self.min_base_quality = min_base_quality
        self.min_blacklist_overlap_length = min_blacklist_overlap_length
        self.min_blacklist_overlap_ratio = min_blacklist_overlap_ratio

        self.remove_all=remove_all
        
        self.contacts = f'{out_prefix}_map3C.pairs.gz'
        self.stats_path = f"{out_prefix}_alignment_stats.txt" 
        self.trimmed_bam = f'{out_prefix}_map3C.bam'

        self.get_aligner()

        self.process_bam()
        
        self.generate_stats()
