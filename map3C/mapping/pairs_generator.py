from Bio.bgzf import BgzfWriter
from .cut_analysis import gap_pair_to_restriction_site
import re
from bisect import bisect_left, bisect_right

class Pair:
    def __init__(self, algn1, 
                 algn2, readID,
                 ct, overlap, rule,
                 cs_locs, dists, pair_index,
                 chrom_orders, flip_pairs=False
                ):
        
        self.algn1 = algn1
        self.algn2 = algn2
        self.readID = readID

        self.passed_filters = True

        self.chrom1 = algn1["chrom"]
        self.chrom2 = algn2["chrom"]
        self.pos1 = algn1["pos"]
        self.pos2 = algn2["pos"]
        self.strand1 = algn1["strand"]
        self.strand2 = algn2["strand"]
        self.type1 = "R" if algn1['type'] == "X" else algn1["type"]
        self.type2 = "R" if algn2['type'] == "X" else algn2["type"]
        self.phase1 = ""
        self.phase2 = ""
        self.rule = rule
        self.dists1 = dists[0]
        self.dists2 = dists[1]
        self.reads = pair_index[1]
        self.contact_class = ct
        self.multimap_overlap = str(overlap)
        self.cut_site_locs = ",".join(cs_locs)

        self._line = {
            "readID" : f"{self.readID}",
            #"readID" : f"{"."}\t",
            "chrom1" : f"\t{self.chrom1}",
            "pos1" : f"\t{self.pos1}",
            "chrom2" : f"\t{self.chrom2}",
            "pos2" : f"\t{self.pos2}",
            "strand1" : f"\t{self.strand1}",
            "strand2" : f"\t{self.strand2}",
            "phase1" : "",
            "phase2" : "",
            "rule" : "",
            "reads" : "",
            "contact_class" : "",
            "multimap_overlap" : "",
            "cut_site_locs" : "",
            "cut_site_dist1" : "",
            "cut_site_dist2" : "",
        }

        
        self.flip = False
        
        if "enzymeless" in ct:
            self.pair_class = "enzymeless"
            if flip_pairs:
                self._evaluate_flip(chrom_orders)
        elif "na" != ct:
            self.pair_class = "enzyme"
            if flip_pairs:
                self._evaluate_flip(chrom_orders)
        else:
            self.pair_class = "na"

        
        if self.chrom1 == self.chrom2:
            self.chrom_type = "intra"
        else:
            self.chrom_type = "inter"
            
    def is_all(self):
        if self._rule == "all":
            return True
        elif self._rule == "mask":
            return False

    def _evaluate_flip(self, chrom_orders):
        
        chrom1_order = chrom_orders[self.chrom1]
        chrom2_order = chrom_orders[self.chrom2]

        if chrom2_order < chrom1_order:
            self.flip = True
    
        elif chrom2_order == chrom1_order:
            if self.pos2 < self.pos1:
                self.flip = True
        
    def add_phase(self, phase1, phase2):
        self.phase1 = phase1
        self.phase2 = phase2

        self._line.update({
            "phase1" : f"\t{self.phase1}",
            "phase2" : f"\t{self.phase2}",
        })


    def add_contact_class(self):
        self._line.update({
            "contact_class" : f"\t{self.contact_class}",
        })        

    def add_metadata(self):

        self._line.update({
            "rule" : f"\t{self.rule}",
            "reads" : f"\t{self.reads}",
            "multimap_overlap" : f"\t{self.multimap_overlap}",
            "cut_site_locs" : f"\t{self.cut_site_locs}",
            "cut_site_dist1" : f"\t{self.dists1}",
            "cut_site_dist2" : f"\t{self.dists2}",
        })
        
    def __str__(self):

        if self.flip:
            line = (
                f"{self._line["readID"]}"
                f"{self._line['chrom2']}"
                f"{self._line['pos2']}"
                f"{self._line['chrom1']}"
                f"{self._line['pos1']}"
                f"{self._line['strand2']}"
                f"{self._line['strand1']}" 
                f"\t{self.type2}{self.type1}"
                f"{self._line['contact_class']}"
                f"{self._line['phase2']}"
                f"{self._line['phase1']}"
                f"{self._line['rule']}"
                f"{self._line['reads']}"
                f"{self._line['multimap_overlap']}"
                f"{self._line['cut_site_locs']}"
                f"{self._line['cut_site_dist2']}"
                f"{self._line['cut_site_dist1']}\n"
            )
        else:
            line = (
                f"{self._line["readID"]}"
                f"{self._line['chrom1']}"
                f"{self._line['pos1']}"
                f"{self._line['chrom2']}"
                f"{self._line['pos2']}"
                f"{self._line['strand1']}"
                f"{self._line['strand2']}"
                f"\t{self.type1}{self.type2}"
                f"{self._line['contact_class']}"
                f"{self._line['phase1']}"
                f"{self._line['phase2']}"   
                f"{self._line['rule']}"
                f"{self._line['reads']}"
                f"{self._line['multimap_overlap']}"
                f"{self._line['cut_site_locs']}"
                f"{self._line['cut_site_dist1']}"
                f"{self._line['cut_site_dist2']}\n"
            )
        return line

class PairsGenerator:

    bam_status = {"C": "alignments_with_conflicting_phase",
                  "P": "alignments_with_phase",
                  "N": "alignments_without_phase"}
    
    def write_header(self, handle):
        
        handle.write("## pairs format v1.0\n")

        if self.flip_pairs:
            handle.write("#shape: upper triangle\n")
        handle.write(f"#genome_assembly: {self.genome}\n")
        
        for chrom in self.chrom_sizes:
            if self.chrom_regex:
                if self.chrom_regex.match(chrom):
                    handle.write(f'#chromsize: {chrom} {str(self.chrom_sizes[chrom])}\n')
            else:
                handle.write(f'#chromsize: {chrom} {str(self.chrom_sizes[chrom])}\n')
        base_columns = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type contact_class"

        if self.read_phaser:
            base_columns += " phase0 phase1"

        if self.full_pairs:
            base_columns += " rule reads multimap_overlap cut_site_locs cut_site_dist1 cut_site_dist2"

        base_columns += "\n"
        handle.write(base_columns)

    def algn_in_blacklist(self, algn):

        chrom = algn["chrom"]
        
        if chrom not in self.blacklist:
            return False
    
        chrom_bl = self.blacklist[chrom]
        
        insert = bisect_left(chrom_bl, algn["reference_start"], key=lambda x : x["start"])
        if insert == len(chrom_bl):
            search_start = len(chrom_bl) - 1
        elif chrom_bl[insert]["start"] == algn["reference_start"]:
            search_start = insert
        else:
            search_start = insert-1
    
        if search_start < 0:
            search_start = 0
        else:
            search_start = chrom_bl[search_start]["idx"]
    
        overlap = 0
        for i in range(search_start, len(chrom_bl)):
            if chrom_bl[i]["start"] >= algn["reference_end"] - 1:
                break
            elif chrom_bl[i]["end"] > algn["reference_start"]:
                min_end = min(chrom_bl[i]["end"], algn["reference_end"])
                max_start = max(chrom_bl[i]["start"], algn["reference_start"])
                if min_end >= max_start:
                    overlap += min_end - max_start
                else:
                    raise Exception("Bug with overlap")
        if (overlap >= (algn["reference_end"] - algn["reference_start"]) * self.min_blacklist_overlap_ratio and 
            overlap >= self.min_blacklist_overlap_length):
            return True

        return False

    
    def classify_pair(self,
                      algn1,
                      algn2, 
                      pair_index,
                      R1_readcut,
                      R2_readcut,
                      rule):

        R1_pairwise_cut_site_options = R1_readcut.pairwise_cut_site_options
        R2_pairwise_cut_site_options = R2_readcut.pairwise_cut_site_options
    
        R1_pairwise_cut_site_assign = R1_readcut.pairwise_cut_site_assign
        R2_pairwise_cut_site_assign = R2_readcut.pairwise_cut_site_assign
    
        R1_pairwise_overlaps = R1_readcut.pairwise_overlaps
        R2_pairwise_overlaps = R2_readcut.pairwise_overlaps

        R1_pairwise_cut_site_dists = R1_readcut.pairwise_cut_site_dists
        R2_pairwise_cut_site_dists = R2_readcut.pairwise_cut_site_dists
        
        overlap = 0
        cs_options = ["na"]
        dists = ("na", "na")
    
        if not algn1["is_mapped"] or not algn1["is_unique"]:
            ct = "na"
            return ct, overlap, cs_options, dists
        if not algn2["is_mapped"] or not algn2["is_unique"]:
            ct = "na"
            return ct, overlap, cs_options, dists

        contact_reads = pair_index[1]

        if rule == "all":
            # Pairtools reports 5' fragment before 3' fragment
            if contact_reads == "R1": 
                idx5 = algn1["idx"]
                idx3 = algn2["idx"]
                cs_key = (idx5, idx3)
                if cs_key not in R1_pairwise_cut_site_assign:
                    ct = "enzymeless_chimera"
                elif R1_pairwise_cut_site_assign[cs_key] == "enzymeless":
                    ct = "enzymeless_chimera"
                    overlap = R1_pairwise_overlaps[cs_key]
                    cs_options = R1_pairwise_cut_site_options[cs_key]
                    dists = R1_pairwise_cut_site_dists[cs_key]
                else:
                    ct = R1_pairwise_cut_site_assign[cs_key]
                    overlap = R1_pairwise_overlaps[cs_key]
                    cs_options = R1_pairwise_cut_site_options[cs_key]
                    dists = R1_pairwise_cut_site_dists[cs_key]
            # Pairtools reports 5' fragment before 3' fragment
            elif contact_reads == "R2": 
                idx5 = algn1["idx"]
                idx3 = algn2["idx"]
                cs_key = (idx5, idx3)
                
                if cs_key not in R2_pairwise_cut_site_assign:
                    ct = "enzymeless_chimera"
                elif R2_pairwise_cut_site_assign[cs_key] == "enzymeless":
                    ct = "enzymeless_chimera"
                    overlap = R2_pairwise_overlaps[cs_key]
                    cs_options = R2_pairwise_cut_site_options[cs_key]
                    dists = R2_pairwise_cut_site_dists[cs_key]
                else:
                    ct = R2_pairwise_cut_site_assign[cs_key]
                    overlap = R2_pairwise_overlaps[cs_key]
                    cs_options = R2_pairwise_cut_site_options[cs_key]
                    dists = R2_pairwise_cut_site_dists[cs_key]
            elif contact_reads in ["R1&2", "R1-2", "comb"]:
                if algn1["mate"] == "R1":
                    alignment1 = R1_readcut.ordered_reads[algn1["idx"]]
                elif algn1["mate"] == "R2":
                    alignment1 = R2_readcut.ordered_reads[algn1["idx"]]

                if algn2["mate"] == "R1":
                    alignment2 = R1_readcut.ordered_reads[algn2["idx"]]
                elif algn2["mate"] == "R2":
                    alignment2 = R2_readcut.ordered_reads[algn2["idx"]]
                
                _, bp_enzyme, _, _ = gap_pair_to_restriction_site(alignment1, alignment2, self.max_cut_site_whole_algn_dist)

                r5_rs = alignment1.gap_cut_site3
                r3_rs = alignment2.gap_cut_site3
                r5_dist = ",".join([f"{i}_{r5_rs[i]["dist"]}" for i in r5_rs])
                r3_dist = ",".join([f"{i}_{r3_rs[i]["dist"]}" for i in r3_rs])
                dists = (r5_dist, r3_dist)
                
                if bp_enzyme != "enzymeless":
                    ct = "gap"
                else:
                    ct = "enzymeless_gap"

        elif contact_reads == "R2_rescue":
            if (0, 1) not in R2_pairwise_cut_site_assign:
                ct = "enzymeless_chimera"
            elif R2_pairwise_cut_site_assign[(0, 1)] == "enzymeless":
                ct = "enzymeless_chimera"
                overlap = R2_pairwise_overlaps[(0, 1)]
                cs_options = R2_pairwise_cut_site_options[(0, 1)]
                dists = R2_pairwise_cut_site_dists[(0, 1)]
            else:
                ct = R2_pairwise_cut_site_assign[(0, 1)]
                overlap = R2_pairwise_overlaps[(0, 1)]
                cs_options = R2_pairwise_cut_site_options[(0, 1)]
                dists = R2_pairwise_cut_site_dists[(0, 1)]
        elif contact_reads == "R1_rescue":
            if (0, 1) not in R1_pairwise_cut_site_assign:
                ct = "enzymeless_chimera"
            elif R1_pairwise_cut_site_assign[(0, 1)] == "enzymeless":
                ct = "enzymeless_chimera"
                overlap = R1_pairwise_overlaps[(0, 1)]
                cs_options = R1_pairwise_cut_site_options[(0, 1)]
                dists = R1_pairwise_cut_site_dists[(0, 1)]
            else:
                ct = R1_pairwise_cut_site_assign[(0, 1)]
                overlap = R1_pairwise_overlaps[(0, 1)]
                cs_options = R1_pairwise_cut_site_options[(0, 1)]
                dists = R1_pairwise_cut_site_dists[(0, 1)]
    
        elif algn1["type"] == "U" and algn2["type"] == "U":
            if algn1["mate"] == "R1":
                alignment1 = R1_readcut.ordered_reads[algn1["idx"]]
            elif algn1["mate"] == "R2":
                alignment1 = R2_readcut.ordered_reads[algn1["idx"]]

            if algn2["mate"] == "R1":
                alignment2 = R1_readcut.ordered_reads[algn2["idx"]]
            elif algn2["mate"] == "R2":
                alignment2 = R2_readcut.ordered_reads[algn2["idx"]]

            _, bp_enzyme, _, _ = gap_pair_to_restriction_site(alignment1, alignment2, self.max_cut_site_whole_algn_dist)

            r5_rs = alignment1.gap_cut_site3
            r3_rs = alignment2.gap_cut_site3
            r5_dist = ",".join([f"{i}_{r5_rs[i]["dist"]}" for i in r5_rs])
            r3_dist = ",".join([f"{i}_{r3_rs[i]["dist"]}" for i in r3_rs])
            dists = (r5_dist, r3_dist)

            if bp_enzyme != "enzymeless":
                ct = "gap"
            else:
                ct = "enzymeless_gap"

        return ct, overlap, cs_options, dists

    def enzyme_pair_is_intra_short(self, pair):

        algn1 = pair.algn1
        algn2 = pair.algn2
        
        if algn1["chrom"] == algn2["chrom"]:
            if algn1["pos"] < algn2["pos"]:
                min_algn = algn1
                max_algn = algn2
            else:
                min_algn = algn2
                max_algn = algn1
    
            if max_algn["strand"] == min_algn["strand"]:
                if (max_algn["pos"] - min_algn["pos"]) < self.min_same_strand_dist_enzyme:
                    pair.passed_filters = False
            elif min_algn["strand"] == "+":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_inward_dist_enzyme:
                    pair.passed_filters = False
            elif min_algn["strand"] == "-":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_outward_dist_enzyme:
                    pair.passed_filters = False

    def enzymeless_pair_is_intra_short(self, pair):

        algn1 = pair.algn1
        algn2 = pair.algn2
        
        if algn1["chrom"] == algn2["chrom"]:
            if algn1["pos"] < algn2["pos"]:
                min_algn = algn1
                max_algn = algn2
            else:
                min_algn = algn2
                max_algn = algn1
    
            if max_algn["strand"] == min_algn["strand"]:
                if (max_algn["pos"] - min_algn["pos"]) < self.min_same_strand_dist_enzymeless:
                    pair.passed_filters = False
            elif min_algn["strand"] == "+":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_inward_dist_enzymeless:
                    pair.passed_filters = False
            elif min_algn["strand"] == "-":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_outward_dist_enzymeless:
                    pair.passed_filters = False

    
    def chrom_regex_pair(self, pair):

        if not self.chrom_regex.match(pair.algn1["chrom"]):
            pair.passed_filters = False

        if not self.chrom_regex.match(pair.algn2["chrom"]):
            pair.passed_filters = False
            
    def blacklist_pair(self, pair):

        algn1_blk = self.algn_in_blacklist(pair.algn1)
        algn2_blk = self.algn_in_blacklist(pair.algn2)

        if algn1_blk or algn2_blk:
            
            pair.passed_filters = True

    def phase_pair(self, pair):

        read1 = pair.algn1["read"]
        read2 = pair.algn2["read"]

        if read1.has_tag("ZP"):
            tag1 = read1.get_tag("ZP").split(",")
            phase1 = tag1[-1]
            status1 = tag1[0]
        else:
            phase1, status1, _, _, _ = self.read_phaser.phase_read(read1)
         
        if read2.has_tag("ZP"):
            tag2 = read2.get_tag("ZP").split(",")
            phase2 = tag2[-1]
            status2 = tag2[0]
        else:
            phase2, status2, _, _, _ = self.read_phaser.phase_read(read2)

        self.phase_stats[PairsGenerator.bam_status[status1]] += 1
        self.phase_stats[PairsGenerator.bam_status[status2]] += 1
        
        pair.add_phase(phase1, phase2)

    def metadata_pair(self, pair):

        pair.add_metadata()

    def contact_class_pair(self, pair):

        pair.add_contact_class()

    def all_pair(self, pair):

        if pair.is_all():
            pair.passed_filters = False

    def write_contact(self, pair):

        if pair.passed_filters:
            self.contacts.write(str(pair))

    def build_enzymeless_funcs(self):
        funcs = []
        funcs.append(self.enzymeless_pair_is_intra_short)
        
        if self.blacklist:
            funcs.append(self.blacklist_pair)
        if self.chrom_regex:
            funcs.append(self.chrom_regex_pair)
        if self.read_phaser:
            funcs.append(self.phase_pair)
        if self.full_pairs:
            funcs.append(self.metadata_pair)
        if self.remove_all:
            funcs.append(self.all_pair)

        funcs.append(self.contact_class_pair)
        funcs.append(self.write_contact)
        
        self.enzymeless_funcs = funcs

    def build_enzyme_funcs(self):
        funcs = []
        funcs.append(self.enzyme_pair_is_intra_short)
        
        if self.blacklist:
            funcs.append(self.blacklist_pair)
        if self.chrom_regex:
            funcs.append(self.chrom_regex_pair)
        if self.read_phaser:
            funcs.append(self.phase_pair)
        if self.full_pairs:
            funcs.append(self.metadata_pair)
        if self.remove_all:
            funcs.append(self.all_pair)

        funcs.append(self.contact_class_pair)
        funcs.append(self.write_contact)
        
        self.enzyme_funcs = funcs
            

    def write_pairs(self, c, readID, R1_readcut, R2_readcut, rule):

        algn1 = c[0]
        algn2 = c[1]
        pair_index = c[2]

        ct, overlap, cs_locs, dists = self.classify_pair(algn1, algn2, pair_index, R1_readcut, R2_readcut, rule)

        if ct == "na":
            return
        pair = Pair(algn1, algn2, readID, ct, overlap, rule, cs_locs, dists, pair_index, self.chrom_orders, self.flip_pairs)

        
        if pair.pair_class == "enzymeless":
            
            self.pair_stats[f"potential_pairs_{pair.chrom_type}_enzymeless"] += 1
            
            for filter in self.enzymeless_funcs:
                filter(pair)

        elif pair.pair_class == "enzyme":

            self.pair_stats[f"potential_pairs_{pair.chrom_type}_enzyme"] += 1

            for filter in self.enzyme_funcs:
                filter(pair)

        # TODO: Add leg orient (cis or trans) stats
    
    def __init__(self, 
                 contacts_path, 
                 chrom_sizes,
                 chrom_orders,
                 blacklist=None,
                 min_blacklist_overlap_length=1,
                 min_blacklist_overlap_ratio=0.5,
                 read_phaser=None,
                 remove_all=True,
                 full_pairs=False,
                 flip_pairs=False,
                 genome="",
                 chrom_regex=None,
                 min_inward_dist_enzyme=1000,
                 min_outward_dist_enzyme=1000,
                 min_same_strand_dist_enzyme=0,
                 min_inward_dist_enzymeless=1000,
                 min_outward_dist_enzymeless=1000,
                 min_same_strand_dist_enzymeless=0,
                 max_cut_site_whole_algn_dist = 500
                ):
        
        self.contacts_path = contacts_path
        self.chrom_sizes = chrom_sizes
        self.chrom_orders = chrom_orders
        
        self.blacklist = blacklist
        self.min_blacklist_overlap_length = min_blacklist_overlap_length
        self.min_blacklist_overlap_ratio = min_blacklist_overlap_ratio

        self.read_phaser = read_phaser
        
        self.remove_all = remove_all
        self.full_pairs = full_pairs
        self.flip_pairs = flip_pairs
        self.genome = genome
        self.chrom_regex = re.compile(chrom_regex) if chrom_regex else None
        
        self.min_inward_dist_enzyme = min_inward_dist_enzyme
        self.min_outward_dist_enzyme = min_outward_dist_enzyme
        self.min_same_strand_dist_enzyme = min_same_strand_dist_enzyme

        self.min_inward_dist_enzymeless = min_inward_dist_enzymeless
        self.min_outward_dist_enzymeless = min_outward_dist_enzymeless
        self.min_same_strand_dist_enzymeless = min_same_strand_dist_enzymeless
        
        self.max_cut_site_whole_algn_dist = max_cut_site_whole_algn_dist

        self.pair_stats = {"potential_pairs_intra_enzymeless" : 0,
                           "potential_pairs_intra_enzyme": 0,
                           "potential_pairs_inter_enzymeless" : 0,
                           "potential_pairs_inter_enzyme": 0}
        
        self.phase_stats = {"alignments_with_conflicting_phase" : 0,
                            "alignments_with_phase" : 0,
                            "alignments_without_phase": 0}
        
        self.build_enzyme_funcs()
        self.build_enzymeless_funcs()
        
    def __enter__(self):
        self.contacts = BgzfWriter(self.contacts_path, 'wb')
        self.write_header(self.contacts)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.contacts.close()
        