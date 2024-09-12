from Bio.bgzf import BgzfWriter
from .read_trimmer import gap_pair_to_restriction_site
import re
from bisect import bisect_left, bisect_right

class Pair:
    def __init__(self, algn1, 
                 algn2, readID,
                 ct, overlap, rule,
                 cs_locs, pair_index,
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
            "cut_site_locs" : ""
        }

        
        self.flip = False
        
        if "artefact" in ct:
            self.pair_class = "artefact"
            if flip_pairs:
                self._evaluate_flip(chrom_orders)
        elif "na" != ct:
            self.pair_class = "contact"
            if flip_pairs:
                self._evaluate_flip(chrom_orders)
        else:
            self.pair_class = "na"
        
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

    def add_metadata(self):

        self._line.update({
            "rule" : f"\t{self.rule}",
            "reads" : f"\t{self.reads}",
            "contact_class" : f"\t{self.contact_class}",
            "multimap_overlap" : f"\t{self.multimap_overlap}",
            "cut_site_locs" : f"\t{self.cut_site_locs}",
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
                f"{self._line['phase2']}"
                f"{self._line['phase1']}"  
                f"{self._line['rule']}"
                f"{self._line['reads']}"
                f"{self._line['contact_class']}"
                f"{self._line['multimap_overlap']}"
                f"{self._line['cut_site_locs']}\n"
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
                f"{self._line['phase1']}"
                f"{self._line['phase2']}"   
                f"{self._line['rule']}"
                f"{self._line['reads']}"
                f"{self._line['contact_class']}"
                f"{self._line['multimap_overlap']}"
                f"{self._line['cut_site_locs']}\n"
            )
        return line

class PairsGenerator:

    bam_status = {"C": "alignments_with_conflicting_phase",
                  "P": "alignments_with_phase",
                  "N": "alignments_without_phase"}
    
    def write_header(self, handle, artefacts=False):
        
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
        base_columns = "#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type"

        if artefacts:
            if self.read_phaser and self.include_artefacts:
                base_columns += " phase0 phase1"
        else:
            if self.read_phaser:
                base_columns += " phase0 phase1"

            
        if self.full_pairs:
            base_columns += " rule reads contact_class multimap_overlap cut_site_locs"

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

    
    def classify_pair(self, algn1,algn2, pair_index,
                         R1_trimmer,R2_trimmer,rule):

        R1_cs_keys = R1_trimmer.cut_site_keys
        R2_cs_keys = R2_trimmer.cut_site_keys
    
        R1_cs_classes = R1_trimmer.cut_site_classes
        R2_cs_classes = R2_trimmer.cut_site_classes
    
        R1_overlap_keys = R1_trimmer.overlap_keys
        R2_overlap_keys = R2_trimmer.overlap_keys
        
        overlap = 0
        cs_locs = ["na"]
    
        if not algn1["is_mapped"] or not algn1["is_unique"]:
            ct = "na"
            return ct, overlap, cs_locs
        if not algn2["is_mapped"] or not algn2["is_unique"]:
            ct = "na"
            return ct, overlap, cs_locs

        contact_reads = pair_index[1]

        if rule == "all":
            # Pairtools reports 5' fragment before 3' fragment
            if contact_reads == "R1": 
                idx5 = algn1["idx"]
                idx3 = algn2["idx"]
                cs_key = (idx5, idx3)
                if cs_key not in R1_cs_classes:
                    ct = "artefact_chimera"
                elif R1_cs_classes[cs_key] == "artefact":
                    ct = "artefact_chimera"
                    overlap = R1_overlap_keys[cs_key]
                    cs_locs = R1_cs_keys[cs_key]
                else:
                    ct = R1_cs_classes[cs_key]
                    overlap = R1_overlap_keys[cs_key]
                    cs_locs = R1_cs_keys[cs_key]
            # Pairtools reports 5' fragment before 3' fragment
            elif contact_reads == "R2": 
                idx5 = algn1["idx"]
                idx3 = algn2["idx"]
                cs_key = (idx5, idx3)
                
                if cs_key not in R2_cs_classes:
                    ct = "artefact_chimera"
                elif R2_cs_classes[cs_key] == "artefact":
                    ct = "artefact_chimera"
                    overlap = R2_overlap_keys[cs_key]
                    cs_locs = R2_cs_keys[cs_key]
                else:
                    ct = R2_cs_classes[cs_key]
                    overlap = R2_overlap_keys[cs_key]
                    cs_locs = R2_cs_keys[cs_key]
            elif contact_reads == "R1&2" or contact_reads == "R1-2":
                bp_class, bp_enzyme, r5_rs, r3_rs = gap_pair_to_restriction_site(algn1["read"], 
                                                                                 algn2["read"], 
                                                                                 self.restriction_sites, 
                                                                                 self.max_cut_site_whole_algn_dist)
                if bp_enzyme != "artefact":
                    ct = "gap"
                else:
                    ct = "artefact_gap"
                    
            elif contact_reads =="comb":
                bp_class, bp_enzyme, r5_rs, r3_rs = gap_pair_to_restriction_site(algn1["read"], 
                                                                                 algn2["read"], 
                                                                                 self.restriction_sites, 
                                                                                 self.max_cut_site_whole_algn_dist)
                if bp_enzyme != "artefact":
                    ct = "gap"
                else:
                    ct = "artefact_gap"
                    
        #elif algn1["type"] == "R":
        elif contact_reads == "R2_rescue":
            if (0, 1) not in R2_cs_classes:
                ct = "artefact_chimera"
            elif R2_cs_classes[(0, 1)] == "artefact":
                ct = "artefact_chimera"
                overlap = R2_overlap_keys[(0, 1)]
                cs_locs = R2_cs_keys[(0, 1)]
            else:
                ct = R2_cs_classes[(0, 1)]
                overlap = R2_overlap_keys[(0, 1)]
                cs_locs = R2_cs_keys[(0, 1)]
        #elif algn2["type"] == "R":
        elif contact_reads == "R1_rescue":
            if (0, 1) not in R1_cs_classes:
                ct = "artefact_chimera"
            elif R1_cs_classes[(0, 1)] == "artefact":
                ct = "artefact_chimera"
                overlap = R1_overlap_keys[(0, 1)]
                cs_locs = R1_cs_keys[(0, 1)]
            else:
                ct = R1_cs_classes[(0, 1)]
                overlap = R1_overlap_keys[(0, 1)]
                cs_locs = R1_cs_keys[(0, 1)]
    
        elif algn1["type"] == "U" and algn2["type"] == "U":
            bp_class, bp_enzyme, r5_rs, r3_rs = gap_pair_to_restriction_site(algn1["read"], 
                                                                             algn2["read"], 
                                                                             self.restriction_sites, 
                                                                             self.max_cut_site_whole_algn_dist)
            if bp_enzyme != "artefact":
                ct = "gap"
            else:
                ct = "artefact_gap"
    
        return ct, overlap, cs_locs

    def contact_pair_is_intra_short(self, pair):

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
                if (max_algn["pos"] - min_algn["pos"]) < self.min_same_strand_dist_contacts:
                    pair.passed_filters = False
            elif min_algn["strand"] == "+":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_inward_dist_contacts:
                    pair.passed_filters = False
            elif min_algn["strand"] == "-":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_outward_dist_contacts:
                    pair.passed_filters = False

    def artefact_pair_is_intra_short(self, pair):

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
                if (max_algn["pos"] - min_algn["pos"]) < self.min_same_strand_dist_artefacts:
                    pair.passed_filters = False
            elif min_algn["strand"] == "+":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_inward_dist_artefacts:
                    pair.passed_filters = False
            elif min_algn["strand"] == "-":
                if (max_algn["pos"] - min_algn["pos"]) < self.min_outward_dist_artefacts:
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

    def all_pair(self, pair):

        if pair.is_all():
            pair.passed_filters = False

    def write_contact(self, pair):

        if pair.passed_filters:
            self.contacts.write(str(pair))

    def write_artefact(self, pair):

        if pair.passed_filters:
            self.artefacts.write(str(pair))

    def build_artefact_funcs(self):
        funcs = []
        funcs.append(self.artefact_pair_is_intra_short)
        
        if self.blacklist:
            funcs.append(self.blacklist_pair)
        if self.chrom_regex:
            funcs.append(self.chrom_regex_pair)
        if self.read_phaser and self.include_artefacts:
            funcs.append(self.phase_pair)
        if self.full_pairs:
            funcs.append(self.metadata_pair)
        if self.remove_all:
            funcs.append(self.all_pair)
        if self.include_artefacts:
            funcs.append(self.write_contact)

        funcs.append(self.write_artefact)
        
        self.artefact_funcs = funcs

    def build_contact_funcs(self):
        funcs = []
        funcs.append(self.contact_pair_is_intra_short)
        
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

        funcs.append(self.write_contact)
        
        self.contact_funcs = funcs
            

    def write_pairs(self, algn1, algn2, readID, pair_index, R1_trimmer, R2_trimmer, rule):

        ct, overlap, cs_locs = self.classify_pair(algn1, algn2, pair_index, R1_trimmer, R2_trimmer, rule)

        if ct == "na":
            return
        pair = Pair(algn1, algn2, readID, ct, overlap, rule, cs_locs, pair_index, self.chrom_orders, self.flip_pairs)
        if pair.pair_class == "artefact":
                
            for filter in self.artefact_funcs:
                filter(pair)

            #if pair.passed_filters:
            #    self.artefacts.write(str(pair))

        elif pair.pair_class == "contact":

            for filter in self.contact_funcs:
                filter(pair)

            #if pair.passed_filters:
            #    self.contacts.write(str(pair))
        
    def __init__(self, 
                 contacts_path, 
                 chrom_sizes,
                 chrom_orders,
                 restriction_sites,
                 artefacts_path, 
                 blacklist=None,
                 min_blacklist_overlap_length=1,
                 min_blacklist_overlap_ratio=0.5,
                 read_phaser=None,
                 remove_all=True,
                 full_pairs=False,
                 flip_pairs=False,
                 genome="",
                 chrom_regex=None,
                 min_inward_dist_contacts=1000,
                 min_outward_dist_contacts=1000,
                 min_same_strand_dist_contacts=0,
                 min_inward_dist_artefacts=1000,
                 min_outward_dist_artefacts=1000,
                 min_same_strand_dist_artefacts=0,
                 include_artefacts=False,
                 max_cut_site_whole_algn_dist = 500
                ):
        
        self.artefacts_path = artefacts_path
        self.contacts_path = contacts_path
        self.chrom_sizes = chrom_sizes
        self.chrom_orders = chrom_orders
        self.restriction_sites = restriction_sites
        
        self.blacklist = blacklist
        self.min_blacklist_overlap_length = min_blacklist_overlap_length
        self.min_blacklist_overlap_ratio = min_blacklist_overlap_ratio

        self.read_phaser = read_phaser
        
        self.remove_all = remove_all
        self.full_pairs = full_pairs
        self.flip_pairs = flip_pairs
        self.genome = genome
        self.chrom_regex = re.compile(chrom_regex) if chrom_regex else None
        
        self.min_inward_dist_contacts = min_inward_dist_contacts
        self.min_outward_dist_contacts = min_outward_dist_contacts
        self.min_same_strand_dist_contacts = min_same_strand_dist_contacts

        self.min_inward_dist_artefacts = min_inward_dist_artefacts
        self.min_outward_dist_artefacts = min_outward_dist_artefacts
        self.min_same_strand_dist_artefacts = min_same_strand_dist_artefacts

        self.include_artefacts = include_artefacts
        
        self.max_cut_site_whole_algn_dist = max_cut_site_whole_algn_dist

        self.phase_stats = {"alignments_with_conflicting_phase" : 0,
                            "alignments_with_phase" : 0,
                            "alignments_without_phase": 0}
        
        self.build_contact_funcs()
        self.build_artefact_funcs()
        
    def __enter__(self):
        self.artefacts = BgzfWriter(self.artefacts_path, 'wb')
        self.write_header(self.artefacts, artefacts=True)
            
        self.contacts = BgzfWriter(self.contacts_path, 'wb')
        self.write_header(self.contacts)
        return self

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.artefacts.close()
        self.contacts.close()
        