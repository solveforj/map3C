from Bio.bgzf import BgzfWriter, BgzfReader
from .utils import process_bed

class PairsLine:

    # readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type rule reads contact_class multimap_overlap cut_site_locs
    def __init__(self, line, column_dict, minimal=True):
        line = line.strip().split()
        self.line = line
        self.readID = line[column_dict["readID"]]
        self.chrom1 = line[column_dict["chrom1"]]
        self.pos1 = int(line[column_dict["pos1"]])
        self.chrom2 = line[column_dict["chrom2"]]
        self.pos2 = int(line[column_dict["pos2"]])
        self.strand1 = line[column_dict["strand1"]]
        self.strand2 = line[column_dict["strand2"]]
        self.pair_type = line[column_dict["pair_type"]]
        self.contact_class = line[column_dict["contact_class"]]

        if not minimal:
            self.rule = line[column_dict["rule"]]
            self.multimap_overlap = line[column_dict["multimap_overlap"]]
            self.cut_site_locs = line[column_dict["cut_site_locs"]]

        self.passed_filters = True
        self.is_artefact = "artefact" in self.contact_class
        self.is_chimera = "gap" not in self.contact_class
        
    def __str__(self):

        return "\t".join(self.line) + "\n"

class ContactFilter:    

    def __init__(self, 
                 input_pairs, 
                 out_prefix, 
                 split_reads=False,
                 min_inward_dist_contacts=0,
                 min_outward_dist_contacts=0,
                 min_same_strand_dist_contacts=0,
                 min_inward_dist_artefacts=0,
                 min_outward_dist_artefacts=0,
                 min_same_strand_dist_artefacts=0,
                 remove_trans_artefacts=False):

        self.input_pairs = input_pairs
        self.output_pairs = f"{out_prefix}.flt.pairs.gz"
        self.split_reads = f"{out_prefix}.split_reads.pairs.gz" if split_reads else None
        self.min_inward_dist_contacts = min_inward_dist_contacts
        self.min_outward_dist_contacts = min_outward_dist_contacts
        self.min_same_strand_dist_contacts = min_same_strand_dist_contacts
        self.min_inward_dist_artefacts = min_inward_dist_artefacts
        self.min_outward_dist_artefacts = min_outward_dist_artefacts
        self.min_same_strand_dist_artefacts = min_same_strand_dist_artefacts
        self.remove_trans_artefacts = remove_trans_artefacts
        
        self.initialize()

        self.process()

        self.close()

    def initialize(self):

        self.pairs_reader = BgzfReader(self.input_pairs, 'rt')
        header = []

        self.pairs_writer = None
        self.sr_writer = None
        
        for line in self.pairs_reader:
            if line[0] == "#":
                header.append(line)
            if line.startswith("#columns: "):
                columns = line.replace("#columns: ", "").strip().split()
                break
        if "contact_class" not in columns:
            self.close()
            raise Exception("Input file is missing metadata")

        self.column_dict = {}
        for i in range(len(columns)):
            self.column_dict[columns[i]] = i

        if "rule" in self.column_dict:
            self.minimal = False
        else:
            self.minimal = True
        
        self.pairs_writer = BgzfWriter(self.output_pairs, 'wb')
        
        for line in header:
            self.pairs_writer.write(line)

        if self.split_reads:
            self.sr_writer = BgzfWriter(self.split_reads, 'wb')
            for line in header:
                self.sr_writer.write(line)

        funcs = []
        funcs.append(self.contact_pair_is_intra_short)
        funcs.append(self.artefact_pair_is_intra_short)        
        if self.remove_trans_artefacts:
            funcs.append(self.is_trans_artefact)
        
        funcs.append(self.write_pair)

        if self.split_reads:
            funcs.append(self.write_sr)
            
        self.pair_funcs = funcs

    def process(self):

        for line in self.pairs_reader:
            pair = PairsLine(line, self.column_dict, self.minimal)

            for filter in self.pair_funcs:
                filter(pair)

    def close(self):

        if self.pairs_reader:
            self.pairs_reader.close()
            
        if self.pairs_writer:
            self.pairs_writer.close()
        
        if self.sr_writer:
            self.sr_writer.close()
            
    def contact_pair_is_intra_short(self, pair):

        if not pair.is_artefact:
            if pair.chrom1 == pair.chrom2:
                if pair.pos1 < pair.pos2:
                    min_pos = pair.pos1
                    min_strand = pair.strand1
                    max_pos = pair.pos2
                    max_strand = pair.strand2
                else:
                    min_pos = pair.pos2
                    min_strand = pair.strand2
                    max_pos = pair.pos1
                    max_strand = pair.strand1
        
                if max_strand == min_strand:
                    if (max_pos - min_pos) < self.min_same_strand_dist_contacts:
                        pair.passed_filters = False
                elif min_strand == "+":
                    if (max_pos - min_pos) < self.min_inward_dist_contacts:
                        pair.passed_filters = False
                elif min_strand == "-":
                    if (max_pos - min_pos) < self.min_outward_dist_contacts:
                        pair.passed_filters = False

    def artefact_pair_is_intra_short(self, pair):

        if pair.is_artefact:
            if pair.chrom1 == pair.chrom2:
                if pair.pos1 < pair.pos2:
                    min_pos = pair.pos1
                    min_strand = pair.strand1
                    max_pos = pair.pos2
                    max_strand = pair.strand2
                else:
                    min_pos = pair.pos2
                    min_strand = pair.strand2
                    max_pos = pair.pos1
                    max_strand = pair.strand1
        
                if max_strand == min_strand:
                    if (max_pos - min_pos) < self.min_same_strand_dist_artefacts:
                        pair.passed_filters = False
                elif min_strand == "+":
                    if (max_pos - min_pos) < self.min_inward_dist_artefacts:
                        pair.passed_filters = False
                elif min_strand == "-":
                    if (max_pos - min_pos) < self.min_outward_dist_artefacts:
                        pair.passed_filters = False

    def is_trans_artefact(self, pair):

        if pair.is_artefact:
            if pair.chrom1 != pair.chrom2:
                    pair.passed_filters = False
    
    def write_pair(self, pair):

        if pair.passed_filters:
            self.pairs_writer.write(str(pair))

    def write_sr(self, pair):

        if pair.is_artefact and pair.is_chimera:
            self.sr_writer.write(str(pair))
                 
