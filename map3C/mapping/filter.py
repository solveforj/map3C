from Bio.bgzf import BgzfWriter, BgzfReader
from .utils import process_bed

class PairsLine:

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
        self.is_enzymeless = "enzymeless" in self.contact_class
        self.is_chimera = "gap" not in self.contact_class
        
    def __str__(self):

        return "\t".join(self.line) + "\n"

class ContactFilter:    

    def __init__(self, 
                 input_pairs, 
                 out_prefix, 
                 enzymeless_split_read_pairs=False,
                 enzyme_pairs=False,
                 enzymeless_pairs=False,
                 min_inward_dist_enzyme=0,
                 min_outward_dist_enzyme=0,
                 min_same_strand_dist_enzyme=0,
                 min_inward_dist_enzymeless=0,
                 min_outward_dist_enzymeless=0,
                 min_same_strand_dist_enzymeless=0,
                 remove_trans_enzymeless=False):

        self.input_pairs = input_pairs
        self.output_pairs = f"{out_prefix}.flt.pairs.gz"
        self.split_reads = f"{out_prefix}.flt.enzymeless_split_reads.pairs.gz" if enzymeless_split_read_pairs else None
        self.enzyme = f"{out_prefix}.flt.enzyme.pairs.gz" if enzyme_pairs else None
        self.enzymeless = f"{out_prefix}.flt.enzymeless.pairs.gz" if enzymeless_pairs else None
        self.min_inward_dist_enzyme = min_inward_dist_enzyme
        self.min_outward_dist_enzyme = min_outward_dist_enzyme
        self.min_same_strand_dist_enzyme = min_same_strand_dist_enzyme
        self.min_inward_dist_enzymeless = min_inward_dist_enzymeless
        self.min_outward_dist_enzymeless = min_outward_dist_enzymeless
        self.min_same_strand_dist_enzymeless = min_same_strand_dist_enzymeless
        self.remove_trans_enzymeless = remove_trans_enzymeless
        
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

        if self.enzyme:
            self.enzyme_writer = BgzfWriter(self.enzyme, 'wb')
            for line in header:
                self.enzyme_writer.write(line)

        if self.enzymeless:
            self.enzymeless_writer = BgzfWriter(self.enzymeless, 'wb')
            for line in header:
                self.enzymeless_writer.write(line)

        funcs = []
        funcs.append(self.contact_pair_is_intra_short)
        funcs.append(self.enzymeless_pair_is_intra_short)        
        if self.remove_trans_enzymeless:
            funcs.append(self.is_trans_enzymeless)
        
        funcs.append(self.write_pair)

        if self.split_reads:
            funcs.append(self.write_sr)

        if self.enzyme:
            funcs.append(self.write_enzyme)

        if self.enzymeless:
            funcs.append(self.write_enzymeless)
            
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

        if self.enzyme_writer:
            self.enzyme_writer.close()

        if self.enzymeless_writer:
            self.enzymeless_writer.close()
            
    def contact_pair_is_intra_short(self, pair):

        if not pair.is_enzymeless:
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
                    if (max_pos - min_pos) < self.min_same_strand_dist_enzyme:
                        pair.passed_filters = False
                elif min_strand == "+":
                    if (max_pos - min_pos) < self.min_inward_dist_enzyme:
                        pair.passed_filters = False
                elif min_strand == "-":
                    if (max_pos - min_pos) < self.min_outward_dist_enzyme:
                        pair.passed_filters = False

    def enzymeless_pair_is_intra_short(self, pair):

        if pair.is_enzymeless:
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
                    if (max_pos - min_pos) < self.min_same_strand_dist_enzymeless:
                        pair.passed_filters = False
                elif min_strand == "+":
                    if (max_pos - min_pos) < self.min_inward_dist_enzymeless:
                        pair.passed_filters = False
                elif min_strand == "-":
                    if (max_pos - min_pos) < self.min_outward_dist_enzymeless:
                        pair.passed_filters = False

    def is_trans_enzymeless(self, pair):

        if pair.is_enzymeless:
            if pair.chrom1 != pair.chrom2:
                    pair.passed_filters = False
    
    def write_pair(self, pair):

        if pair.passed_filters:
            self.pairs_writer.write(str(pair))

    def write_sr(self, pair):

        if pair.is_enzymeless and pair.is_chimera and pair.passed_filters:
            self.sr_writer.write(str(pair))

    def write_enzyme(self, pair):

        if not pair.is_enzymeless and pair.passed_filters:
            self.enzyme_writer.write(str(pair))

    def write_enzymeless(self, pair):

        if pair.is_enzymeless and pair.passed_filters:
            self.enzymeless_writer.write(str(pair))
                 
