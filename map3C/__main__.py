from map3C import __version__

import argparse
import inspect
import subprocess
import sys
import logging
import os


log = logging.getLogger()

DESCRIPTION = """
Pipeline for mapping 3C/Hi-C data

"""

EPILOG = ''

class NiceFormatter(logging.Formatter):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Do not prefix "INFO:" to info-level log messages (but do it for all other
    levels).
    Based on http://stackoverflow.com/a/9218261/715090 .
    """

    def format(self, record):
        if record.levelno != logging.INFO:
            record.msg = '{}: {}'.format(record.levelname, record.msg)
        return super().format(record)

def setup_logging(stdout=False, quiet=False, debug=False):
    """
    From Cutadapt https://github.com/marcelm/cutadapt
    Attach handler to the global logger object
    """
    # Due to backwards compatibility, logging output is sent to standard output
    # instead of standard error if the -o option is used.
    stream_handler = logging.StreamHandler(sys.stdout if stdout else sys.stderr)
    stream_handler.setFormatter(NiceFormatter())
    # debug overrides quiet
    if debug:
        level = logging.DEBUG
    elif quiet:
        level = logging.ERROR
    else:
        level = logging.INFO
    stream_handler.setLevel(level)
    log.setLevel(level)
    log.addHandler(stream_handler)

def prepare_demultiplex_register_subparser(subparser):
    parser = subparser.add_parser('prepare-demultiplex',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Setup demultiplexing")
    
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--config', type=str, default=None, required=True,
                        help='Path to config YML file')

    parser_opt = parser.add_argument_group("optional arguments")
    
    parser_opt.add_argument('--snakemake-params', type=str, default="", 
                            help="Snakemake-specific parameters")

def prepare_mapping_register_subparser(subparser):
    parser = subparser.add_parser('prepare-mapping',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Setup mapping")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument('--config', type=str, default=None, required=True,
                        help='Path to config YML file')

    parser_opt = parser.add_argument_group("optional arguments")
    
    parser_opt.add_argument('--snakemake-params', type=str, default="", 
                            help="Snakemake-specific parameters")

def contamination_filter_register_subparser(subparser):
    parser = subparser.add_parser('contamination-filter',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Select only reads that have low CH methylation and/or small numbers of CH sites")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')

    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output filtered bam')

    parser_opt = parser.add_argument_group("optional arguments")

    parser_opt.add_argument('--min-mapq', type=int, default=30, 
                            help="MAPQ threshold for considering a read's CH methylation")

    parser_opt.add_argument('--max-mc-ch', type=float, default=0.7,  
                            help="Methylated fraction of CH sites threshold")

    parser_opt.add_argument('--max-ch-sites', type=int, default=3,  
                            help="CH sites threshold")
                                

def call_contacts_register_subparser(subparser):
    parser = subparser.add_parser('call-contacts',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Call contacts from BAM file and trim split alignments to remove within-mate multimapping")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output files')

    parser_req.add_argument('--chrom-sizes', type=str, default=None, required=True,
                            help='Path to chromosome sizes file')

    parser_req.add_argument('--reference-name', type=str, default=None, required=True,
                            help='Name of reference genome (i.e. hg38 or mm10)')

    parser_req.add_argument('--restriction-sites', type=str, action="append", nargs="+", default=[], required=True,
                            help="""Paths to restriction sites files. For multiple files, either list sequentially separated by 
                                    spaces or specify this argument multiple times followed by one file.""")

    parser_req.add_argument('--restriction-enzymes', type=str, action="append", nargs="+", default=[], required=True,
                            help="""Names of restriction enzymes. For multiple enzymes, either list sequentially separated by 
                                    spaces or specify this argument multiple times followed by one enzyme. Should be in same order
                                    as --restriction-sites argument.""")

    parser_opt = parser.add_argument_group("optional arguments")
    
    parser_opt.add_argument('--keep-duplicates', action="store_true",
                        help='Keep reads marked as duplicates in output BAM file.')

    parser_opt.add_argument('--no-output-bam', action="store_true",
                        help='Do not write output bam file')
    
    parser_opt.add_argument('--min-mapq', type=int, default=30,
                        help='Minimum MAPQ to consider alignment (Pairtools parameter)')

    parser_opt.add_argument('--max-molecule-size', type=int, default=750,
                        help="""The maximal size of a Hi-C molecule; used to rescue single ligations 
                                (from molecules with three alignments) and to rescue complex ligations.
                                Used for walks-policy mask, not walks-policy all (Pairtools parameter)""")
    
    parser_opt.add_argument('--max-inter-align-gap', type=int, default=20,
                      help="""Read segments that are not covered by any alignment and longer than the 
                              specified value are treated as “null” alignments. These null alignments 
                              convert otherwise linear alignments into walks, and affect how they get reported 
                              as a Hi-C pair (Pairtools parameter)""")

    parser_opt.add_argument('--trim-reads', action="store_true",
                        help="""Remove multimapped region from both split alignments""")
    
    parser_opt.add_argument('--trim-reporting', type=str, default="minimal", choices=['minimal', 'full'],
                            help="""If set, output BAM files have 4 extra tags for split alignments. ZU and ZD report the 
                                    number of basepairs trimmed off of the 5' and 3' end of the alignment, respectively. 
                                    ZL reports if an alignment was determined to have a cut site at the 5' end (U), 3' 
                                    end (D), both (B), or neither (N)""")
    
    parser_opt.add_argument('--min-inward-dist-enzyme', type=int, default=1000,
                            help="""Minimum distance for intrachromosomal contacts with +/- strandedness 
                                    (downstream read/upstream read). Filters out WGS-like reads, such as those due 
                                    to dangling ends.""")
    
    parser_opt.add_argument('--min-outward-dist-enzyme', type=int, default=1000,
                            help="""Minimum distance for intrachromosomal contacts with -/+ strandedness 
                                    (downstream read/upstream read). Filters out self-ligations.""")

    parser_opt.add_argument('--min-same-strand-dist-enzyme', type=int, default=0,
                            help="""Minimum distance for intrachromosomal contacts with +/+ or -/- strandedness 
                                    (downstream read/upstream read). Typically, no cutoff is needed for Hi-C/3C.""")

    parser_opt.add_argument('--min-inward-dist-enzymeless', type=int, default=1000,
                            help="""Minimum distance for intrachromosomal artefacts with +/- strandedness 
                                    (downstream read/upstream read).""")
    
    parser_opt.add_argument('--min-outward-dist-enzymeless', type=int, default=1000,
                            help="""Minimum distance for intrachromosomal artefacts with -/+ strandedness 
                                    (downstream read/upstream read).""")

    parser_opt.add_argument('--min-same-strand-dist-enzymeless', type=int, default=0,
                            help="""Minimum distance for intrachromosomal contacts with +/+ or -/- strandedness 
                                    (downstream read/upstream read).""")

    parser_opt.add_argument('--read-type', type=str, default="bsdna", choices=['bsdna', 'dna'],
                            help='Indicates that reads were bisulfite converted (bsdna) or not bisulfite converted (dna)')

    parser_opt.add_argument('--manual-mate-annotation', action="store_true",
                            help="""If set, input bam file is understood to have mates manually added as suffixes 
                                    to read names (i.e. @readname_1 or @readname_2). This is done in the case of SE 
                                    alignment to distinguish the mates.""")

    parser_opt.add_argument('--max-cut-site-split-algn-dist', type=int, default = 20,
                            help="""Max allowed distance (bp) from nearest cut site to split alignment to be considered ligation event""")

    parser_opt.add_argument('--max-cut-site-whole-algn-dist', type=int, default = 20,
                            help="""Max allowed distance (bp) from nearest cut site to whole alignment to be considered ligation event""")

    parser_opt.add_argument('--pairs-reporting', type=str, default = "minimal", choices=["full", "minimal"],
                            help="""If full, extra metadata will be reported as additional columns in pairs files""")

    parser_opt.add_argument('--variants', type=str, default = None,
                            help="""Path to phased variants file""")

    parser_opt.add_argument('--min-base-quality', type=int, default = 20,
                            help="""If variants are specified for contact phasing, this specifies the minimum base quality for a read nucleotide
                                    to be considered to empirically determine the phase of a relevant read.""")

    parser_opt.add_argument('--chrom-regex', type=str, default = None,
                            help="""Regex to select specific contacts and artefacts from specific chromosomes""")

    parser_opt.add_argument('--blacklist', type=str, default = None,
                            help="""Path to blacklist BED file. Contacts and artefacts that intersect with blacklisted regions will
                                    not be reported""")
    
    parser_opt.add_argument('--min-blacklist-overlap-length', type=int, default = 1,
                            help="""If blacklist is provided, sets minimum number of bp that overlap with blacklisted region for any read
                                    in a potential contact pair""")

    parser_opt.add_argument('--min-blacklist-overlap-ratio', type=float, default = 0.5,
                            help="""If blacklist is provided, sets minimum fraction of bp that overlap with blacklisted region for any read
                                    in a potential contact pair""")
    
    parser_opt.add_argument('--remove-all', action="store_true",
                            help="""If set, then contacts and artefacts discovered exclusively by the pairtools all algorithm will
                                    not be reported. This eliminates detection of multiple ligation events by this tool and is not
                                    recommended to be set.""")

    parser_opt.add_argument('--pair-combinations', action="store_true",
                            help="""If set, then for walk pairs (called using pairtools all algorithm), all combinations of alignments within a given
                                    read pair will be called as pairs. If an alignment combination is not adjacent on the same read, they will be identified 
                                    as artefacts using the --max-cut-site-whole-algn-dist parameter.
                                    """)
    
    parser_opt.add_argument('--no-flip', action="store_true",
                            help="""If set, then contacts and artefacts will be left in their original orientation. By default, they are flipped to
                                    create an upper triangular contact matrix.""")

    parser_opt.add_argument('--phase-bam', action="store_true",
                            help="""If set, then all reported alignments will be phased.Otherwise, only alignments involved in pairs will be phased.
                                    """)


def mask_overlaps_register_subparser(subparser):
    parser = subparser.add_parser('mask-overlaps',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Masks overlapping bases between mates of the same read pair")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam', type=str, default=None, required=True,
                            help='Path to input bam')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output bam file')

    parser_opt = parser.add_argument_group("optional arguments")

    parser_opt.add_argument('--min-mapq', type=int, default=30,
                            help='Minimum MAPQ to consider alignment')

    parser_opt.add_argument('--manual-mate-annotation', action="store_true",
                            help="""If set, input bam file is understood to have mates manually added as suffixes 
                                    to read names (i.e. @readname_1 or @readname_2). This is done in the case of SE 
                                    alignment to distinguish the mates.""")

def bam_to_allc_register_subparser(subparser):
    parser = subparser.add_parser('bam-to-allc',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="Convert a biscuit-derived BAM file to ALLC format")

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--bam-path', type=str, default=None, required=True,
                            help='Path to input bam file')

    parser_req.add_argument('--reference-fasta', type=str, default=None, required=True,
                            help='Path to reference fasta file')
    
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output files')

    parser_opt = parser.add_argument_group("optional arguments")

    parser_opt.add_argument('--control-chrom', type=str, default=None,
                            help='Control chromosome')

    parser_opt.add_argument('--num-upstr-bases', type=int, default=0,
                            help='Number of upstream bases for context')

    parser_opt.add_argument('--num-downstr-bases', type=int, default=2,
                            help='Number of downstream bases for context')
    
    parser_opt.add_argument('--min-mapq', type=int, default=30,
                            help='Minimum MAPQ score for including aligned reads')

    parser_opt.add_argument('--min-base-quality', type=int, default=20,
                            help='Minimum base quality for including aligned nucleotides')
    
    parser_opt.add_argument('--compress-level', type=int, default=5,
                            help='Compression level')
    
    parser_opt.add_argument('--save-count-df', action="store_true",
                            help='If set, save context count summary file')

def pairtools_stats_register_subparser(subparser):
    parser = subparser.add_parser('pairtools-stats',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="""
                                        Compute stats for pairs files for contacts and non-ligation artefacts 
                                        """
                                 )
    # Required arguments
    parser_req = parser.add_argument_group("required arguments")

    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                        help='Path including name prefix for output stats file')

    parser_req.add_argument('--input-pairs', type=str, default=None, required=True,
                            help='Path to contacts pairs file')

    parser_opt = parser.add_argument_group("optional arguments")

    parser_opt.add_argument('--pairs-dedup-stats', type=str, nargs="?", default=None, required=False,
                            help='Path to pairtools dedup stats for contacts pairs file')
    
    parser_opt.add_argument('--pairs-filterbycov-stats', type=str, nargs="?", default=None, required=False,
                            help='Path to pairtools filterbycov stats for contacts pairs file')
                            
    
def aggregate_qc_stats_register_subparser(subparser):
    parser = subparser.add_parser('aggregate-qc-stats',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="""
                                        Aggregate QC stats generated during trimming, filtering, mapping, and contact/methylation calling.
                                        """
                                 )
    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--job', type=str, default=None, required=True,
                            help='Name of job (i.e. cell or experiment name)')
                            
    parser_req.add_argument('--out-prefix', type=str, default=None, required=True,
                            help='Path including name prefix for output stats file')

    parser_req.add_argument('--mode', type=str, default=None, choices=["bsdna", "dna"], required=True,
                        help='Mode')

def restriction_sites_register_subparser(subparser):
    parser = subparser.add_parser('restriction-sites',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="""
                                        Identify restriction enzyme cut sites for a reference genome.
                                        """
                                 )
    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--cut-seqs', type=str, default=[], nargs='+', required=True,
                            help='Cut site sequence for enzyme')
                            
    parser_req.add_argument('--reference', type=str, default=None, required=True,
                            help='Path to reference genome FASTA')

    parser_req.add_argument('--output', type=str, default=None, required=True,
                        help='Full path to output file')

def filter_pairs_register_subparser(subparser):
    parser = subparser.add_parser('filter-pairs',
                                  formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                  description="""
                                        Filter bgzipped pairs file based on various criteria. All filtered reads are added into an output 
                                        bgzipped pairs file. There are also options to additionally create even more bgzipped pairs files 
                                        that only contain certain groups of pairs. 
                                        """
                                 )

    # Required arguments
    parser_req = parser.add_argument_group("required arguments")
    
    parser_req.add_argument('--input-pairs', type=str, required=True,
                            help='Input bgzipped pairs file')

    parser_req.add_argument('--out-prefix', type=str, required=True,
                            help='Output bgzipped pairs file')
                            
    parser_subset = parser.add_argument_group("subsetting pairs (optional arguments)",
                                            """These flags each subset the specified group of pairs to a separate file, given
                                                that the pairs in the respective group pass all other filtering criteria in this function.""")

    parser_subset.add_argument('--enzymeless-split-read-pairs', action="store_true",
                            help="""Generate bgzipped pairs file that only reports pairs resulting from split reads not 
                                    proximal to restriction site""")

    parser_subset.add_argument('--enzymeless-pairs', action="store_true",
                            help='Generate bgzipped pairs file that only reports pairs not proximal to restriction site')
    
    parser_subset.add_argument('--enzyme-pairs', action="store_true",
                            help='Generate bgzipped pairs file that only reports pairs proximal to restriction site')

    parser_filter = parser.add_argument_group("filtering pairs (optional arguments)",
                                            """These parameters filter out specific pairs""")
    
    parser_filter.add_argument('--remove-trans-enzymeless', action="store_true",
                            help='Remove trans artefacts')

    parser_filter.add_argument('--min-inward-dist-enzyme', type=int, default=0,
                            help="""Minimum distance for intrachromosomal contacts with +/- strandedness 
                                    (downstream read/upstream read). Filters out WGS-like reads, such as those due 
                                    to dangling ends.""")
    
    parser_filter.add_argument('--min-outward-dist-enzyme', type=int, default=0,
                            help="""Minimum distance for intrachromosomal contacts with -/+ strandedness 
                                    (downstream read/upstream read). Filters out self-ligations.""")

    parser_filter.add_argument('--min-same-strand-dist-enzyme', type=int, default=0,
                            help="""Minimum distance for intrachromosomal contacts with +/+ or -/- strandedness 
                                    (downstream read/upstream read). Typically, no cutoff is needed for Hi-C/3C.""")

    parser_filter.add_argument('--min-inward-dist-enzymeless', type=int, default=0,
                            help="""Minimum distance for intrachromosomal artefacts with +/- strandedness 
                                    (downstream read/upstream read).""")
    
    parser_filter.add_argument('--min-outward-dist-enzymeless', type=int, default=0,
                            help="""Minimum distance for intrachromosomal artefacts with -/+ strandedness 
                                    (downstream read/upstream read).""")

    parser_filter.add_argument('--min-same-strand-dist-enzymeless', type=int, default=0,
                            help="""Minimum distance for intrachromosomal contacts with +/+ or -/- strandedness 
                                    (downstream read/upstream read).""")
    

def main():
    parser = argparse.ArgumentParser(description=DESCRIPTION,
                                     epilog=EPILOG,
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     )
    subparsers = parser.add_subparsers(
        title="functions",
        dest="command",
        metavar="",
        required=True
    )

    # add subparsers
    current_module = sys.modules[__name__]
    # get all functions in parser
    for name, register_subparser_func in inspect.getmembers(current_module, inspect.isfunction):
        if 'register_subparser' in name:
            register_subparser_func(subparsers)

    # initiate
    args = None
    if len(sys.argv) > 1:
        # print out version
        if sys.argv[1] in ['-v', '--version']:
            print(__version__)
            exit()
        else:
            args = parser.parse_args()
    else:
        # print out help
        parser.parse_args(["-h"])
        exit()

    # set up logging
    if not logging.root.handlers:
        setup_logging(stdout=True,
                      quiet=False)
    # execute command
    args_vars = vars(args)
    for k, v in args_vars.items():
        log.info(f'{k}\t{v}')

    cur_command = args_vars.pop('command').lower().replace('_', '-')
    # Do real import here:
    if cur_command in ['prepare-demultiplex']:
        from .demultiplex import PrepareDemultiplex as func
    elif cur_command in ['prepare-mapping']:
        from .mapping import PrepareMapping as func
    elif cur_command in ['contamination-filter']:
        from .mapping import ContaminationFilter as func
    elif cur_command in ['call-contacts']:
        from .mapping import ContactGenerator as func
    elif cur_command in ['pairtools-stats']:
        from .mapping import pairtools_stats as func
    elif cur_command in ['mask-overlaps']:
        from .mapping import OverlapMask as func
    elif cur_command in ['bam-to-allc']:
        from .mapping import bam_to_allc as func
    elif cur_command in ['aggregate-qc-stats']:
        from .mapping import aggregate_qc_stats as func
    elif cur_command in ['restriction-sites']:
        from .mapping import ComputeRestrictionSites as func
    elif cur_command in ['filter-pairs']:
        from .mapping import ContactFilter as func
    else:
        log.debug(f'{cur_command} is not an valid sub-command')
        parser.parse_args(["-h"])
        return

    # run the command
    log.info(f"Executing {cur_command}...")
    func(**args_vars)
    log.info(f"{cur_command} finished.")
    return
    

if __name__ == '__main__':
    main()