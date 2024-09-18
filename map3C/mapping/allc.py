"""
This file is modified from methylpy https://github.com/yupenghe/methylpy.

Author: Yupeng He

Notes added on 07/03/2022 - Difference between bismark and hisat-3n
Read mpileup doc first: http://www.htslib.org/doc/samtools-mpileup.html

For Bismark SE mapping:
Bismark converted the reads orientation based on their conversion type.
C to T conversion reads are always map to the forward strand, regardless of the R1 R2
G to A conversion reads are always map to the reverse strand, regardless of the R1 R2
Therefore, in the original bam-to-allc function, we simply consider the strandness in
the pileup format to distinguish which C needs to be counted.
There are two situations:
1. If the ref base is C, the read need to map to forward strand in order to be counted
[.ATCGN] corresponding to the forward strand
2. If the ref base is G, the read need to map to reverse strand in order to be counted
[,atcgn] corresponding to the reverse strand

For Hisat-3n PE or Biskarp PE mapping:
R1 R2 are mapped as their original orientation, therefore,
both the C to T and G to A conversion reads can have forward and reverse strand alignment.
We can not distinguish the conversion type by strand in input bam.
Here I add a check using the YZ tag of hisat-3n BAM file or XG tag of bismark BAM file.
If the YZ tag is "+" or XG tag is "CT", the read is C to T conversion, I change the flag to forward mapping
no matter R1 or R2 by read.is_forward = True
If the YZ tag is "-" or XG tag is "GA", the read is G to A conversion, I change the flag to reverse mapping
no matter R1 or R2 by read.is_forward = False
In this case, the read orientation is the same as bismark bam file, and the following base
count code no need to change.

"""

import argparse
import collections
import logging
import pathlib
import shlex
import subprocess

import pandas as pd
import pysam
import numpy as np

# logger
log = logging.getLogger(__name__)
log.addHandler(logging.NullHandler())

import codecs
import gzip
import os
import sys
import time
from subprocess import PIPE, Popen, run

from .utils import compute_genome_coverage, compute_mapped_nucleotides

try:
    run(["pigz", "--version"], capture_output=True)
    PIGZ = True
except OSError:
    PIGZ = False

try:
    run(["bgzip", "--version"], capture_output=True)
    BGZIP = True
except OSError:
    BGZIP = False


class Closing:
    """
    Closeable class

    Inherit from this class and implement a close() method to offer context
    manager functionality.
    """

    def close(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, *exc_info):
        self.close()

    def __del__(self):
        try:
            self.close()
        except OSError:
            pass
            
class PipedGzipWriter(Closing):
    """
    Piped Gzip Writer.

    Write gzip-compressed files by running an external gzip or pigz process and
    piping into it. On Python 2, this is faster than using gzip.open(). On
    Python 3, it allows to run the compression in a separate process and can
    therefore also be faster.
    """

    def __init__(self, path, mode="wt", compresslevel=6, threads=1):
        """
        Pipe a file to an external gzip or pigz process.

        mode -- one of 'w', 'wt', 'wb', 'a', 'at', 'ab'
        compresslevel -- gzip compression level
        threads (int) -- number of pigz threads (None means to let pigz decide)
        """
        if mode not in ("w", "wt", "wb", "a", "at", "ab"):
            raise ValueError(f"Mode is '{mode}', but it must be 'w', 'wt', 'wb', 'a', 'at' or 'ab'")

        self.outfile = open(path, mode)
        self.devnull = open(os.devnull, mode)
        self.closed = False
        self.name = path

        kwargs = {"stdin": PIPE, "stdout": self.outfile, "stderr": self.devnull}
        # Setting close_fds to True in the Popen arguments is necessary due to
        # <http://bugs.python.org/issue12786>.
        # However, close_fds is not supported on Windows. See
        # <https://github.com/marcelm/cutadapt/issues/315>.
        if sys.platform != "win32":
            kwargs["close_fds"] = True

        # only apply to gzip and pigz, bgzip use -l
        if "w" in mode and compresslevel != 6:
            extra_args = ["-" + str(compresslevel)]
        else:
            extra_args = []

        if "b" not in mode:
            encoding = None
        else:
            encoding = "utf8"
        kwargs["encoding"] = encoding

        try:
            if BGZIP:
                bgzip_args = ["bgzip"]
                if threads is not None and threads > 0:
                    bgzip_args += ["-@", str(threads), "-l", str(compresslevel)]
                self.process = Popen(bgzip_args, **kwargs)
                self.program = "bgzip"
            elif PIGZ:
                pigz_args = ["pigz"]
                if threads is not None and threads > 0:
                    pigz_args += ["-p", str(threads)]
                self.process = Popen(pigz_args + extra_args, **kwargs)
                self.program = "pigz"
            else:
                # pigz not found, try regular gzip
                self.process = Popen(["gzip"] + extra_args, **kwargs)
                self.program = "gzip"
        except OSError:
            self.outfile.close()
            self.devnull.close()
            raise
        self._file = codecs.getwriter("utf-8")(self.process.stdin)

    def write(self, arg):
        self._file.write(arg)

    def close(self):
        self.closed = True
        self._file.close()
        return_code = self.process.wait()
        self.outfile.close()
        self.devnull.close()
        if return_code != 0:
            raise OSError(f"Output {self.program} process terminated with exit code {return_code}")


class PipedGzipReader(Closing):
    # decompression can't be parallel even in pigz, so there is not thread/cpu parameter
    def __init__(self, path, region=None, mode="r"):
        if mode not in ("r", "rt", "rb"):
            raise ValueError(f"Mode is {mode}, but it must be 'r', 'rt' or 'rb'")
        if "b" not in mode:
            encoding = "utf8"
        else:
            encoding = None

        if region is None:
            self.process = Popen(["gzip", "-cd", path], stdout=PIPE, stderr=PIPE, encoding=encoding)
        else:
            self.process = Popen(
                ["tabix", path] + region.split(" "),
                stdout=PIPE,
                stderr=PIPE,
                encoding=encoding,
            )

        self.name = path
        self._file = self.process.stdout
        self._stderr = self.process.stderr
        self.closed = False
        # Give gzip a little bit of time to report any errors
        # (such as a non-existing file)
        time.sleep(0.01)
        self._raise_if_error()

    def close(self):
        self.closed = True
        return_code = self.process.poll()
        if return_code is None:
            # still running
            self.process.terminate()
        self._raise_if_error()

    def __iter__(self):
        yield from self._file
        self.process.wait()
        self._raise_if_error()

    def readline(self):
        return self._file.readline()

    def _raise_if_error(self):
        """
        Raise an exception if the gzip process has exited with an error.

        Raise IOError if process is not running anymore and the
        exit code is nonzero.
        """
        return_code = self.process.poll()
        if return_code is not None and return_code != 0:
            message = self._stderr.read().strip()
            raise OSError(message)

    def read(self, *args):
        data = self._file.read(*args)
        if len(args) == 0 or args[0] <= 0:
            # wait for process to terminate until we check the exit code
            self.process.wait()
        self._raise_if_error()
        return data
        

def open_gz(file_path, mode="r", compresslevel=3, threads=1, region=None):
    if region is not None:
        if not isinstance(region, str):
            raise TypeError("region parameter need to be string.")

    if "r" in mode:
        try:
            return PipedGzipReader(file_path, region=region, mode=mode)
        except OSError:
            # gzip not installed
            return gzip.open(file_path, mode)
    else:
        try:
            return PipedGzipWriter(file_path, mode, compresslevel, threads=threads)
        except OSError:
            return gzip.open(file_path, mode, compresslevel=compresslevel)

class PipedBamReader(Closing):
    # decompression can't be parallel even in pigz, so there is not thread/cpu parameter
    def __init__(self, path, region=None, mode="r", include_header=True, samtools_parms_str=None):
        if mode not in ("r", "rt", "rb"):
            raise ValueError(f"Mode is {mode}, but it must be 'r', 'rt' or 'rb'")
        if "b" not in mode:
            encoding = "utf8"
        else:
            encoding = None

        command_list = ["samtools", "view"]
        if include_header:
            command_list.append("-h")
        if samtools_parms_str is not None:
            command_list.extend(samtools_parms_str.split(" "))

        if region is None:
            self.process = Popen(command_list + [path], stdout=PIPE, stderr=PIPE, encoding=encoding)
        else:
            self.process = Popen(
                command_list + [path] + region.split(" "),
                stdout=PIPE,
                stderr=PIPE,
                encoding=encoding,
            )

        self.name = path
        self.file = self.process.stdout
        self._stderr = self.process.stderr
        self.closed = False
        # Give gzip a little bit of time to report any errors
        # (such as a non-existing file)
        time.sleep(0.01)
        self._raise_if_error()

    def close(self):
        self.closed = True
        return_code = self.process.poll()
        if return_code is None:
            # still running
            self.process.terminate()
        self._raise_if_error()

    def __iter__(self):
        yield from self.file
        self.process.wait()
        self._raise_if_error()

    def readline(self):
        return self.file.readline()

    def _raise_if_error(self):
        """Raise IOError if process is not running anymore and the exit code is nonzero."""
        return_code = self.process.poll()
        if return_code is not None and return_code != 0:
            message = self._stderr.read().strip()
            raise OSError(message)

    def read(self, *args):
        data = self.file.read(*args)
        if len(args) == 0 or args[0] <= 0:
            # wait for process to terminate until we check the exit code
            self.process.wait()
        self._raise_if_error()
        return data


class PipedBamWriter(Closing):
    def __init__(self, path, mode="wt", threads=1):
        """
        Open a pipe to samtools view, which will write to the given path.

        mode -- one of 'w', 'wt', 'wb', 'a', 'at', 'ab'
        threads (int) -- number of samtools threads
        """
        if mode not in ("w", "wt", "wb", "a", "at", "ab"):
            raise ValueError(f"Mode is '{mode}', but it must be 'w', 'wt', 'wb', 'a', 'at' or 'ab'")

        self.outfile = open(path, mode)
        self.devnull = open(os.devnull, mode)
        self.closed = False
        self.name = path

        kwargs = {"stdin": PIPE, "stdout": self.outfile, "stderr": self.devnull}
        # Setting close_fds to True in the Popen arguments is necessary due to
        # <http://bugs.python.org/issue12786>.
        # However, close_fds is not supported on Windows. See
        # <https://github.com/marcelm/cutadapt/issues/315>.
        if sys.platform != "win32":
            kwargs["close_fds"] = True

        try:
            samtools_args = ["samtools", "view", "-b", "-@", str(threads)]
            self.process = Popen(samtools_args, **kwargs)
            self.program = "samtools"
        except OSError:
            self.outfile.close()
            self.devnull.close()
            raise
        self._file = codecs.getwriter("utf-8")(self.process.stdin)
        return

    def write(self, arg):
        self._file.write(arg)

    def close(self):
        self.closed = True
        self._file.close()
        return_code = self.process.wait()
        self.outfile.close()
        self.devnull.close()
        if return_code != 0:
            raise OSError(f"Output {self.program} process terminated with exit code {return_code}")

def has_bai(filename):
    index_path = filename + ".bai"
    if os.path.exists(index_path):
        return True
    else:
        return False

def open_bam(
    file_path,
    mode="r",
    region=None,
    include_header=True,
    samtools_parms_str=None,
    threads=1,
):
    if "r" in mode:
        if region is not None:
            if not has_bai(file_path):
                raise FileNotFoundError(f".bai index file not found for {file_path}. Index before region query.")
        try:
            return PipedBamReader(
                file_path,
                region=region,
                mode=mode,
                include_header=include_header,
                samtools_parms_str=samtools_parms_str,
            )
        except OSError as e:
            raise e
    else:
        return PipedBamWriter(file_path, mode, threads=threads)

            
def open_allc(file_path, mode="r", compresslevel=3, threads=1, region=None):
    """
    Open a .allc file.

    A replacement for the "open" function that can also open files that have
    been compressed with gzip, bzip2 or xz. If the file_path is '-', standard
    output (mode 'w') or input (mode 'r') is returned.

    The file type is determined based on the file_path: .gz is gzip, .bz2 is bzip2 and .xz is
    xz/lzma.

    When writing a gzip-compressed file, the following methods are tried in order to get the
    best speed 1) using a pigz (parallel gzip) subprocess; 2) using a gzip subprocess;
    3) gzip.open. A single gzip subprocess can be faster than gzip.open because it runs in a
    separate process.

    Uncompressed files are opened with the regular open().

    mode can be: 'rt', 'rb', 'at', 'ab', 'wt', or 'wb'. Also, the 't' can be omitted,
    so instead of 'rt', 'wt' and 'at', the abbreviations 'r', 'w' and 'a' can be used.

    threads is the number of threads for pigz. If None, then the pigz default is used.
    multi-thread only apply to writer, reader (decompression) can't be paralleled
    """
    file_path = pathlib.Path(file_path).resolve()
    file_exist = file_path.exists()
    if "w" not in mode:
        if not file_exist:
            raise FileNotFoundError(f"{file_path} does not exist.")
    else:
        if not file_path.parent.exists():
            raise OSError(f"File directory {file_path.parent} does not exist.")
    file_path = str(file_path)

    if mode in ("r", "w", "a"):
        mode += "t"
    if mode not in ("rt", "rb", "wt", "wb", "at", "ab"):
        raise ValueError(f"mode '{mode}' not supported")
    if compresslevel not in range(1, 10):
        raise ValueError("compresslevel must be between 1 and 9")

    if file_path.endswith("gz"):
        return open_gz(file_path, mode, compresslevel, threads, region=region)
    else:
        return open(file_path, mode)

def _is_read_ct_conversion_biscuit(read):
    return read.get_tag("YD") == "f"

def _convert_bam_strandness(in_bam_path, out_bam_path):
    with pysam.AlignmentFile(in_bam_path) as in_bam, pysam.AlignmentFile(
        out_bam_path, header=in_bam.header, mode="wb"
    ) as out_bam:
        is_ct_func = None
        for read in in_bam:
            if is_ct_func is None:
                if read.has_tag("YD"):
                    is_ct_func = _is_read_ct_conversion_biscuit
                else:
                    raise ValueError(
                        "The bam file reads has no conversion type tag "
                        "(XG by bismark or YZ by hisat-3n). Please note that this function can "
                        "only process bam files generated by bismark or hisat-3n."
                    )
            if is_ct_func(read):
                read.is_forward = True
                if read.is_paired:
                    read.mate_is_forward = True
            else:
                read.is_forward = False
                if read.is_paired:
                    read.mate_is_forward = False
            out_bam.write(read)
    return

def get_stats(output_path, mCH_command, mCG_command):
    command1 = f"zcat {output_path}"
    
    p1 = subprocess.Popen(shlex.split(command1), stdout=subprocess.PIPE)
    p2 = subprocess.run(shlex.split(mCH_command), stdin=p1.stdout, capture_output=True, text=True)
    
    out = p2.stdout.strip().split("\t")
    mch = int(out[0])
    ch = int(out[1])
    
    if ch > 0:
        mch_fraction = mch / ch
    else:
        mch_fraction = np.nan
    
    p1 = subprocess.Popen(shlex.split(command1), stdout=subprocess.PIPE)
    p2 = subprocess.run(shlex.split(mCG_command), stdin=p1.stdout, capture_output=True, text=True)
    
    out = p2.stdout.strip().split("\t")
    mcg = int(out[0])
    cg = int(out[1])
    
    if cg > 0:
        mcg_fraction = mcg / cg
    else:
        mcg_fraction = np.nan

    return mcg_fraction, mch_fraction

def global_stats(output_path, control_chrom):

    if control_chrom:
        mCH_command = f"awk -v OFS='\t' -v mCH=0 -v CH=0 '{{ if (($1 != \"{control_chrom}\") && ($4 ~ /^C[ATC]/)) {{mCH += $5; CH += $6 }}}} END {{print mCH, CH }} '"
        mCG_command = f"awk -v OFS='\t' -v mCG=0 -v CG=0 '{{ if (($1 != \"{control_chrom}\") && ($4 ~ /^C[G]/)) {{mCG += $5; CG += $6 }}}} END {{print mCG, CG }} '"
    else:
        mCH_command = f"awk -v OFS='\t' -v mCH=0 -v CH=0 '{{ if (($4 ~ /^C[ATC]/)) {{mCH += $5; CH += $6 }}}} END {{print mCH, CH }} '"
        mCG_command = f"awk -v OFS='\t' -v mCG=0 -v CG=0 '{{ if (($4 ~ /^C[G]/)) {{mCG += $5; CG += $6 }}}} END {{print mCG, CG }} '" 

    mcg_fraction, mch_fraction = get_stats(output_path, mCH_command, mCG_command)

    return mcg_fraction, mch_fraction
    
def control_stats(output_path, control_chrom):

    if control_chrom == "chrL":
        mCH_command = "awk -v OFS='\t' -v mCH=0 -v CH=0 '{ if (($4 ~ /^[ATCG]C[ATC]/) && ($4 != \"CCAGG\") && ($4 != \"CCTGG\")) {mCH += $5; CH += $6 }} END {print mCH, CH } '"
        mCG_command = "awk -v OFS='\t' -v mCG=0 -v CG=0 '{ if (($4 ~ /^[ATCG]C[G]/) && ($4 != \"CCAGG\") && ($4 != \"CCTGG\")) {mCG += $5; CG += $6 }} END {print mCG, CG } '"
    else:
        mCH_command = f"awk -v OFS='\t' -v mCH=0 -v CH=0 '{{ if (($4 ~ /^C[ATC]/)) {{mCH += $5; CH += $6 }}}} END {{print mCH, CH }} '"
        mCG_command = f"awk -v OFS='\t' -v mCG=0 -v CG=0 '{{ if (($4 ~ /^C[G]/)) {{mCG += $5; CG += $6 }}}} END {{print mCG, CG }} '" 

    mcg_fraction, mch_fraction = get_stats(output_path, mCH_command, mCG_command)

    return mcg_fraction, mch_fraction

def _read_faidx(faidx_path):
    """
    Read fadix of reference fasta file.

    samtools fadix ref.fa
    """
    return pd.read_csv(
        faidx_path,
        index_col=0,
        header=None,
        sep="\t",
        names=["NAME", "LENGTH", "OFFSET", "LINEBASES", "LINEWIDTH"],
    )


def _get_chromosome_sequence_upper(fasta_path, fai_df, query_chrom):
    """Read a whole chromosome sequence into memory."""
    chrom_pointer = fai_df.loc[query_chrom, "OFFSET"]
    tail = fai_df.loc[query_chrom, "LINEBASES"] - fai_df.loc[query_chrom, "LINEWIDTH"]
    seq = ""
    with open(fasta_path) as f:
        f.seek(chrom_pointer)
        for line in f:
            if line[0] == ">":
                break
            seq += line[:tail]  # trim \n
    return seq.upper()


def _get_bam_chrom_index(bam_path):
    result = subprocess.run(["samtools", "idxstats", bam_path], stdout=subprocess.PIPE, encoding="utf8").stdout

    chrom_set = set()
    for line in result.split("\n"):
        chrom = line.split("\t")[0]
        if chrom not in ["", "*"]:
            chrom_set.add(chrom)
    return pd.Index(chrom_set)


def _bam_to_allc_worker(
    bam_path,
    reference_fasta,
    fai_df,
    output_path,
    region=None,
    num_upstr_bases=0,
    num_downstr_bases=2,
    buffer_line_number=100000,
    min_mapq=0,
    min_base_quality=1,
    compress_level=5,
    tabix=True,
    save_count_df=False,
):
    """None parallel bam_to_allc worker function, call by bam_to_allc."""
    # mpileup
    if region is None:
        mpileup_cmd = f"samtools mpileup -Q {min_base_quality} " f"-q {min_mapq} --ff QCFAIL,DUP -B -f {reference_fasta} {bam_path}"
        pipes = subprocess.Popen(
            shlex.split(mpileup_cmd),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )
    else:
        bam_handle = open_bam(
            bam_path,
            region=region,
            mode="r",
            include_header=True,
            samtools_parms_str=None,
        )
        mpileup_cmd = f"samtools mpileup -Q {min_base_quality} " f"-q {min_mapq} --ff QCFAIL,DUP -B -f {reference_fasta} -"
        pipes = subprocess.Popen(
            shlex.split(mpileup_cmd),
            stdin=bam_handle.file,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            universal_newlines=True,
        )

    result_handle = pipes.stdout

    output_file_handler = open_allc(output_path, mode="w", compresslevel=compress_level)

    # initialize variables
    complement = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    mc_sites = {"C", "G"}
    context_len = num_upstr_bases + 1 + num_downstr_bases
    cur_chrom = ""
    line_counts = 0
    total_line = 0
    out = ""
    seq = None  # whole cur_chrom seq
    chr_out_pos_list = []
    cur_out_pos = 0
    cov_dict = collections.defaultdict(int)  # context: cov_total
    mc_dict = collections.defaultdict(int)  # context: mc_total
    
    # process mpileup result
    for line in result_handle:
        total_line += 1
        fields = line.split("\t")
        fields[2] = fields[2].upper()
        # if chrom changed, read whole chrom seq from fasta
        if fields[0] != cur_chrom:
            cur_chrom = fields[0]
            chr_out_pos_list.append((cur_chrom, str(cur_out_pos)))
            # get seq for cur_chrom
            seq = _get_chromosome_sequence_upper(reference_fasta, fai_df, cur_chrom)

        if fields[2] not in mc_sites:
            continue

        # deal with indels
        read_bases = fields[4]
        incons_basecalls = read_bases.count("+") + read_bases.count("-")
        if incons_basecalls > 0:
            read_bases_no_indel = ""
            index = 0
            prev_index = 0
            while index < len(read_bases):
                if read_bases[index] == "+" or read_bases[index] == "-":
                    # get insert size
                    indel_size = ""
                    ind = index + 1
                    while True:
                        try:
                            int(read_bases[ind])
                            indel_size += read_bases[ind]
                            ind += 1
                        except Exception:
                            break
                    try:
                        # sometimes +/- does not follow by a number and
                        # it should be ignored
                        indel_size = int(indel_size)
                    except Exception:
                        index += 1
                        continue
                    read_bases_no_indel += read_bases[prev_index:index]
                    index = ind + indel_size
                    prev_index = index
                else:
                    index += 1
            read_bases_no_indel += read_bases[prev_index:index]
            fields[4] = read_bases_no_indel

        # count converted and unconverted bases
        if fields[2] == "C":
            # mpileup pos is 1-based, turn into 0 based
            pos = int(fields[1]) - 1
            try:
                context = seq[(pos - num_upstr_bases) : (pos + num_downstr_bases + 1)]
            except Exception:  # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(".")
            converted_c = fields[4].count("T")
            cov = unconverted_c + converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                data = (
                    "\t".join(
                        [
                            cur_chrom,
                            str(pos + 1),
                            "+",
                            context,
                            str(unconverted_c),
                            str(cov),
                            "1",
                        ]
                    )
                    + "\n"
                )
                cov_dict[context] += cov
                mc_dict[context] += unconverted_c
                out += data
                cur_out_pos += len(data)

        elif fields[2] == "G":
            pos = int(fields[1]) - 1
            try:
                context = "".join(
                    [
                        complement[base]
                        for base in reversed(seq[(pos - num_downstr_bases) : (pos + num_upstr_bases + 1)])
                    ]
                )
            except Exception:  # complete context is not available, skip
                continue
            unconverted_c = fields[4].count(",")
            converted_c = fields[4].count("a")
            cov = unconverted_c + converted_c
            if cov > 0 and len(context) == context_len:
                line_counts += 1
                data = (
                    "\t".join(
                        [
                            cur_chrom,
                            str(pos + 1),  # ALLC pos is 1-based
                            "-",
                            context,
                            str(unconverted_c),
                            str(cov),
                            "1",
                        ]
                    )
                    + "\n"
                )
                cov_dict[context] += cov
                mc_dict[context] += unconverted_c
                out += data
                cur_out_pos += len(data)

        if line_counts > buffer_line_number:
            output_file_handler.write(out)
            line_counts = 0
            out = ""

    if line_counts > 0:
        output_file_handler.write(out)
    result_handle.close()
    output_file_handler.close()

    if tabix:
        subprocess.run(shlex.split(f"tabix -b 2 -e 2 -s 1 {output_path}"), check=True)

    count_df = pd.DataFrame({"mc": mc_dict, "cov": cov_dict})
    count_df["mc_rate"] = count_df["mc"] / count_df["cov"]

    total_genome_length = fai_df["LENGTH"].sum()
    count_df["genome_cov"] = total_line / total_genome_length

    if save_count_df:
        count_df.to_csv(output_path + ".count.csv")
        return None
    else:
        return count_df


def _aggregate_count_df(count_dfs):
    total_df = pd.concat(count_dfs)
    total_df = total_df.groupby(total_df.index).sum()
    total_df["mc_rate"] = total_df["mc"] / total_df["cov"]
    total_df["mc"] = total_df["mc"].astype(int)
    total_df["cov"] = total_df["cov"].astype(int)
    return total_df


def bam_to_allc(
    bam_path,
    reference_fasta,
    out_prefix,
    control_chrom=None,
    num_upstr_bases=0,
    num_downstr_bases=2,
    min_mapq=10,
    min_base_quality=20,
    compress_level=5,
    save_count_df=False,
    convert_bam_strandness=True,
):
    """\
    Generate 1 ALLC file from 1 position sorted BAM file via samtools mpileup.

    Parameters
    ----------
    bam_path
        Path to 1 position sorted BAM file
    reference_fasta
        {reference_fasta_doc}
    out_prefix
        Prefix for output ALLC file
    num_upstr_bases
        Number of upstream base(s) of the C base to include in ALLC context column,
        usually use 0 for normal BS-seq, 1 for NOMe-seq.
    num_downstr_bases
        Number of downstream base(s) of the C base to include in ALLC context column,
        usually use 2 for both BS-seq and NOMe-seq.
    min_mapq
        Minimum MAPQ for a read being considered, samtools mpileup parameter, see samtools documentation.
    min_base_quality
        Minimum base quality for a base being considered, samtools mpileup parameter,
        see samtools documentation.
    compress_level
        {compress_level_doc}
    save_count_df
        If true, save an ALLC context count table next to ALLC file.
    convert_bam_strandness
        {convert_bam_strandness_doc}

    Returns
    -------
    count_df
        a pandas.DataFrame for overall mC and cov count separated by mC context.
    """
    buffer_line_number = 100000
    tabix = True

    stats = {}
    genome_cov_dup = compute_genome_coverage(bam_path, min_mapq, min_base_quality, keep_dup=True)
    genome_cov_dedup = compute_genome_coverage(bam_path, min_mapq, min_base_quality, keep_dup=False)
    mapped_nuc_dup = compute_mapped_nucleotides(bam_path, min_mapq, min_base_quality, keep_dup=True)
    mapped_nuc_dedup = compute_mapped_nucleotides(bam_path, min_mapq, min_base_quality, keep_dup=False)

    stats.update(
        {f"allc_input_genome_coverage_dup" : genome_cov_dup,
         f"allc_input_genome_coverage_dedup" : genome_cov_dedup,
         f"allc_input_mapped_bp_dup" : mapped_nuc_dup,
         f"allc_input_mapped_bp_dedup" : mapped_nuc_dedup}
    )
    
    # Check fasta index
    if not pathlib.Path(reference_fasta).exists():
        raise FileNotFoundError(f"Reference fasta not found at {reference_fasta}.")
    if not pathlib.Path(reference_fasta + ".fai").exists():
        raise FileNotFoundError("Reference fasta not indexed. Use samtools faidx to index it and run again.")
    fai_df = _read_faidx(pathlib.Path(reference_fasta + ".fai"))

    if convert_bam_strandness:
        temp_bam_path = f"{out_prefix}.temp.bam"
        _convert_bam_strandness(bam_path, temp_bam_path)
        bam_path = temp_bam_path

    if not pathlib.Path(bam_path + ".bai").exists():
        subprocess.check_call(shlex.split("samtools index " + bam_path))

    # check chromosome between BAM and FASTA
    # samtools have a bug when chromosome not match...
    bam_chroms_index = _get_bam_chrom_index(bam_path)
    unknown_chroms = [i for i in bam_chroms_index if i not in fai_df.index]
    if len(unknown_chroms) != 0:
        unknown_chroms = " ".join(unknown_chroms)
        raise IndexError(
            f"BAM file contain unknown chromosomes: {unknown_chroms}\n"
            "Make sure you use the same genome FASTA file for mapping and bam-to-allc."
        )
    
    # Output path
    output_path = f"{out_prefix}.allc.tsv.gz"

    result = _bam_to_allc_worker(
        bam_path,
        reference_fasta,
        fai_df,
        output_path,
        region=None,
        num_upstr_bases=num_upstr_bases,
        num_downstr_bases=num_downstr_bases,
        buffer_line_number=buffer_line_number,
        min_mapq=min_mapq,
        min_base_quality=min_base_quality,
        compress_level=compress_level,
        tabix=tabix,
        save_count_df=save_count_df,
    )

    if control_chrom:
        control_output_path = f"{out_prefix}_control.allc.tsv.gz"
        
        result = _bam_to_allc_worker(
            bam_path,
            reference_fasta,
            fai_df,
            control_output_path,
            region=control_chrom,
            num_upstr_bases=1,
            num_downstr_bases=3,
            buffer_line_number=buffer_line_number,
            min_mapq=min_mapq,
            min_base_quality=min_base_quality,
            compress_level=compress_level,
            tabix=False,
            save_count_df=False,
        )

        control_mcg_fraction, control_mch_fraction = control_stats(control_output_path, control_chrom)

        stats.update(
            {f"mCG/CG_{control_chrom}" : control_mcg_fraction,
             f"mCH/CH_{control_chrom}" : control_mch_fraction}
        )

        subprocess.check_call(["rm", "-f", f"{control_output_path}"])

    global_mcg_fraction, global_mch_fraction = global_stats(output_path, control_chrom)
    
    stats.update(
        {f"mCG/CG_global" : global_mcg_fraction,
         f"mCH/CH_global" : global_mch_fraction}
    )

    stats_path = f"{out_prefix}_methylation_stats.txt"
    
    stats_df = pd.DataFrame.from_dict(stats, orient="index").T
    stats_df.to_csv(stats_path, index=False, sep="\t")

    # clean up temp bam
    if convert_bam_strandness:
        # this bam path is the temp file path
        subprocess.check_call(["rm", "-f", bam_path])
        subprocess.check_call(["rm", "-f", f"{bam_path}.bai"])

    

    return result
