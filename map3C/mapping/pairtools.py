"""
This file contains functions (some modified) from pairtools https://github.com/open2c/pairtools (v1.0.3)
"""

UNMAPPED_CHROM = "!"
UNMAPPED_POS = 0
UNMAPPED_STRAND = "-"

def cigar_dict(cigartuples, cigarstring):
    """Parse cigar tuples reported as cigartuples of pysam read entry.
    Reports alignment span, clipped nucleotides and more.
    See https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.cigartuples
    """
    matched_bp = 0
    algn_ref_span = 0
    algn_read_span = 0
    read_len = 0
    clip5_ref = 0
    clip3_ref = 0

    if cigartuples is not None:
        for operation, length in cigartuples:
            if operation == 0:  # M, match
                matched_bp += length
                algn_ref_span += length
                algn_read_span += length
                read_len += length
            elif operation == 1:  # I, insertion
                algn_read_span += length
                read_len += length
            elif operation == 2:  # D, deletion
                algn_ref_span += length
            elif (
                    operation == 4 or operation == 5
            ):  # S and H, soft clip and hard clip, respectively
                read_len += length
                if matched_bp == 0:
                    clip5_ref = length
                else:
                    clip3_ref = length

    return {
        "clip5_ref": clip5_ref,
        "clip3_ref": clip3_ref,
        "cigar": cigarstring,
        "algn_ref_span": algn_ref_span,
        "algn_read_span": algn_read_span,
        "read_len": read_len,
        "matched_bp": matched_bp
    }

def parse_read(read, 
               span,
               idx,
               report_3_alignment_end=True,
               min_mapq=30):

    cigarstring = read.cigarstring
    cigartuples = read.cigartuples

    mate = read.query_name.split("_")[-1]
    flag = read.flag
    
    is_unique = read.mapping_quality >= min_mapq
    is_mapped = (flag & 0x04) == 0
    is_linear = not read.has_tag("SA")
    cigar = cigar_dict(cigartuples, cigarstring)
    
    if is_mapped:
        if (flag & 0x10) == 0:
            strand = "+"
            dist_to_5 = cigar["clip5_ref"]
            dist_to_3 = cigar["clip3_ref"]
        else:
            strand = "-"
            dist_to_5 = cigar["clip3_ref"]
            dist_to_3 = cigar["clip5_ref"]

        if is_unique:
            chrom = read.reference_name
            if strand == "+":
                # Note that pysam output is zero-based, thus add +1:
                pos5 = read.reference_start + 1
                pos3 = read.reference_start + cigar["algn_ref_span"]
            else:
                pos5 = read.reference_start + cigar["algn_ref_span"]
                # Note that pysam output is zero-based, thus add +1:
                pos3 = read.reference_start + 1
        else:
            chrom = UNMAPPED_CHROM
            strand = UNMAPPED_STRAND
            pos5 = UNMAPPED_POS
            pos3 = UNMAPPED_POS
    else:
        chrom = UNMAPPED_CHROM
        strand = UNMAPPED_STRAND
        pos5 = UNMAPPED_POS
        pos3 = UNMAPPED_POS

        dist_to_5 = 0
        dist_to_3 = 0

    algn = {
        "chrom": chrom,
        "pos5": pos5,
        "pos3": pos3,
        "strand": strand,
        "mapq": read.mapping_quality,
        "is_mapped": is_mapped,
        "is_unique": is_unique,
        "is_linear": is_linear,
        "dist_to_5": dist_to_5,
        "dist_to_3": dist_to_3,
        "type": ("N" if not is_mapped else ("M" if not is_unique else "U")),
        "span":span,
        "mate":mate,
        "idx": idx,
        "read" : read,
        "reference_start": read.reference_start,
        "reference_end": read.reference_end
    }

    algn.update(cigar)

    algn["pos"] = algn["pos3"] if report_3_alignment_end else algn["pos5"]

    return algn
    

### Additional functions for complex walks rescue ###
def partial_overlap(algn1, algn2, max_insert_size=500, dedup_max_mismatch=5):
    """
    Two ends of alignments overlap if:
     1) they are from the same chromosome,
     2) map in the opposite directions,
     3) the distance between the outer ends of the two alignments is below the specified max_insert_size,
     4) the distance between the outer ends of the two alignments is above the maximum alignment size.
    (4) guarantees that the alignments point towards each other on the chromosomes.

    Allowed mismatch between intramolecular alignments to detect readthrough duplicates.

    Return: 1 if the alignments overlap or both have troubles with unique mapping,
            0 if they do not overlap or if we do not have enough information
            (e.g. only one of the alignments have troubles with being mapped).
    """

    # Alignments with no match or with multiple matches are counted as overlaps
    if not (algn1["is_mapped"] and algn1["is_unique"]):
        if not (algn2["is_mapped"] and algn2["is_unique"]):
            return 1

    # We assume that successful alignment cannot be an overlap with unmapped or multi-mapped region
    if not (algn1["is_mapped"] and algn1["is_unique"]):
        return 0
    if not (algn2["is_mapped"] and algn2["is_unique"]):
        return 0

    # Both alignments are mapped and unique
    do_overlap = True

    do_overlap &= algn1["chrom"] == algn2["chrom"]
    do_overlap &= algn1["strand"] != algn2["strand"]

    if algn1["strand"] == "+":
        min_algn_size = max(
            algn1["pos3"] - algn1["pos5"], algn2["pos5"] - algn2["pos3"]
        )
        distance_outer_ends = algn2["pos5"] - algn1["pos5"]
    else:
        min_algn_size = max(
            algn1["pos5"] - algn1["pos3"], algn2["pos3"] - algn2["pos5"]
        )
        distance_outer_ends = algn1["pos5"] - algn2["pos5"]

    do_overlap &= distance_outer_ends <= max_insert_size + dedup_max_mismatch
    do_overlap &= distance_outer_ends >= min_algn_size - dedup_max_mismatch

    if do_overlap:
        return 1
    return 0


def pairs_overlap(algns1, algns2, dedup_max_mismatch=3):
    """
    We assume algns1 originate from left read, and algns2 originate from right read:
    left read:                             right read:
    ---------------------------->     <----------------------------
                algns1                             algns2
    5------------3_5------------3     3------------5_3------------5'
    left_5'-algn    left_3'-algn      right_3'-algn   right_5'-algn

    Two pairs of alignments overlap if:
    1) chromosomes/mapping/strand of left_5'-algn and right_3'-algn are the same,
    2) chromosomes/mapping/strand of left_3'-algn and right_5'-algn are the same,
    3) pos3 of left_5'-algn is close to pos5 of right_3'-algn (with dedup_max_mismatch), and
    4) pos5 of left_3'-algn is close to pos3 of right_5'-algn.

    Return: 1 of the pairs of alignments overlap, 0 otherwise.
    """
    left5_algn = algns1[0]
    left3_algn = algns1[1]
    right5_algn = algns2[0]
    right3_algn = algns2[1]

    # We assume that successful alignment cannot be an overlap with unmapped or multi-mapped region:
    mapped_left5_algn = left5_algn["is_mapped"] and left5_algn["is_unique"]
    mapped_left3_algn = left3_algn["is_mapped"] and left3_algn["is_unique"]
    mapped_right5_algn = right5_algn["is_mapped"] and right5_algn["is_unique"]
    mapped_right3_algn = right3_algn["is_mapped"] and right3_algn["is_unique"]

    if not mapped_left5_algn and not mapped_right3_algn:
        left_overlap = True
    elif not mapped_left5_algn and mapped_right3_algn:
        left_overlap = False
    elif mapped_left5_algn and not mapped_right3_algn:
        left_overlap = False
    else:
        left_overlap = True
        left_overlap &= left5_algn["chrom"] == right3_algn["chrom"]
        left_overlap &= left5_algn["strand"] != right3_algn["strand"]

    if not mapped_left3_algn and not mapped_right5_algn:
        right_overlap = True
    elif not mapped_left3_algn and mapped_right5_algn:
        right_overlap = False
    elif mapped_left3_algn and not mapped_right5_algn:
        right_overlap = False
    else:
        right_overlap = True
        right_overlap &= left3_algn["chrom"] == right5_algn["chrom"]
        right_overlap &= left3_algn["strand"] != right5_algn["strand"]

    same_pair = True
    same_pair &= abs(left5_algn["pos3"] - right3_algn["pos5"]) <= dedup_max_mismatch
    same_pair &= abs(left3_algn["pos5"] - right5_algn["pos3"]) <= dedup_max_mismatch

    if left_overlap & right_overlap & same_pair:
        return 1
    else:
        return 0

def flip_position(hic_algn):
    """
    Flip ends of a single alignment
    :param hic_algn: Alignment to be modified
    :return:
    """
    hic_algn = dict(hic_algn)  # overwrite the variable with the copy of dictionary
    hic_algn["pos5"], hic_algn["pos3"] = hic_algn["pos3"], hic_algn["pos5"]
    return hic_algn

def flip_orientation(hic_algn):
    """
    Flip orientation of a single alignment
    :param hic_algn: Alignment to be modified
    :return:
    """
    hic_algn = dict(hic_algn)  # overwrite the variable with the copy of dictionary
    hic_algn["strand"] = "+" if (hic_algn["strand"] == "-") else "-"
    return hic_algn

def format_pair(
    hic_algn1,
    hic_algn2,
    pair_index,
    report_position="outer",
    report_orientation="pair",
    algn1_pos5=None,
    algn1_pos3=None,
    algn2_pos5=None,
    algn2_pos3=None,
):
    """
    Return a triplet: pair of formatted alignments and pair_index in a walk

    :param hic_algn1: Left alignment forming a pair
    :param hic_algn2: Right alignment forming a pair
    :param algns1: All left read alignments for formal reporting
    :param algns2: All right read alignments for formal reporting
    :param pair_index: Index of the pair
    :param algn1_pos5: Replace reported 5'-position of the alignment 1 with this value
    :param algn1_pos3: Replace reported 3'-position of the alignment 1 with this value
    :param algn2_pos5: Replace reported 5'-position of the alignment 2 with this value
    :param algn2_pos3: Replace reported 3'-position of the alignment 2 with this value

    """
    # Make sure the original data is not modified:
    hic_algn1, hic_algn2 = dict(hic_algn1), dict(hic_algn2)

    # Adjust the 5' and 3'-ends:
    hic_algn1["pos5"] = algn1_pos5 if not algn1_pos5 is None else hic_algn1["pos5"]
    hic_algn1["pos3"] = algn1_pos3 if not algn1_pos3 is None else hic_algn1["pos3"]
    hic_algn2["pos5"] = algn2_pos5 if not algn2_pos5 is None else hic_algn2["pos5"]
    hic_algn2["pos3"] = algn2_pos3 if not algn2_pos3 is None else hic_algn2["pos3"]

    hic_algn1["type"] = (
        "N"
        if not hic_algn1["is_mapped"]
        else "M"
        if not hic_algn1["is_unique"]
        else "U"
    )

    hic_algn2["type"] = (
        "N"
        if not hic_algn2["is_mapped"]
        else "M"
        if not hic_algn2["is_unique"]
        else "U"
    )

    # Change orientation and positioning of pair for reporting:
    # AVAILABLE_REPORT_POSITION    = ["outer", "pair", "read", "walk"]
    # AVAILABLE_REPORT_ORIENTATION = ["pair", "pair", "read", "walk"]
    pair_type = pair_index[1]

    if report_orientation == "read":
        pass
    elif report_orientation == "walk":
        if pair_type == "R2":
            hic_algn1 = flip_orientation(hic_algn1)
            hic_algn2 = flip_orientation(hic_algn2)
        elif pair_type == "R1-2":
            hic_algn2 = flip_orientation(hic_algn2)
    elif report_orientation == "pair":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn2 = flip_orientation(hic_algn2)
        elif pair_type == "R2":
            hic_algn1 = flip_orientation(hic_algn1)
    elif report_orientation == "junction":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn1 = flip_orientation(hic_algn1)
        elif pair_type == "R2":
            hic_algn2 = flip_orientation(hic_algn2)
        else:
            hic_algn1 = flip_orientation(hic_algn1)
            hic_algn2 = flip_orientation(hic_algn2)

    if report_position == "read":
        pass
    elif report_position == "walk":
        if pair_type == "R2":
            hic_algn1 = flip_position(hic_algn1)
            hic_algn2 = flip_position(hic_algn2)
        elif pair_type == "R1-2":
            hic_algn2 = flip_position(hic_algn2)
    elif report_position == "outer":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn2 = flip_position(hic_algn2)
        elif pair_type == "R2":
            hic_algn1 = flip_position(hic_algn1)
    elif report_position == "junction":
        if pair_type == "R1" or pair_type == "R1&R2":
            hic_algn1 = flip_position(hic_algn1)
        elif pair_type == "R2":
            hic_algn2 = flip_position(hic_algn2)
        else:
            hic_algn1 = flip_position(hic_algn1)
            hic_algn2 = flip_position(hic_algn2)

    return [hic_algn1, hic_algn2, pair_index]




def parse_complex_walk(
    algns1,
    algns2,
    max_insert_size,
    report_position="outer",
    report_orientation="pair",
    dedup_max_mismatch=3,
):
    """
    Parse a set of ligations that appear as a complex walk.
    This procedure is equivalent to intramolecular deduplication that preserved pair order in a walk.

    :param algns1: List of sequential lefts alignments
    :param algns2: List of sequential right alignments
    :param max_insert_size: maximum insert size when searching for overlapping ends of R1 and R2
    :param report_position: one of "outer", "junction", "read", "walk"; sets pos5 and pos3
    :param report_orientation: one of "pair", "junction", "read", "walk"; sets strand
    :param dedup_max_mismatch: allowed mismatch between intramolecular alignments to detect readthrough duplicates
    :param expand: perform combinatorial expansion of pairs or not
    :param max_expansion_depth: maximum depth (number of segments separating pair). All by default.

    :return: iterator with parsed pairs

    **Intramolecular deduplication**

     Forward read (left):                       right read (right):
    5'------------------------->3'     3'<--------------------------5'
             algns1                              algns2
    <5---3><5---3><5---3><5---3>        <3---5><3---5><3---5><3---5>
       l0     l1    l2     l3              r3     r2     r1    r0

    Alignment - bwa mem reported hit or alignment after gaps conversion.
    Left and right alignments (algns1: [l0, l1, l2, l3], algns2: [r0, r1, r2, r3])
    - alignments on left and right reads reported from 5' to 3' orientation.

    Intramolecular deduplication consists of two steps:
    I. iterative search of overlapping alignment pairs (aka overlap),
    II. if no overlaps or search not possible (less than 2 alignments on either sides),
    search for overlap of end alignments (aka partial overlap).
    III. report pairs before the overlap, deduplicated pairs of overlap and pairs after that.

    Iterative search of overlap is in fact scanning of the right read pairs for the hit
    with the 3'-most pair of the left read:
        1. Initialize.
            Start from 3' of left and right reads. Set `current_left_pair` and `current_right_pair` pointers
        2. Initial compare.
            Compare pairs l2-l3 and r3-r2 by `pairs_overlap`.
                If successful, we found the overlap, go to reporting.
                If unsuccessful, continue search.
        3. Increment.
            Shift `current_right_pair` pointer by one (e.g., take the pair r2-r1).
        4. Check.
            Check that this pair can form a potential overlap with left alignments:
            the number of pairs downstream from l2-l3 on left read should not be less than
            the number of pairs upstream from r2-r1 on right read.
                If overlap cannot be formed, no other overlap in this complex walk is possible, safely exit.
                If the potential overlap can be formed, continue comparison.
        5. Compare.
            Compare the current pair of pairs on left and right reads.
                If comparison fails, go to step 3.
                If comparison is successful, go to 6.
        6. Verify.
            Check that downstream pairs on the left read overlap with the upstream pairs on the right read.
                 If yes, exit.
                 If not, we do not have an overlap, go to step 3.
    """

    AVAILABLE_REPORT_POSITION = ["outer", "junction", "read", "walk"]
    assert report_position in AVAILABLE_REPORT_POSITION, (
        f"Cannot report position {report_position}, as it is not implemented"
        f'Available choices are: {", ".join(AVAILABLE_REPORT_POSITION)}'
    )

    AVAILABLE_REPORT_ORIENTATION = ["pair", "junction", "read", "walk"]
    assert report_orientation in AVAILABLE_REPORT_ORIENTATION, (
        f"Cannot report orientation {report_orientation}, as it is not implemented"
        f'Available choices are: {", ".join(AVAILABLE_REPORT_ORIENTATION)}'
    )

    output_pairs = []

    # Initialize (step 1).
    n_algns1 = len(algns1)
    n_algns2 = len(algns2)
    current_left_pair = current_right_pair = 1
    remaining_left_pairs = (
        n_algns1 - 1
    )  # Number of possible pairs remaining on left read
    remaining_right_pairs = (
        n_algns2 - 1
    )  # Number of possible pairs remaining on right read
    checked_right_pairs = (
        0  # Number of checked pairs on right read (from the end of read)
    )
    is_overlap = False

    # I. Iterative search of overlap, at least two alignments on each side:
    if (n_algns1 >= 2) and (n_algns2 >= 2):
        #print(">=2 on each side")
        # Iteration includes check (step 4):
        while (remaining_left_pairs > checked_right_pairs) and (
            remaining_right_pairs > 0
        ):
            pair1 = (algns1[-current_left_pair - 1], algns1[-current_left_pair])
            pair2 = (algns2[-current_right_pair - 1], algns2[-current_right_pair])
            # Compare (initial or not, step 2 or 5):
            is_overlap = pairs_overlap(
                pair1, pair2, dedup_max_mismatch=dedup_max_mismatch
            )
            # print(pair1)
            # print(pair2)
            # print(is_overlap)
            # print(remaining_right_pairs)
            # print(remaining_left_pairs)
            if is_overlap:
                last_idx_left_temp = current_left_pair
                last_idx_right_temp = current_right_pair
                checked_right_temp = checked_right_pairs
                # Verify (step 6):
                while is_overlap and (checked_right_temp > 0):
                    last_idx_left_temp += 1
                    last_idx_right_temp -= 1
                    pair1 = (
                        algns1[-last_idx_left_temp - 1],
                        algns1[-last_idx_left_temp],
                    )
                    pair2 = (
                        algns2[-last_idx_right_temp - 1],
                        algns2[-last_idx_right_temp],
                    )
                    is_overlap &= pairs_overlap(
                        pair1, pair2, dedup_max_mismatch=dedup_max_mismatch
                    )
                    checked_right_temp -= 1
                if is_overlap:  # exit
                    current_right_pair += 1
                    break

            # Increment pointers (step 3)
            current_right_pair += 1
            checked_right_pairs += 1
            remaining_right_pairs -= 1

        # No overlap found, roll the current_idx_right back to the initial value:
        if not is_overlap:
            #print("No overlap found")
            current_right_pair = 1

        #print("overlap" if is_overlap else "no overlap")
    #else:
        #print("not >=2 on each side")

    # II. Search of partial overlap if there are less than 2 alignments at either sides, or no overlaps found
    if current_right_pair == 1:
        
        last_reported_alignment_left = last_reported_alignment_right = 1
        if partial_overlap(
            algns1[-1],
            algns2[-1],
            max_insert_size=max_insert_size,
            dedup_max_mismatch=dedup_max_mismatch,
        ):
            #print("Partial overlap")
            if (
                n_algns1 >= 2
            ):  # single alignment on right read and multiple alignments on left
                #print("algns1 greater than 2")
                pair_index = (len(algns1) - 1, "R1")
                output_pairs.append(
                    format_pair(
                        algns1[-2],
                        algns1[-1],
                        pair_index=pair_index,
                        algn2_pos3=algns2[-1]["pos5"],
                        report_position=report_position,
                        report_orientation=report_orientation,
                    )
                )
                last_reported_alignment_left = 2  # set the pointer for reporting

            if (
                n_algns2 >= 2
            ):  # single alignment on left read and multiple alignments on right
                #print("algns2  greater than 2")
                pair_index = (len(algns1), "R2")
                output_pairs.append(
                    format_pair(
                        algns2[-1],
                        algns2[-2],
                        pair_index=pair_index,
                        algn1_pos3=algns1[-1]["pos5"],
                        report_position=report_position,
                        report_orientation=report_orientation,
                    )
                )
                last_reported_alignment_right = 2  # set the pointer for reporting

            # Note that if n_algns1==n_algns2==1 and alignments overlap, then we don't need to check,
            # it's a non-ligated DNA fragment that we don't report.

        else:  # end alignments do not overlap, report regular pair:
            #print("No partial overlap")
            pair_index = (len(algns1), "R1-2")
            output_pairs.append(
                format_pair(
                    algns1[-1],
                    algns2[-1],
                    pair_index=pair_index,
                    report_position=report_position,
                    report_orientation=report_orientation,
                )
            )

    else:  # there was an overlap, set some pointers:
        last_reported_alignment_left = (
            last_reported_alignment_right
        ) = current_right_pair

    # overlap same, no overlap, can also be same 
    #print(last_reported_alignment_left)
    #print(current_right_pair)

    # III. Report all remaining alignments.
    # Report all unique alignments on left read (sequential):
    for i in range(0, n_algns1 - last_reported_alignment_left):
        pair_index = (i + 1, "R1")
        output_pairs.append(
            format_pair(
                algns1[i],
                algns1[i + 1],
                pair_index=pair_index,
                report_position=report_position,
                report_orientation=report_orientation,
            )
        )

    # Report the pairs where both left alignments overlap right:
    for i_overlapping in range(current_right_pair - 1):
        idx_left = n_algns1 - current_right_pair + i_overlapping
        idx_right = n_algns2 - 1 - i_overlapping
        pair_index = (idx_left + 1, "R1&2")
        output_pairs.append(
            format_pair(
                algns1[idx_left],
                algns1[idx_left + 1],
                pair_index=pair_index,
                algn2_pos3=algns2[idx_right - 1]["pos5"],
                report_position=report_position,
                report_orientation=report_orientation,
            )
        )

    # Report all the sequential chimeric pairs in the right read, but not the overlap:
    #reporting_order = range(
    #    0, min(current_right_pair, n_algns2 - last_reported_alignment_right)
    #)

    reporting_order = range(
        0, n_algns2 - max(current_right_pair, last_reported_alignment_right)
    )
    #print(reporting_order)
    for i in reporting_order:
        # Determine the pair index depending on what is the overlap:
        shift = -1 if current_right_pair > 1 else 0
        pair_index = (
            (
                n_algns1
                #+ min(current_right_pair, n_algns2 - last_reported_alignment_right)
                + (n_algns2 - max(current_right_pair, last_reported_alignment_right))
                - i
                + shift
            ),
            "R2",
        )
        output_pairs.append(
            format_pair(
                algns2[i + 1],
                algns2[i],
                pair_index=pair_index,
                report_position=report_position,
                report_orientation=report_orientation,
            )
        )

    # Sort the pairs according to the pair index:
    output_pairs.sort(key=lambda x: int(x[-1][0]))

    return iter(output_pairs)

def empty_alignment():
    return {
        "chrom": UNMAPPED_CHROM,
        "pos5": UNMAPPED_POS,
        "pos3": UNMAPPED_POS,
        "pos": UNMAPPED_POS,
        "strand": UNMAPPED_STRAND,
        "dist_to_5": 0,
        "dist_to_3": 0,
        "mapq": 0,
        "is_unique": False,
        "is_mapped": False,
        "is_linear": True,
        "cigar": "*",
        "type": "N",
        "span":"NA",
        "mate": "NA",
        "idx" : "NA",
        "reference_start": "NA",
        "reference_end": "NA"
    }

def _convert_gaps_into_alignments(sorted_algns, max_inter_align_gap=20):
    """
    Inplace conversion of gaps longer than max_inter_align_gap into alignments
    """
    if (len(sorted_algns) == 1) and (not sorted_algns[0]["is_mapped"]):
        return

    last_5_pos = 0
    for i in range(len(sorted_algns)):
        algn = sorted_algns[i]
        if algn["dist_to_5"] - last_5_pos > max_inter_align_gap:
            new_algn = empty_alignment()
            new_algn["dist_to_5"] = last_5_pos
            new_algn["algn_read_span"] = algn["dist_to_5"] - last_5_pos
            new_algn["read_len"] = algn["read_len"]
            new_algn["dist_to_3"] = new_algn["read_len"] - algn["dist_to_5"]
            new_algn["gap"] = True # CHANGE

            last_5_pos = algn["dist_to_5"] + algn["algn_read_span"]

            sorted_algns.insert(i, new_algn)
            i += 2
        else:
            last_5_pos = max(last_5_pos, algn["dist_to_5"] + algn["algn_read_span"])
            i += 1


def rescue_walk(algns1, algns2, max_molecule_size):
    # If both sides have one alignment or none, no need to rescue!
    n_algns1 = len(algns1)
    n_algns2 = len(algns2)

    if (n_algns1 <= 1) and (n_algns2 <= 1):
        return None

    # Can rescue only pairs with one chimeric alignment with two parts.
    if not (
        ((n_algns1 == 1) and (n_algns2 == 2)) or ((n_algns1 == 2) and (n_algns2 == 1))
    ):
        return None

    first_read_is_chimeric = n_algns1 > 1
    chim5_algn = algns1[0] if first_read_is_chimeric else algns2[0]
    chim3_algn = algns1[1] if first_read_is_chimeric else algns2[1]
    linear_algn = algns2[0] if first_read_is_chimeric else algns1[0]

    # the linear alignment must be uniquely mapped
    if not (linear_algn["is_mapped"] and linear_algn["is_unique"]):
        return None

    can_rescue = True
    # we automatically rescue chimeric alignments with null and non-unique
    # alignments at the 3' side
    if chim3_algn["is_mapped"] and chim5_algn["is_unique"]:
        # 1) in rescued walks, the 3' alignment of the chimeric alignment must be on
        # the same chromosome as the linear alignment on the opposite side of the
        # molecule
        can_rescue &= chim3_algn["chrom"] == linear_algn["chrom"]

        # 2) in rescued walks, the 3' supplemental alignment of the chimeric
        # alignment and the linear alignment on the opposite side must point
        # towards each other
        can_rescue &= chim3_algn["strand"] != linear_algn["strand"]
        if linear_algn["strand"] == "+":
            can_rescue &= linear_algn["pos5"] < chim3_algn["pos5"]
        else:
            can_rescue &= linear_algn["pos5"] > chim3_algn["pos5"]

        # 3) in single ligations appearing as walks, we can infer the size of
        # the molecule and this size must be smaller than the maximal size of
        # Hi-C molecules after the size selection step of the Hi-C protocol
        if linear_algn["strand"] == "+":
            molecule_size = (
                chim3_algn["pos5"]
                - linear_algn["pos5"]
                + chim3_algn["dist_to_5"]
                + linear_algn["dist_to_5"]
            )
        else:
            molecule_size = (
                linear_algn["pos5"]
                - chim3_algn["pos5"]
                + chim3_algn["dist_to_5"]
                + linear_algn["dist_to_5"]
            )

        can_rescue &= molecule_size <= max_molecule_size

    if can_rescue:
        # changing the type of the 3' alignment on side 1, does not show up in the output:
        if first_read_is_chimeric:

            algns1[1]["type"] = "X"
            algns2[0]["type"] = "R"
            return 1
        # changing the type of the 3' alignment on side 2, does not show up in the output:
        else:
            algns1[0]["type"] = "R"
            algns2[1]["type"] = "X"
            return 2
    else:
        return None

def contact_iter(R1, R2, min_mapq=30, max_molecule_size=750, max_inter_align_gap=None):

    R1_parsed = []
    for i in R1:
        read_data = R1[i]
        read = read_data["trimmed_read"]
        span = read_data["adjusted_span"]
        R1_parsed.append(parse_read(read, span, i, min_mapq))
    
    R2_parsed = []
    for i in R2:
        read_data = R2[i]
        read = read_data["trimmed_read"]
        span = read_data["adjusted_span"]
        R2_parsed.append(parse_read(read, span, i, min_mapq))

    if len(R1_parsed) > 0:
        R1_parsed = sorted(R1_parsed, key=lambda algn: algn["dist_to_5"])
    else:
        R1_parsed = [empty_alignment()]  # Empty alignment dummy
    if len(R2_parsed) > 0:
        R2_parsed = sorted(R2_parsed, key=lambda algn: algn["dist_to_5"])
    else:
        R2_parsed = [empty_alignment()]  # Empty alignment dummy
        
    if max_inter_align_gap is not None:
        _convert_gaps_into_alignments(R1_parsed, max_inter_align_gap)
        _convert_gaps_into_alignments(R2_parsed, max_inter_align_gap)
    
    hic_algn1 = R1_parsed[0]
    hic_algn2 = R2_parsed[0]
    pair_index = (1, "R1-2")

    # Define the type of alignment on each side:
    is_chimeric_1 = len(R1_parsed) > 1
    is_chimeric_2 = len(R2_parsed) > 1

    output_iter = iter([(hic_algn1, hic_algn2, pair_index)])
    
    rule = "mask"
    if is_chimeric_1 or is_chimeric_2:

        rescued_linear_side = rescue_walk(R1_parsed, R2_parsed, max_molecule_size)
        # Walk was rescued as a simple walk:
        if rescued_linear_side == 1 and "gap" not in hic_algn1:
            pair_index = (1, "R1")
            output_iter = iter([(hic_algn1, hic_algn2, pair_index)])
        elif rescued_linear_side == 2 and "gap" not in hic_algn2:
            pair_index = (1, "R2")
            output_iter = iter([(hic_algn1, hic_algn2, pair_index)])
        # Walk is unrescuable:
        else:
            rule = "all"
            output_iter = parse_complex_walk(
                    R1_parsed,
                    R2_parsed,
                    max_molecule_size,
                    report_position="outer",
                    report_orientation="pair")
    return output_iter, rule 
