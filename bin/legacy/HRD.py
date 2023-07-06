##!/usr/bin/python3
import sys

segment_filename = sys.argv[1]
HRD_outfile = open( sys.argv[2], "w")


def read_cnvkit_segments(input_filename, segments):
    segment_file = open(input_filename)
    next(segment_file)
    for line in segment_file:
        columns = line.strip().split("\t")
        if len(columns) < 13:
            continue
        chrom = columns[0]
        start_pos = int(columns[1])
        end_pos = int(columns[2])
        exons = columns[3].split(",")
        nr_exons = len(exons)
        length = end_pos - start_pos + 1
        log2 = float(columns[4])
        baf = 0.5
        if columns[5] != "":
            baf = float(columns[5])
        cn = int(columns[8])
        cn1 = ""
        if columns[9] != "":
            cn1 = int(columns[9])
        cn2 = ""
        if columns[10] != "":
            cn2 = int(columns[10])
        depth = float(columns[11])
        nr_probes = int(columns[12])
        if chrom == "chrY":
            continue
        if chrom not in segments:
            segments[chrom] = []
        segments[chrom].append({
            "start_pos": start_pos, "end_pos": end_pos, "exons": exons, "nr_exons": nr_exons, "length": length,
            "log2": log2, "baf": baf, "cn": cn, "cn1": cn1, "cn2": cn2, "depth": depth, "nr_probes": nr_probes
        })
    return segments


def filter_merge_segments(segments, min_size):
    filtered_merged_segments = {}
    for chrom in segments:
        filtered_merged_segments[chrom] = []
        prev_segment = []
        for segment in segments[chrom]:
            if segment["nr_exons"] <= 1 or segment["nr_probes"] < 20 or segment["length"] < min_size:
                continue
            if prev_segment != []:
                # if prev_segment["cn"] == segment["cn"]:
                if (
                    (prev_segment["cn"] > 2 and segment["cn"] > 2) or
                    (prev_segment["cn"] < 2 and segment["cn"] < 2) or
                    (prev_segment["cn"] == 2 and segment["cn"] == 2)
                ):
                    prev_segment["end_pos"] = segment["end_pos"]
                    prev_segment["length"] = prev_segment["end_pos"] - prev_segment["start_pos"] + 1
                else:
                    filtered_merged_segments[chrom].append(prev_segment)
                    prev_segment = {
                        "start_pos": segment["start_pos"], "end_pos": segment["end_pos"],
                        "cn": segment["cn"], "cn1": segment["cn1"], "cn2": segment["cn2"], "length": segment["length"]
                    }
            else:
                prev_segment = {
                    "start_pos": segment["start_pos"], "end_pos": segment["end_pos"],
                    "cn": segment["cn"], "cn1": segment["cn1"], "cn2": segment["cn2"], "length": segment["length"]
                }
        if prev_segment != []:
            filtered_merged_segments[chrom].append(prev_segment)
    return filtered_merged_segments


def count_LoH_score(filtered_merged_segments):
    LoH_score = 0
    for chrom in filtered_merged_segments:
        nr_segments = 0
        nr_LOH_segments = 0
        for segment in filtered_merged_segments[chrom]:
            nr_segments += 1
            if segment["cn"] != 2 and segment["length"] > 15000000:
                nr_LOH_segments += 1
        '''No LoH_score for LoH of entire chromosome'''
        if nr_segments > 2 or (nr_segments == 2 and nr_LOH_segments == 1):
            LoH_score += nr_LOH_segments
    return LoH_score


def count_LST_score(filtered_merged_segments):
    LST_score = 0
    for chrom in filtered_merged_segments:
        cn_sequence = ""
        for segment in filtered_merged_segments[chrom]:
            length = segment["length"]
            if length < 10000000:
                cn_sequence += "s"
            else:
                cn_sequence += "l"
        i = 0
        max_i = len(cn_sequence) - 1 - 2
        while i <= max_i:
            '''Two adjacent large (l) regions with a gap lower than 3M bp that has lower CN'''
            if cn_sequence[i] == "l" and cn_sequence[i+2] == "l":
                gap = filtered_merged_segments[chrom][i+1]["start_pos"] - filtered_merged_segments[chrom][i]["end_pos"]
                cn1 = filtered_merged_segments[chrom][i]["cn"]
                cn2 = filtered_merged_segments[chrom][i+1]["cn"]
                cn3 = filtered_merged_segments[chrom][i+2]["cn"]
                if gap < 3000000 and cn2 < cn1 and cn2 < cn3:
                    LST_score += 1
            i += 1
    return LST_score


def count_TAI_score(filtered_merged_segments):
    TAI_score = 0
    for chrom in filtered_merged_segments:
        nr_segments = 0
        nr_TAI_segments = 0
        i = 0
        last_segment = len(filtered_merged_segments[chrom]) - 1
        for segment in filtered_merged_segments[chrom]:
            nr_segments += 1
            if (
                (segment["cn"] != 2 or (segment["cn"] == 2 and segment["cn1"] != segment["cn2"]))
                and (i == 0 or i == last_segment)
            ):
                nr_TAI_segments += 1
            i += 1
        '''No TAI_score for CNV of entire chromosome'''
        if nr_segments > 2 or (nr_segments == 2 and nr_TAI_segments == 1):
            TAI_score += nr_TAI_segments
    return TAI_score


'''Read segments from cnvkit'''
'''Skip chr X?'''
segments = {}
segments = read_cnvkit_segments(segment_filename, segments)

'''Filter out small regions and merge regions with same (direction of) cn after filtering'''
filtered_merged_segments = filter_merge_segments(segments, 100000)

'''Calculate LoH score: nr segments > 15Mb with LoH'''
LoH_score = count_LoH_score(filtered_merged_segments)

'''Calculate LST score: nr of adjacent (gap < 3Mb) segments > 10Mb with different cn'''
LST_score = count_LST_score(filtered_merged_segments)

'''Should cn != 2 of entire chromosome be counted as 2,1 or 0 for TAI? Now it is 2.'''
'''Calculate TAI score: nr end segments with cn != 2'''
TAI_score = count_TAI_score(filtered_merged_segments)


HRD_score = LoH_score + LST_score + TAI_score
HRD_outfile.write("HRD_score\tLoH_score\tLST_score\tTAI_score\n")
print ("HRD_score\tLoH_score\tLST_score\tTAI_score\n")
if segments == {}:
    HRD_outfile.write("The HRD score could no be calculated! Probably due to lack of SNPs\n")
else:
    print (str(HRD_score) + "\t" + str(LoH_score) + "\t" + str(LST_score) + "\t" + str(TAI_score) + "\n")
    HRD_outfile.write(str(HRD_score) + "\t" + str(LoH_score) + "\t" + str(LST_score) + "\t" + str(TAI_score) + "\n")
