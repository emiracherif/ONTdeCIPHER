#!/usr/bin/env python
import json
import re
import sys
from collections import OrderedDict

from .vcftagprimersites import read_bed_file

# Alignment_Length_Threshold drops binned reads that are <X% of amplicon length)
Alignment_Length_Threshold = 0.95

# Amplicon_Dropout_Val will report amplicon dropout in any amplicon which has fewer than X reads
Amplicon_Dropout_Val = 50

# Template for the amplicon plot data
amplicon_plot_template = {
    "id": "custom_data_lineplot",
    "section_name": "ARTIC: Amplicon Coverage",
    "description": "This plot summarises the number of reads that were assigned to each amplicon in the primer scheme.\nWe use the align_trim report file from the ARTIC pipeline and group each read by its assigned amplicon.\nIf the length of alignment between read and reference is <{}% of the amplicon length, the read discarded from the coverage plot.\nIf the total number of reads assigned to an amplicon is below {} (red dashed line),\nthe amplicon is marked as dropped out." .format(Alignment_Length_Threshold, Amplicon_Dropout_Val),
    "plot_type": "linegraph",
    "pconfig": {
        "id": "custom_data_linegraph",
        "title": "",
        "categories": "True",
        "yDecimals": "False",
        "xDecimals": "False",
        "ylab": "# reads",
        "xlab": "amplicon",
        "yPlotLines": [{
            "color": "#FF0000",
            "width": 2,
            "dashStyle": "LongDash",
            "label": "Amplicon dropout",
            "value": Amplicon_Dropout_Val
        }]
    },
    "data": {}
}

# Template for the stats table data
amplicon_stats_template = {
    "id": "custom_data_json_table",
    "section_name": "ARTIC: General Stats",
    "description": "A summary of stats from the consensus genome pipeline.",
    "plot_type": "table",
    "pconfig": {
        "id": "custom_data_json_table_table",
        "title": "",
        "min": 0,
        "scale": "RdYlGn-rev",
        "format": "{:,.0f}"
    },
    "data": {}
}

def getSchemeAmplicons(schemeFile):
    """Get the expected amplicon names from the provided scheme.

    Parameters
    ----------
    schemeFile : string
        The filename of the primer scheme
    
    Returns
    -------
    dict
        A dict of amplicon names -> zeroed counter
    """
    amplicons = {}
    primer_scheme = read_bed_file(schemeFile)
    for primer in primer_scheme:
        amplicon = ""
        if primer["direction"] == "+":
            amplicon = primer["Primer_ID"].split("_LEFT")[0]
        else:
            amplicon = primer["Primer_ID"].split("_RIGHT")[0]
        if amplicon not in amplicons:
            amplicons[amplicon] = 0
        amplicons[amplicon] += 1
    named_amplicons = {}
    for amplicon in amplicons:
        if amplicons[amplicon] != 2:
            print("in correct numbers of primer for {}" .format(amplicon), file=sys.stderr)
            raise SystemExit(1)
        named_amplicons[("{}_LEFT_{}_RIGHT" .format(amplicon, amplicon))] = 0
    return named_amplicons

def getAmpliconCounts(amplicons, align_trim_report):
    """Get the read counts per amplicon.

    Parameters
    ----------
    amplicons : list
        Dict of amplicon names found in scheme, linked to a zeroed counter

    align_trim_report: string
        File path to the align_trim report

    Returns
    -------
    dict
        Dict of amplicon names -> populated read counts
    """
    # process the align_trim report
    with open(align_trim_report, "r") as fh:

        # skip the first line (header)
        fh.readline()

        # process each line and add to counts
        for l in fh:
            fields = l.rstrip().split('\t')

            # check read is from a properly paired amplicon
            if int(fields[12]) != 1:
                continue

            # check the read alignment length covers enough of the amplicon
            aLen = int(fields[11]) - int(fields[10])
            rLen = int(fields[2]) - int(fields[1])
            if aLen < (Alignment_Length_Threshold * rLen):
                continue

            # increment the read count for this amplicon
            if fields[3] not in amplicons:
                print("amplicon in align_trim report but not in primer scheme {}" .format(fields[3]), file=sys.stderr)
                raise SystemExit(1)
            amplicons[fields[3]] += 1
    return amplicons

def getVCFreportInfo(vcf_report):
    """Get the read counts per amplicon.

    Parameters
    ----------
    vcf_report: string
        File path to the vcf_report

    Returns
    -------
    dict
        Dict of vcf stats -> values
    """
    # Read vcfcheck report and get important stuff out (NOTE: more to be added in next release)
    stats = dict()
    total_vars = 0
    passed_vars = 0
    with open(vcf_report, "r") as fh:
        for l in fh:
            match = re.search(r'.*\t(\d+)\svariant\srecords\sprocessed', l)
            if match:
                total_vars = int(match.group(1))
            match = re.search(r'.*\t(\d+)\svariant\srecords\spassed\schecks', l)
            if match:
                passed_vars = int(match.group(1))
        stats["# overlap var. fails"] = total_vars - passed_vars
    return stats

def run(args):
    """Collect stats from ARTIC pipeline output and generate files for use by MultiQC.
    """
    # get a list of expected amplicon names
    amplicons = getSchemeAmplicons(args.scheme)

    # open align trim output and count reads per amplicon in scheme
    amplicon_counts = getAmpliconCounts(amplicons, args.align_report)

    # replace amplicon names with ints and count number of dropouts
    dropouts = 0
    amplicon_renamed_counts = dict()
    for amplicon, count in amplicon_counts.items():
        amplicon_renamed_counts[int(amplicon.split('_')[1])] = count
        if count < Amplicon_Dropout_Val:
            dropouts += 1
    
    # add counts to multiqc amplicon plot template
    amplicon_plot_template["data"][args.sample] = amplicon_renamed_counts

    # write the amplicon plot output
    with open("{}.amplicon_plot_data_mqc.json" .format(args.sample), "w") as amplicon_plot_mqc_file:
        json.dump(amplicon_plot_template, amplicon_plot_mqc_file, indent=4, sort_keys=False)
    amplicon_plot_mqc_file.close()

    # add counts to multiqc stats template
    amplicon_stats_template["data"][args.sample] = dict()
    amplicon_stats_template["data"][args.sample]["# low cov. amplicons"] = dropouts

    # parse VCF report if provided and add to the stats template
    if args.vcf_report:
        for stat, value in getVCFreportInfo(args.vcf_report).items():
            amplicon_stats_template["data"][args.sample][stat] = value

    # write the stats output
    with open("{}.amplicon_stats_data_mqc.json" .format(args.sample), "w") as amplicon_stats_mqc_file:
        json.dump(amplicon_stats_template, amplicon_stats_mqc_file, indent=4, sort_keys=False)
    amplicon_stats_mqc_file.close()

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Collect stats from ARTIC pipeline output and generate files for use by MultiQC')
    parser.add_argument('--scheme', required=True, type=str, help='the amplicon scheme used')
    parser.add_argument('--align-report', required=True, type=str, help='the report file from align_trim (*.alignreport.txt')
    parser.add_argument('--vcf-report', required=False, type=str, help='the report file from vcf_check (*.vcfreport.txt')
    parser.add_argument('sample', type=str, help='the sample name')
    args = parser.parse_args()
    run(args)

if __name__ == "__main__":
    main()

