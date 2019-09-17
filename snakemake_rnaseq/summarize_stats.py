from pandas.io.json import json_normalize
import pandas as pd
from ast import literal_eval


parser = argparse.ArgumentParser(description='summarize_stats.py: script to fix rmd files for better rendering',
                                 epilog='For questions or comments, please contact Matt Settles <settles@ucdavis.edu> '
                                        'or Keith Mitchell <kgmitchell@ucdavis.edu\n', add_help=True)
parser.add_argument('-s', '--samples', help="Path to the file 'samples.txt' where each row is the name of a file as "
                                            "follows '01-HTS_Preproc/{SAMPLE}/{SAMPLE}_htsStats.log'")


# parser.add_argument('-s', '--samples', help="Input .md rendered for Rmd to be fixed for better presentation.")


# parser.add_argument('-o', '--output', help="Output .md fixed for better presentation in template.",
#                     type=str)

options = parser.parse_args()


# filenames = ['01-HTS_Preproc/SampleAC1/SampleAC1_htsStats.log', '01-HTS_Preproc/SampleAC2/SampleAC2_htsStats.log', ]

samples = []
with open(options['samples']) as samples_file:
    for line in samples_file:
        samples.append(line.strip('\n'))

filenames = [f'01-HTS_Preproc/{sample}/{sample}_htsStats.log' for sample in samples]
print(filenames)
json_total = ''
apps = ["hts_Stats", "hts_SeqScreener", "hts_SuperDeduper", "hts_AdapterTrimmer", "hts_QWindowTrim", "hts_NTrimmer",
        "hts_CutTrim", ]


# TODO CROSS CHECK EACH FILE FOR THE RIGHT APPS AND SEQUENCE
# TODO add sample names a column
check_apps = ["htsStats1", "htsStats2", "htsSeqScreener1", "htsSeqScreener2", "htsSuperDeduper", "htsAdapterTrimmer",
              "htsQWindowTrim", "htsNTrimmer", "htsCutTrim", ]

# MERGE JSON FILES
for fname in filenames:
    stats_dup = iter(['1', '2', ])  # pre stats and post stats
    screen_dup = iter(['1', '2', ])  # Phix screen and seq screen
    with open(fname) as infile:
        for line in infile:
            if any(app in line for app in apps):
                split_line = line.split('_')

                # RENAME DUPLICATES IF STATS OR SCREEN
                if "hts_Stats" in line:
                    rejoin_line = split_line[:-1] + [next(stats_dup) + '": {\n']
                elif "hts_SeqScreen" in line:
                    rejoin_line = split_line[:-1] + [next(screen_dup) + '": {\n']
                else:
                    rejoin_line = split_line[:-1] + ['": {\n']
                json_total += ''.join(rejoin_line)
            else:
                json_total += line
        json_total = json_total.rstrip()
        json_total += ","

# LOAD JSONS MERGED INTO PANDAS DATAFRAME
data = [i for i in literal_eval(json_total)]
df = pd.DataFrame.from_dict(json_normalize(data), orient='columns')

keys = [
    # RAW STATS
    "Raw_Reads", "Raw_Bp", "Raw_R1_PercentQ30", "Raw_R2_PercentQ30", "Raw_Percent_CG", "Raw_Chars_N",
    # PHIX SCREEN
    "PhiX_IN", "PhiX_Discard", "PhiX_Percent_Discard",
    # rRNA SCREEN
    "rRNA_IN", "rRNA_Identified", "rRNA_Percent_Identified",
    # SUPER DEDUPER
    "SuperDeduper_IN", "SuperDeduper_Ignored", "SuperDeduper_Duplicate", "SuperDeduper_Percent_Duplicate",
    # ADAPTER TRIM
    "AdapterTrimmed_IN", "AdapterTrimmed_Reads", "AdapterTrimmed_Percent_Reads", "AdapterTrimmed_BP",
    # QWINDOW TRIM
    "QwindowTrimmed_IN", "QwindowTrimmed_R1_LeftBpTrim", "QwindowTrimmed_R1_RightBpTrim", "QwindowTrimmed_R2_LeftBpTrim", "QwindowTrimmed_R2_RightBpTrim", "QwindowTrimmed_Discard",
    # NCHAR
    "NcharTrimmed_IN", "NcharTrimmed_R1_LeftBpTrim", "NcharTrimmed_R1_RightBpTrim", "NcharTrimmed_R2_LeftBpTrim", "NcharTrimmed_R2_RightBpTrim", "NcharTrimmed_Discard",
    # CUT TRIM
    "MinLen_IN", "MinLen_Discard",
    # POST STATS
    "Proc_Reads", "Proc_Bp", "Proc_R1_PercentQ30", "Proc_R2_PercentQ30", "Proc_Percent_CG", "Proc_Chars_N",
    # FINAL
    "Final_Percent_Read", "Final_Percent_Bp"
]


# RAW STATS
df['Raw_Reads'] = df['htsStats1.totalFragmentsInput']
df['Raw_Bp'] = df['htsStats1.Paired_end.R1_bpLen'] + df['htsStats1.Paired_end.R2_bpLen']
df['Raw_R1_PercentQ30'] = (df['htsStats1.Paired_end.R1_bQ30']/df['htsStats1.Paired_end.R1_bpLen']) * 100
df['Raw_R2_PercentQ30'] = (df['htsStats1.Paired_end.R2_bQ30']/df['htsStats1.Paired_end.R2_bpLen']) * 100
df['Raw_Percent_CG'] = ((df['htsStats1.Base_composition.C']+df['htsStats1.Base_composition.G'])/df['Raw_Bp']) * 100
df['Raw_Chars_N'] = df['htsStats1.Base_composition.N']

# PHIX SCREEN
df['PhiX_IN'] = df['htsSeqScreener1.totalFragmentsInput']
df['PhiX_Discard'] = df['htsSeqScreener1.Paired_end.PE_hits']
df['PhiX_Percent_Discard'] = (df['PhiX_Discard']/df['PhiX_IN']) * 100

# rRNA SCREEN
df['rRNA_IN'] = df['htsSeqScreener2.totalFragmentsInput']
df['rRNA_Identified'] = df['htsSeqScreener2.Paired_end.PE_hits']
df['rRNA_Percent_Identified'] = (df['rRNA_Identified']/df['rRNA_IN']) * 100

# SUPER DEDUPER
df['SuperDeduper_IN'] = df['htsSuperDeduper.totalFragmentsInput']
df['SuperDeduper_Ignored'] = df['htsSuperDeduper.ignored']
df['SuperDeduper_Duplicate'] = df['htsSuperDeduper.duplicate']
df['SuperDeduper_Percent_Duplicate'] = (df['SuperDeduper_Duplicate']/df['SuperDeduper_IN']) * 100

# ADAPTER TRIM
df['AdapterTrimmed_IN'] = df['htsAdapterTrimmer.totalFragmentsInput']
df['AdapterTrimmed_Reads'] = df['htsAdapterTrimmer.Paired_end.PE_adapterTrim']
df['AdapterTrimmed_Percent_Reads'] = (df['AdapterTrimmed_Reads']/df['AdapterTrimmed_IN']) * 100
df['AdapterTrimmed_BP'] = df['htsAdapterTrimmer.Paired_end.PE_adapterBpTrim']

# QWINDOW TRIM
df['QwindowTrimmed_IN'] = df['htsQWindowTrim.totalFragmentsInput']
df['QwindowTrimmed_R1_LeftBpTrim'] = df['htsQWindowTrim.Paired_end.R1_leftTrim']
df['QwindowTrimmed_R2_LeftBpTrim'] = df['htsQWindowTrim.Paired_end.R2_leftTrim']
df['QwindowTrimmed_R1_RightBpTrim'] = df['htsQWindowTrim.Paired_end.R1_rightTrim']
df['QwindowTrimmed_R2_RightBpTrim'] = df['htsQWindowTrim.Paired_end.R2_rightTrim']
df['QwindowTrimmed_Discard'] = df['htsQWindowTrim.Paired_end.PE_discarded']

# N TRIM
df['NcharTrimmed_IN'] = df['htsNTrimmer.totalFragmentsInput']
df['NcharTrimmed_R1_LeftBpTrim'] = df['htsNTrimmer.Paired_end.R1_leftTrim']
df['NcharTrimmed_R1_RightBpTrim'] = df['htsNTrimmer.Paired_end.R1_rightTrim']
df['NcharTrimmed_R2_LeftBpTrim'] = df['htsNTrimmer.Paired_end.R2_leftTrim']
df['NcharTrimmed_R2_RightBpTrim'] = df['htsNTrimmer.Paired_end.R2_rightTrim']
df['NcharTrimmed_Discard'] = df['htsNTrimmer.Paired_end.PE_discarded']

# CUT TRIM
df['MinLen_IN'] = df["htsCutTrim.totalFragmentsInput"]
df['MinLen_Discard'] = df["htsCutTrim.Paired_end.PE_discarded"]

# POST STATS
df['Proc_Reads'] = df['htsStats2.totalFragmentsInput']
df['Proc_Bp'] = df['htsStats2.Paired_end.R1_bpLen'] + df['htsStats2.Paired_end.R2_bpLen']
df['Proc_R1_PercentQ30'] = (df['htsStats2.Paired_end.R1_bQ30']/df['htsStats2.Paired_end.R1_bpLen']) * 100
df['Proc_R2_PercentQ30'] = (df['htsStats2.Paired_end.R2_bQ30']/df['htsStats2.Paired_end.R2_bpLen']) * 100
df['Proc_Percent_CG'] = ((df['htsStats2.Base_composition.C']+df['htsStats2.Base_composition.G'])/df['Proc_Bp']) * 100
df['Proc_Chars_N'] = df['htsStats2.Base_composition.N']

# FINAL PERCENTS
df['Final_Percent_Read'] = (df['Proc_Reads']/df['Raw_Reads']) * 100
df['Final_Percent_Bp'] = (df['Proc_Bp']/df['Raw_Bp']) * 100


df.to_csv("testing.csv", columns=keys)