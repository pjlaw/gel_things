# -*- coding: utf-8 -*-
"""
Created on Wed Oct 11 10:03:50 2023

@author: plaw
"""
import pandas as pd
from scipy import stats

results_path="/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/perturbSeq/annotated_enhancers"
ht29 = pd.read_csv(f"{results_path}/HT29_enhancers_annotated_newTF.txt", sep="\t")
sw480 = pd.read_csv(f"{results_path}/SW480_enhancers_annotated_newTF.txt", sep="\t")
results_file = ht29[["chr","start","end","TF_bound","ABC"]].merge(sw480[["chr","start","end","TF_bound","ABC"]], on=["chr","start","end"], suffixes=('_ht29', '_sw480'), how="outer")

datafiles={"ActiveDriverWGSR":"activedriver.tsv", "MutSpot_indel":"mutspot_indel.tsv", "MutSpot_snv":"mutspot_snv.tsv", "OncodriveCLUSTL":"oncodriveClustl.tsv", "OncodriveFML":"oncodriveFML.tsv"}


def find_overlap(this_dat, this_name, enh_chrom, enh_start, enh_end):
    overlaps = this_dat.query(f"chrom=='{enh_chrom}' and start<={enh_end} and end>={enh_start}")
    if this_name == "OncodriveFML":
        #oncodrivefml
        wanted_cols = ["MUTS", "MUTS_RECURRENCE", "SAMPLES", "P_VALUE", "Q_VALUE"]
    elif this_name == "OncodriveCLUSTL":
        #oncodriveclustl
        wanted_cols = ['TOTAL_MUT', 'CLUSTERED_MUT', 'P_TOPCLUSTER','Q_TOPCLUSTER']
    elif this_name == "ActiveDriverWGSR":
        #activedriver
        wanted_cols = ["pp_element","element_muts_obs","fdr_element"]
    else:
        #mutspot
        wanted_cols = ["pval", "fdr", "k"]

    #add the name
    renamed_cols = [f"{this_name}_{s}" for s in wanted_cols]

    if overlaps.shape[0]==0:
        output={k:[""] for k in renamed_cols}
        return(pd.DataFrame(output))
    else:
        wanted_data = overlaps[wanted_cols]
        wanted_data.columns = renamed_cols
        if wanted_data.shape[0]>1:
            collapsed_row = wanted_data.apply(lambda row: ";".join(row.astype(str)), axis=0)

            # Create a new DataFrame with the concatenated rows
            collapsed_df = pd.DataFrame(collapsed_row).transpose()
            return(collapsed_df)
        else:
            return(wanted_data)

all_results=[results_file]
datapath="/data/scratch/DGE/DUDGE/MOPOPGEN/plaw/perturbSeq/CRC_noncoding_for_James/CRCpeterbseq"
for toolname, filename in datafiles.items():
    print(toolname)
    this_data = pd.read_csv(f"{datapath}/{filename}", sep="\t", index_col=0)
    if not toolname.startswith("MutSpot"):
        if toolname.startswith("Onco"):
            colname="SYMBOL"
        else:
            colname="id"
        this_data[["chrom", "start", "end"]]=this_data[colname].str.split("_", expand=True) #r"[:|-]"
        this_data[["start", "end"]] = this_data[["start", "end"]].apply(pd.to_numeric)
        #this_data["enhancer"]=this_data.apply(lambda row: "_".join(row[["chr", "start", "end"]]), axis=1)
    tool_results = results_file.apply(lambda row: find_overlap(this_data, toolname, row["chr"], row["start"], row["end"]), axis=1)
    # Concatenate the resulting DataFrames
    tool_results = pd.concat([tool_results[i] for i in tool_results.index], axis=0)
    tool_results.reset_index(drop=True, inplace=True)
    all_results.append(tool_results)

output_file = pd.concat(all_results, axis=1)
output_file.to_csv(f"{results_path}/gel_annotated_enhancers.txt", sep="\t", index=False)





