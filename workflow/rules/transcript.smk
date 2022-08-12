# Input: trimmed fastq files, reference indice
# Output: salmon alignment output
# Description: Traanscript level alignment
rule Salmon:
    input:
          r1,
          r2
    output:
          counts = out + "/transcript/{sample}/quant.sf",
          txt = out + "/transcript/{sample}/complete.txt"
    threads: 16
    params:
          indice = salmon_indice,
          salmon = salmon
    log: "log/transcript/{sample}_salmon.err"
    shell:
          """
          salmon quant -i {params.indice} \
                       -p {threads} \
                       -l A \
                       -1 {input[0]} \
                       -2 {input[1]} \
                       --validateMappings \
                       -o {out}/transcript/{wildcards.sample} \
                       {params.salmon} \
                       2>{log}
          touch {output.txt}
          """

# Input: transcript counts table per sample
# Output: transcipt counts table for all samples
# Description: merge transcipt counts tables
rule Counts_Table_Transcirpt:
    input:
          expand(rules.Salmon.output.counts,sample = samples)
    output:
          out + "/quantification/" + run_id + "_" + build + "_counts_by_transcript.csv",
          out + "/quantification/" + run_id + "_" + build + "_TPM_by_transcript.csv",
          out + "/quantification/" + run_id + "_" + build + "_transcripts_detected.csv",
          transcript_count,
          tpm
    threads: 2
    run:
          shell("mkdir -p {out}/quantification")
          df1 = None
          df2 = None
          df3 = None
          for file in input:
              id = file.split("transcript/")[-1].split("/quant.sf")[0]
              if ((df1 is None) and (df2 is None) and (df3 is None)):
                 df1 = pd.read_table(file, sep="\t", usecols=[0,4], skiprows=1, names=["Transcript_ID", id])
                 count = str(len(df1[df1[id] > 4]))
                 df2 = pd.read_table(file, sep="\t", usecols=[0,3], skiprows=1, names=["Transcript_ID", id])
                 df3 = pd.DataFrame({"Sample":[id],"#_Transcripts":[count]})
              elif ((df1 is not None) and (df2 is not None) and (df3 is not None)):
                 df_new1 = pd.read_table(file, sep="\t", usecols=[0,4], skiprows=1, names=["Transcript_ID", id])
                 df1 = pd.merge(df1, df_new1)
                 count = str(len(df_new1[df_new1[id] > 4]))
                 df_new2 = pd.read_table(file, sep="\t", usecols=[0,3], skiprows=1, names=["Transcript_ID", id])
                 df2 = pd.merge(df2, df_new2)
                 df_new3 = pd.DataFrame({"Sample":[id],"#_Transcripts":[count]})
                 df3 = pd.concat([df3,df_new3])
          df1.to_csv(str(output[0]), index=False)
          df2.to_csv(str(output[1]), index=False)
          df3.to_csv(str(output[2]), index=False)
          shell("cp {output[0]} {output[3]}; cp {output[1]} {output[4]}")
