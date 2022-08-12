if pre_qc == "Y":
  # Input: merged fastq files
  # Output: pre-trimming fastqc report
  # Description: pre-trimming fastqc
  rule Pretrim_Qc:
      input:
            fastq + "/" + fastq_name + "_R1.fastq.gz",
            fastq + "/" + fastq_name + "_R2.fastq.gz"
      output:
            out + "/pretrim_qc/" + fastq_name + "_R1_fastqc.zip",
            out + "/pretrim_qc/" + fastq_name + "_R2_fastqc.zip"
      threads: 8
      log: "log/pretrim_qc/{sample}_pretrimqc.err"
      shell:
            """
            fastqc {input} \
                   -o {out}/pretrim_qc \
                   -f fastq \
                   --noextract \
                   -t {threads} \
                   2>{log}
            """

# Input: merged fastq files
# Output: trimmed fastq files
# Description: trimming
rule Trim:
    input:
          fastq + "/" + fastq_name + "_R1.fastq.gz",
          fastq + "/" + fastq_name + "_R2.fastq.gz"
    output:
          r1 = out + "/trimmed/{sample}/{sample}_trimmed_R1.fastq.gz",
          r2 = out + "/trimmed/{sample}/{sample}_trimmed_R2.fastq.gz",
          json = out + "/trimmed/{sample}/{sample}_trimming.fastp.json",
          html = out + "/trimmed/{sample}/{sample}_trimming.fastp.html"
    params:
          trim = trim
    threads: 8
    log: "log/trim/{sample}_trim.err"
    shell:
          """
          fastp -i {input[0]} -I {input[1]} \
                -o {output.r1} -O {output.r2} \
                -w {threads} \
                -j {output.json} \
                -h {output.html} \
                {params.trim} \
                2>{log}
          """

# Input: trimmed fastq files
# Output: fastqc report
# Description: post-trimming fastqc
rule Posttrim_Qc:
    input:
          rules.Trim.output.r1,
          rules.Trim.output.r2
    output:
          out + "/posttrim_qc/{sample}_trimmed_R1_fastqc.zip",
          out + "/posttrim_qc/{sample}_trimmed_R2_fastqc.zip",
    threads: 8
    log: "log/posttrim_qc/{sample}_posttrim_qc.err"
    shell:
          """
          fastqc {input} \
                 -o {out}/posttrim_qc \
                 -f fastq \
                 --noextract \
                 -t {threads} \
                 2>{log}
          """


# Input: trimmed fastq files, reference indice, gtf file
# Output: star alignment output
# Description: star alignment
rule Align:
    input:
          r1,
          r2
    output:
          bam = out + "/align/{sample}/{sample}_Aligned.sortedByCoord.out.bam",
          counts = out + "/align/{sample}/{sample}_ReadsPerGene.out.tab",
          txt = out + "/align/{sample}/complete.txt"
    threads: 16
    params:
          indice = star_indice,
          gtf = gtf,
          star = star
    log:
          log1 = "log/align/{sample}_align.err",
          log2 = "log/align/{sample}_index.err"
    shell:
          """
          STAR --runThreadN {threads} \
               --genomeDir {params.indice} \
               --sjdbGTFfile {params.gtf} \
               --readFilesIn {input[0]} {input[1]} \
               --readFilesCommand zcat \
               --outFileNamePrefix {out}/align/{wildcards.sample}/{wildcards.sample}_ \
               --outSAMtype BAM SortedByCoordinate \
               --quantMode GeneCounts \
               --twopassMode Basic \
               {params.star} \
               2>{log.log1}
          samtools index {output.bam} 2>{log.log2}
          touch {output.txt}
          """

# Input: gene counts table per sample
# Output: gene counts table for all samples
# Description: merge gene counts table
rule Counts_Table:
    input:
          expand(rules.Align.output.counts,sample = samples)
    output:
          out + "/quantification/" + run_id + "_" + build + "_counts_by_gene.csv",
          gene_count
    params:
          column = count_column
    threads: 2
    log: "log/quantification/counts_by_gene.err"
    run:
          shell("mkdir -p {out}/quantification")
          df = None
          for file in input:
              id = file.split("/")[-1].split("_ReadsPerGene.out.tab")[0]
              if df is None:
                 df = pd.read_table(file, sep="\t", usecols=[0,count_column], names=["Gene_ID", id])
              else:
                 df_new = pd.read_table(file, sep="\t", usecols=[0,count_column], names=["Gene_ID", id])
                 df = df.merge(df_new, on="Gene_ID")
          df.to_csv(str(output[0]), index=False)
          shell("cp {output[0]} {output[1]}")

# Input: star output bam file
# Output: QC metrics
# Description: Alingment QC
rule Qualimap_Bamqc:
    input:
          rules.Align.output.txt,
          rules.Align.output.bam
    output:
          out + "/qualimap_bamqc/{sample}/complete.txt"
    params:
          strand = qualimap_strand,
          bamqc = bamqc
    threads: 16
    log: "log/qualimap/{sample}_bamqc.err"
    shell:
          """
          mkdir -p ./tmp
          export JAVA_OPTS="-Djava.io.tmpdir=./tmp"
          qualimap bamqc -bam {input[1]} \
                         --java-mem-size={java_mem} \
                         -c \
                         -nt {threads} \
                         -outdir {out}/qualimap_bamqc/{wildcards.sample} \
                         -p {params.strand} \
                         {params.bamqc} \
                         2>{log}
          touch {output}
          """

# Input: star output bam file, gtf file
# Output: QC metrics
# Description: RNAseq QC
rule Qualimap_Rnaseq:
    input:
          rules.Align.output.txt,
          rules.Align.output.bam
    output:
          out + "/qualimap_rnaseq/{sample}/complete.txt"
    params:
          strand = qualimap_strand,
          rnaseq = rnaseq,
          gtf = gtf
    threads: 16
    log: "log/qualimap/{sample}_rnaseq.err"
    shell:
          """
          mkdir -p ./tmp
          export JAVA_OPTS="-Djava.io.tmpdir=./tmp"
          qualimap rnaseq -bam {input[1]} \
                         --java-mem-size=20G \
                         -gtf  {params.gtf} \
                         -outdir {out}/qualimap_rnaseq/{wildcards.sample} \
                         -p {params.strand} \
                         -pe \
                         {params.rnaseq} \
                         2>{log}
          touch {output}
          """

# Input: star output bam file, gencode collapsed gtf file
# Output: QC metrics
# Description: RNAseq QC
rule Rnaseqc:
    input:
          rules.Align.output.txt,
          rules.Align.output.bam
    output:
          out + "/rnaseqc/{sample}/complete.txt"
    params:
          strand = rnaseqc_strand,
          rnaseqc = rnaseqc,
          gtf = gencode_gtf
    threads: 16
    log: "log/rnaseqc/{sample}_rnaseqc.err"
    shell:
          """
          rnaseqc {params.gtf} \
                  {input[1]} \
                  {out}/rnaseqc/{wildcards.sample} \
                  -s {wildcards.sample} \
                  {params.strand} \
                  --coverage \
                  {params.rnaseqc} \
                  2>{log}
          touch {output}
          """

################# QC on targeted regions if applied #############################
# Input: star output bam file, target bed file
# Output: targeted bam file
# Description: filtering off-target reads
rule Target:
    input:
          rules.Align.output.txt,
          rules.Align.output.bam
    output:
          bam = out + "/target/bam/{sample}_target.bam",
          bai = out + "/target/bam/{sample}_target.bam.bai"
    params:
          target = target
    threads: 8
    log:
          log1 = "log/target/{sample}_target.err",
          log2 = "log/target/{sample}_index.err"
    shell:
          """
          samtools view -b -h -L {params.target} {input[1]} | \
          samtools sort -O bam -o {output.bam} 2>{log.log1}
          samtools index {output.bam} 2>{log.log2}
          """

# Input: targeted bam file, target bed file
# Output: QC metrics
# Description: Alingment QC on targeted regions
rule Qualimap_Bamqc_Target:
    input:
          rules.Target.output.bam
    output:
          out + "/qualimap_bamqc/{sample}_target/complete.txt"
    params:
          strand = qualimap_strand,
          bamqc = bamqc,
          bamqc_target = bamqc_target
    threads: 16
    log: "log/qualimap/{sample}_target_bamqc.err"
    shell:
          """
          mkdir -p ./tmp
          export JAVA_OPTS="-Djava.io.tmpdir=./tmp"
          qualimap bamqc -bam {input} \
                         --java-mem-size={java_mem} \
                         -c \
                         -nt {threads} \
                         -outdir {out}/qualimap_bamqc/{wildcards.sample}_target \
                         -p {params.strand} \
                         {params.bamqc} \
                         {params.bamqc_target} \
                         2>{log}
          touch {output}
          """

# Input: targeted bam file, gtf file
# Output: QC metrics
# Description: RNAseq QC on targeted regions
rule Qualimap_Rnaseq_Target:
    input:
          rules.Target.output.bam
    output:
          out + "/qualimap_rnaseq/{sample}_target/complete.txt"
    params:
          strand = qualimap_strand,
          rnaseq = rnaseq,
          gtf = gtf
    threads: 16
    log: "log/qualimap/{sample}_target_rnaseq.err"
    shell:
          """
          mkdir -p ./tmp
          export JAVA_OPTS="-Djava.io.tmpdir=./tmp"
          qualimap rnaseq -bam {input} \
                         --java-mem-size=20G \
                         -gtf  {params.gtf} \
                         -outdir {out}/qualimap_rnaseq/{wildcards.sample}_target \
                         -p {params.strand} \
                         -pe \
                         {params.rnaseq} \
                         2>{log}
          touch {output}
          """

# Input: targeted bam file, gencode collapsed gtf file
# Output: QC metrics
# Description: RNAseq QC on targeted regions
rule Rnaseqc_Target:
    input:
          rules.Target.output.bam
    output:
          out + "/rnaseqc/{sample}_target/complete.txt"
    params:
          strand = rnaseqc_strand,
          rnaseqc = rnaseqc,
          gtf = gencode_gtf
    threads: 16
    log: "log/rnaseqc/{sample}_target_rnaseqc.err"
    shell:
          """
          rnaseqc {params.gtf} \
                  {input} \
                  {out}/rnaseqc/{wildcards.sample}_target \
                  -s {wildcards.sample}_target \
                  {params.strand} \
                  --coverage \
                  {params.rnaseqc} \
                  2>{log}
          touch {output}
          """

# Input: targeted bam file, target bed file
# Output: enrichment metircs
# Description: plotEnrichment from deeptools
rule Plot_Enrichment:
    input:
          expand(rules.Align.output.bam,sample = samples)
    output:
          out + "/target/plot_enrichment/complete.txt"
    params:
          target = target
    threads: 8
    log: "log/target/plotenrichment.err"
    shell:
          """
          plotEnrichment -b {input} \
                         --BED {params.target} \
                         --regionLabels "on target" \
                         --smartLabels \
                         -T "On Target Rate" \
                         --outRawCounts {out}/target/plot_enrichment/{run_id}_{build}_enrichment.txt \
                         -p {threads} \
                         -o {out}/target/plot_enrichment/{run_id}_{build}_enrichment.png \
                         2>{log}
          touch {output}
          """

# Input: enrichment metrics table from output of rule Plot_Enrichment
# Output: enrichment metircs table for multiqc report
# Description: data processing
rule Enrichment_File:
    input:
          rules.Plot_Enrichment.output
    output:
          out + "/target/plot_enrichment/" + run_id + "_" + build + "_onTarget.csv"
    threads: 2
    run:
          file = str(out + "/target/plot_enrichment/" + run_id + "_" + build + "_enrichment.txt")
          df = pd.read_table(file, sep="\t", usecols=[0,2,3], skiprows=1, names=["Sample", "onTarget_Rate_percent", "onTarget_ReadCount"])
          df["Sample"] = df["Sample"].str.replace("_Aligned.sortedByCoord.out","")
          df.to_csv(str(output[0]), index=False)

#################################################################################

# Input: all QC output
# Output: multiqc report
# Description: multiqc
rule Multiqc:
    input:
          lambda wildcards: expand(rules.Pretrim_Qc.output, sample = samples) if pre_qc == "Y" else [],
          lambda wildcards: expand(rules.Posttrim_Qc.output, sample = samples) if trimming == "Y" else [],
          expand(rules.Qualimap_Bamqc.output, sample = samples),
          expand(rules.Qualimap_Rnaseq.output, sample = samples),
          expand(rules.Rnaseqc.output, sample = samples),
          lambda wildcards: expand(rules.Salmon.output, sample = samples) if run_transcript == "Y" else [],
          lambda wildcards: expand(rules.Gffcompare.output, sample = samples) if run_splicing == "Y" else [],
          lambda wildcards: expand(rules.Qualimap_Bamqc_Target.output, sample = samples) if target != "" else [],
          lambda wildcards: expand(rules.Qualimap_Rnaseq_Target.output, sample = samples) if target != "" else [],
          lambda wildcards: expand(rules.Rnaseqc_Target.output, sample = samples) if target != "" else [],
          rules.Enrichment_File.output if target != "" else []
    output:
          out + "/multiqc/" + run_id + "_" + build +"_multiqc_report.html",
          multiqc
    params:
          multiqc_config = multiqc_config
    threads: 8
    log: "log/multiqc/multiqc.err"
    shell:
          """
          set +e
          touch {out}/multiqc/dummy_mqc.csv
          multiqc {out} \
                  --title {run_id}_{build} \
                  -o {out}/multiqc \
                  --ignore *STARpass1* --ignore *fi_workdir* \
                  -f \
                  -c {params.multiqc_config} \
                  2>{log}
          cp {output[0]} {output[1]}
          cp config/config.yaml {out}/archive/
          """
