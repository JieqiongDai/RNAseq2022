# Input: star output bam file or targeted bam file, gtf file
# Output: stringTie output
# Description: splicing analysis
rule Stringtie:
    input:
          rules.Align.output.txt,
          rules.Align.output.bam,
          rules.Target.output.bam if target != "" else []
    output:
          gtf = out + "/splicing/{sample}/{sample}.gtf",
          splicing_gtf = splicing_gtf
    params:
          stringtie = config["stringtie"],
          gtf = gtf,
          strand = stringtie_strand
    threads: 16
    log: "log/splicing/{sample}_stringtie.err"
    run:
         if target != "":
             bam = {input[2]}
         else:
             bam = {input[1]}
         shell("""
               stringtie {params.strand} \
                         -G {params.gtf} \
                         -p {threads} \
                         -A {out}/splicing/{wildcards.sample}/{wildcards.sample}_gene_abund.tab \
                         {params.stringtie} \
                         -o {output.gtf} \
                         {bam} \
                         2>{log} ;
               cp {output.gtf} {output[1]}
               """)

# Input: stringtie output gtf file, reference sequence, gtf file
# Output: gffcompare output
# Description: gff compare
rule Gffcompare:
    input:
          rules.Stringtie.output.gtf
    output:
          gtf = out +"/splicing/{sample}/{sample}_gffcmp/{sample}_gffcmp.annotated.gtf",
          stat = out +"/splicing/{sample}/{sample}_gffcmp/{sample}_gffcmp.stats",
          splicing_gtf_anno = splicing_gtf_anno,
          splicing_stat = splicing_stat
    params:
          gffcmp = gffcmp,
          gtf = gtf,
          ref = ref
    threads: 2
    shell:
          """
          mkdir -p {out}/splicing/{wildcards.sample}/{wildcards.sample}_gffcmp
          cd {out}/splicing/{wildcards.sample}/{wildcards.sample}_gffcmp
          gffcompare -r {params.gtf} \
                     -R -Q \
                     -s {params.ref} \
                     {params.gffcmp} \
                     {run}{input} \
                     -o {wildcards.sample}_gffcmp \
                     2>{run}log/splicing/{wildcards.sample}_gffcompare.err

          cp {run}{output.gtf} {run}{output[2]}
          cp {run}{output.stat} {run}{output[3]}
          """

# Input: stringtie output gtf file, reference sequence
# Output: fasta file
# Description: convert gtf to fasta sequence
rule Gffread:
    input:
          rules.Stringtie.output.gtf
    output:
          out + "/splicing/{sample}/{sample}_gffread/{sample}_gffread.fa",
          splicing_fasta
    params:
          ref = ref
    threads: 2
    log: "log/splicing/{sample}_gffread.err"
    shell:
          """
          gffread {input} \
                  -g {params.ref} \
                  -w {output[0]} \
                  2>{log}
          cp {output[0]} {output[1]}
          """
