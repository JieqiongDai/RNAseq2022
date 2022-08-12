# Input: trimmed fastq files, reference indice
# Output: star fusion output
# Description: fusion calling
rule Star_Fusion:
    input:
          r1,
          r2
    output:
          file = out + "/fusion/{sample}/star-fusion.fusion_predictions.tsv",
          txt = out + "/fusion/{sample}/complete.txt",
          fusions = fusions
    threads: 16
    params:
          star_fusion = star_fusion,
          genome_lib = genome_lib,
          bind = bind,
          image = image
    log: "log/fusion/{sample}_star_fusion.err"
    shell:
          """
          SINGULARITY_BINDPATH={params.bind}
          export SINGULARITY_BINDPATH
          SINGULARITY_CACHEDIR=./ \
          singularity exec -e {params.image} \
          STAR-Fusion --left_fq {input[0]} \
                      --right_fq {input[1]} \
                      --genome_lib_dir {params.genome_lib} \
                      --CPU {threads} \
                      --output_dir {out}/fusion/{wildcards.sample} \
                      --FusionInspector validate \
                      --examine_coding_effect \
                      --denovo_reconstruct
                      {params.star_fusion} \
                      2>{log}
          touch {output.txt}
          cp {output.file} {output[2]}
          """
