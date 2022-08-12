# Short-Reads RNA-Sequencing data QC and analysis

##### Author: Jieqiong Dai

## Description
This snakemake pipeline is for short-reds RNAseq data QC and analysis. The pipeline may be run on an HPC or in a local environment.

Major functions in this pipeline include:
* Standard data processing and QC, including trimming (optional), gene-level alignment and quantification, and QC
* Optional transcript-level alignment, quantification and QC
* Optional fusion calling
* Optional splicing analysis

The final deliverables include:
* A multiQC report of all QC metrics with all samples (standard)
* A table of read counts at gene level with all samples (standard)
* A table of read counts at transcript level with all samples, and a table of transcript TPM with all samples (optional)
* Files of fusion calling results for each sample (optional)
* Files of splicing analysis results for each sample (optional)

## Software requirement
* [Conda](https://conda.io/projects/conda/en/latest/user-guide/install/index.html)
* [Snakemake](https://snakemake.readthedocs.io/en/stable/)
* [Fastp](https://github.com/OpenGene/fastp)
* [FastqQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [STAR](https://github.com/alexdobin/STAR)
* [Samtools](http://www.htslib.org/)
* [QualiMap](http://qualimap.conesalab.org/)
* [RNASeQC](https://github.com/getzlab/rnaseqc)
* [DeepTools](https://deeptools.readthedocs.io/en/develop/index.html)
* [MultiQC](https://multiqc.info/)
* [Salmon](https://salmon.readthedocs.io/en/latest/)
* [StringTie](https://github.com/gpertea/stringtie)
* [GffCompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
* [GffRead](https://github.com/gpertea/gffread)
* [STAR-Fusion](https://github.com/STAR-Fusion/STAR-Fusion/wiki)
* [Singularity](https://sylabs.io/guides/3.0/user-guide/quick_start.html)

## User's guide
### I. Input requirement
* Edited config/config.yaml
* Demultiplexed fastq files
* Reference genome file and annotation files
* STAR index
* Salmon index (optional, required to run transcript-level analysis)
* STAR-Fusion genome library (optional, required to run fusion calling)
* Bed file including targeted regions (optional, required to run target enrichment analysis)

### II. Editing the config.yaml
#### Basic (need to be modified case by case)
* `run_id`: Run ID
* `run_transcript`: If run transcript level analysis, Y or N, default N
* `run_fusion`: If run fusion calling, Y or N, default N
* `run_splicing`: If run splicing analysis, Y or N, default N
* `trimming`: If perform trimming or not, Y or N, default Y
* `pre_qc`: If needs pre-trimming QC, Y or N, default N
* `strand`: Library strandedness, non, F, or R. F refers to the 1st read strand aligned with RNA (htseq-count option -s yes); R refers to the 2nd read strand aligned with RNA (htseq-count option -s reverse)
* `fastq`: Path to the merged fastq files stored directory.
* `name_split`: Sample ID delimitation. eg. _R, _sub_R
* `fastq_name`: Naming format of fastq files. eg. {sample}, {sample}_sub
* `out`: Path to the output directory. Default output/ in working directory
* `build`: Reference genome build, eg. hg19, hg38
* `reference`: Path to reference genome fasta file
* `gtf`: Path to reference genome annotation GTF file
* `star_indice`: Path to STAR index
* `gencode_gtf`: Path to Gencode collapsed annotation file, for rRNA estimation. Use [collapse_annotation.py](https://github.com/broadinstitute/gtex-pipeline/blob/master/gene_model/collapse_annotation.py) to generate the file if not available.
* `java_mem`: Java memory, default 20G

#### Optional
* `salmon_indice`: Path to Salmon index, required to run transcript-level analysis
* `target`: Path to target enrichment bed file if applied
* `bind`: Path to Singularity bind-mounted directories, required to run fusion calling
* `image`: Path to STAR-Fusion Singularity image, required to run fusion calling
* `trim`: Optional Fastp trimming parameters
* `star`: Optional STAR alingment parameters
* `bamqc`: Optional Qualimap BamQC Parameters
* `rnaseq`: Optional Qualimap RNAseq parameters
* `rnaseqc`: Optional RNA-SeQC parameters
* `dt_enrich`: Optional deepTools plotEnrichment parameters
* `salmon`: Optional Salmon parameters
* `stringtie`: Optional StringTie parameters
* `gffcmp`: Optional GffCompare parameters
* `star_fusion`: Optional STAR-Fusion parameters
* `genome_lib`: Path to STAR-Fusion genome_lib, required to run fusion calling

### III. To run
#### Clone repository
* Clone the repository to your working directory
```bash
git clone git@github.com:JieqiongDai/RNAseq2022.git
```

#### Environment
* Install Conda if not pre-installed
* If not pre-installed, create a conda environment in the appropriate directory using the provided yaml file from the git cloned directory `workflow/envs/` and activate it after installation:
```bash
conda env create -f RNAseq_env.yaml
conda activate RNAseq
```
* Install Singularity if not pre-installed
* Download STAR-FUSION v10.0.1 full package and generate the Singularity image

#### Running the code
* Edit and save `config/config.yaml`
* To run on an HPC using slurm job scheduler:
  Edit `config/cluster_slurm_config.yaml` according to your HPC information;
  Run `wrapper/run_sbatch.sh` to initiate running of the pipeline
```bash
bash wrapper/run_sbatch.sh
```  
* To run on an HPC using SGE/UGE job scheduler:
  Edit `config/cluster_SGE_config.yaml` according to your HPC information;
  Run `wrapper/run_qsub.sh` to initiate running of the pipeline  
```bash
bash wrapper/run_qsub.sh
```  
* To run on a local server:
  Run `wrapper/snakemake.bash` to initiate running of the pipeline
```bash
bash wrapper/snakemake.batch
```  
* Look in log directory for logs for each rule
* To view the snakemake rule graph:
```bash
snakemake --rulegraph | dot -T png > RNAseq.png
```

### IV. Example output
```bash
.  user/defined/output_dir
├── align
│   ├── {sample_A}
│   │   ├── {sample_A}_Aligned.sortedByCoord.out.bam
│   │   ├── {sample_A}_Aligned.sortedByCoord.out.bam.bai
│   │   ├── {sample_A}_ReadsPerGene.out.tab
│   │   └── ...
│   ├── {sample_B}
│   │   ├── {sample_B}_Aligned.sortedByCoord.out.bam
│   │   ├── {sample_B}_Aligned.sortedByCoord.out.bam.bai
│   │   ├── {sample_B}_ReadsPerGene.out.tab
│   │   └── ...
│   ├── {sample_C}
│   │   ├── {sample_C}_Aligned.sortedByCoord.out.bam
│   │   ├── {sample_C}_Aligned.sortedByCoord.out.bam.bai
│   │   ├── {sample_C}_ReadsPerGene.out.tab
│   │   └── ...
├── archive
│   ├── config.yaml
│   └── deliverbale
│       ├── fusion (optional)
│       │   ├── {build}_{sample_A}_star-fusion.fusion_predictions.tsv
│       │   ├── {build}_{sample_B}_star-fusion.fusion_predictions.tsv
│       │   ├── {build}_{sample_C}_star-fusion.fusion_predictions.tsv
│       ├── splicing (optional)
│       │   ├── {build}_{sample_A}_gffcmp.annotated.gtf
│       │   ├── {build}_{sample_A}_gffcmp.stats
│       │   ├── {build}_{sample_A}_gffread.fa
│       │   ├── {build}_{sample_A}.gtf
│       │   ├── {build}_{sample_B}_gffcmp.annotated.gtf
│       │   ├── {build}_{sample_B}_gffcmp.stats
│       │   ├── {build}_{sample_B}_gffread.fa
│       │   ├── {build}_{sample_B}.gtf
│       │   ├── {build}_{sample_C}_gffcmp.annotated.gtf
│       │   ├── {build}_{sample_C}_gffcmp.stats
│       │   ├── {build}_{sample_C}_gffread.fa
│       │   └── {build}_{sample_C}.gtf
│       ├── {run_id}_{build}_counts_by_gene.csv
│       ├── {run_id}_{build}_counts_by_transcript.csv (optional)
│       ├── {run_id}_{build}_multiqc_report.html
│       └── {run_id}_{build}_TPM_by_transcript.csv (optional)
├── fusion (optional)
├── multiqc
│   ├── {run_id}_{build}_multiqc_report_data
│   └── {run_id}_{build}_multiqc_report.html
├── posttrim_qc (optional)
├── pretrim_qc (optional)
├── qualimap_bamqc
│   ├── {sample_A}
│   ├── {sample_A}_target (optional)
│   ├── {sample_B}
│   ├── {sample_B}_target (optional)
│   ├── {sample_C}
│   └── {sample_C}_target (optional)
├── qualimap_rnaseq
│   ├── {sample_A}
│   ├── {sample_A}_target (optional)
│   ├── {sample_B}
│   ├── {sample_B}_target (optional)
│   ├── {sample_C}
│   └── {sample_C}_target (optional)
├── quantification
│   ├── {run_id}_{build}_counts_by_gene.csv
│   ├── {run_id}_{build}_counts_by_transcript.csv
│   ├── {run_id}_{build}_TPM_by_transcript.csv
│   └── {run_id}_{build}_transcripts_detected.csv
├── rnaseqc
│   ├── {sample_A}
│   ├── {sample_A}_target (optional)
│   ├── {sample_B}
│   ├── {sample_B}_target (optional)
│   ├── {sample_C}
│   └── {sample_C}_target (optional)
├──target (optional)
│   ├── bam
│   └── plot_enrichment
├── splicing (optional)
├── transcript (optional)
└── trimmed (optional)
```
