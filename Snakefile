import os
import re

#Location of raw files
raw_seq_folder = "data/bam"

#Make sample list
file_list = [f for f in sorted(os.listdir(raw_seq_folder)) if (str(f))[-4:] == ".bam"]
samples = []
for file in file_list:
    dir_and_file = raw_seq_folder + "/" + file
    if re.sub(".bam", "",file) not in samples:
      samples.append(re.sub(".bam", "",file))

localrules: all, clean

#Make all results listed as inputs to this rule
rule all:
    input:
        expand("data/fastx/{strain}.fasta.gz", strain=samples),
        expand("data/assembly/{strain}_flye/{strain}_polished.fasta", strain=samples)

rule bam2fastx: #Convert PacBio BAM files to fasta for assembler
    input:
        bam="data/bam/{strain}.bam"
    output:
        fasta = "data/fastx/{strain}.fasta.gz"
    conda:
          "code/bam2fastx.yaml"
    shell:
        """
        #-c for bam2fastx sets the gzip compression level, lower = faster but less compressed
        bam2fasta -c 9 -o data/fastx/{wildcards.strain} {input.bam}
        """

rule flye_assemble: #assemble genomes using Flye assembler, have found to be more succesfull than HGAP3/FALCON
    input: "data/fastx/{strain}.fasta.gz"
    output: "data/assembly/{strain}_flye/assembly.fasta"
    conda:
        "code/flye.yaml"
    resources: cpus=10, mem_mb=100000, time_min=120
    shell:
        """
        mkdir -p data/assembly/{wildcards.strain}_flye
        flye --pacbio-raw {input} --out-dir data/assembly/{wildcards.strain}_flye --threads {resources.cpus}
        """

rule fix_start: #Fix the start point of the assembly to be DNAa using circlator's fixstart tool, which doesn't require prior annotation
    input: "data/assembly/{strain}_flye/assembly.fasta"
    output: "data/assembly/{strain}_flye/{strain}_fixstart.fasta"
    conda:
        "code/circlator.yaml"
    shell:
        """
        circlator fixstart {input} data/assembly/{wildcards.strain}_flye/{wildcards.strain}_fixstart
        """

rule align: #align raw reads to assembly for polishing using PacBio's minimap2 aligner
    input:
        assembly = "data/assembly/{strain}_flye/{strain}_fixstart.fasta",
        bam = "data/bam/{strain}.bam"
    output:
        bam_out = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam",
        bai = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam.bai"
    conda:
        "code/pbmm2.yaml"
    resources: cpus=20, mem_mb=800000, time_min=120
    shell:
        """
        pbmm2 index {input.assembly} {input.assembly}.mmi
        pbmm2 align {input.assembly} {input.bam} {output.bam_out} --sort --sort-memory 50G -j {resources.cpus}
        """

rule polish: #polish the genome using arrow algorithm in pacbio gcpp
    input:
        alignment = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam",
        assembly = "data/assembly/{strain}_flye/{strain}_fixstart.fasta",
        bai = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam.bai"
    output:
        fasta = "data/assembly/{strain}_flye/{strain}_polished.fasta"
    conda:
        "code/pbgcpp.yaml"
    resources: cpus=64, mem_mb=800000, time_min=1000
    shell:
        """
        # --log-level TRACE or --log-level DEBUG can also be used
        gcpp --log-level INFO -j {resources.cpus} -r {input.assembly} -o {output.fasta} {input.alignment}
        """

#Dangerous!
#Cleanup all results
rule clean:
    shell:
        """
        rm -rf data/assembly/*
        rm -rf data/fastx/*
        """

#BLASR can also be used for alignment, although it is slower and less accurate
#Was having issues polishing and thought it might be due to incompatibilities with mm2 and pbgcpp
#But, they seem to be working reliably once a previous version (1.0) of gcpp is used

# rule align_blasr:
#     input:
#         assembly = "data/assembly/{strain}_flye/{strain}_fixstart.fasta",
#         bam = "data/bam/{strain}.bam"
#     output:
#         bam_out = "data/assembly/{strain}_flye/{strain}_aligned.bam"
#     conda:
#         "code/blasr.yaml"
#     resources: cpus=16, mem_mb=500000, time_min=1440
#     shell:
#         """
#         blasr {input.bam} {input.assembly} --bam --placeGapConsistently --nproc {resources.cpus} --out {output.bam_out}
#         """

##If blasr is used, the bam file has to then be sorted and indexed
# rule make_bai:
#     input: "data/assembly/{strain}_flye/{strain}_aligned.bam"
#     output:
#         bai = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam.bai",
#         bamSorted = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam"
#     conda:
#         "code/samtools.yaml"
#     resources: cpus=8, mem_mb=200000, time_min=1440
#     shell:
#         """
#         samtools sort -@ {resources.cpus} -m 16G -o {output.bamSorted} -l 9 {input}
#         samtools index {output.bamSorted}
#         """

# rule make_pbi: #make pacbio index for genomic_consensus
#     input: "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam"
#     output:
#         pbi = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam.pbi"
#     conda:
#         "code/pbbam.yaml"
#     threads: 4
#     shell:
#         """
#         pbindex {input}
#         """

##Was having issues with gcpp, so tried arrow in the genomicconsensus conda package
##Issues have been resolved, so using gcpp instead for polishing
# rule arrow_polish:
#     input:
#         alignment = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam",
#         assembly = "data/assembly/{strain}_flye/{strain}_fixstart.fasta",
#         bai = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam.bai",
#         pbi = "data/assembly/{strain}_flye/{strain}_aligned_sorted.bam.pbi"
#     output:
#         fasta = "data/assembly/{strain}_flye/{strain}_arrow_polished.fasta"
#     conda:
#         "code/genomic_consensus.yaml"
#     resources: cpus=8, mem_mb=800000, time_min=360
#     shell:
#         """
#         arrow -j {resources.cpus} -r {input.assembly} -o {output.fasta} {input.alignment}
#         """
