# preprocess.py

import os

def run_fastqc(input_fastq, output_dir):
    os.system(f"fastqc {input_fastq} -o {output_dir}")

def run_trimming(input_fastq, output_fastq):
    # Example using fastp for trimming
    os.system(f"fastp -i {input_fastq} -o {output_fastq} --detect_adapter_for_pe -w 4")

def run_alignment(fastq, ref_genome, output_bam):
    # Example BWA-MEM2 alignment pipeline
    sam_file = output_bam.replace(".bam", ".sam")
    os.system(f"bwa mem {ref_genome} {fastq} > {sam_file}")
    os.system(f"samtools view -bS {sam_file} > {output_bam}")
    os.system(f"samtools sort {output_bam} -o {output_bam}")
    os.system(f"samtools index {output_bam}")

def run_variant_calling(bam_file, ref_genome, output_vcf):
    # Example FreeBayes
    os.system(f"freebayes -f {ref_genome} {bam_file} > {output_vcf}")
