# variant_calling.py

import os

def run_freebayes_variant_calling(bam_file, reference_genome, output_vcf):
    """
    Call variants using FreeBayes.

    Parameters:
    - bam_file: Path to input sorted BAM file
    - reference_genome: Path to reference genome FASTA file
    - output_vcf: Path to output VCF file
    """
    print(f"🔬 Running FreeBayes on {bam_file}")
    command = f"freebayes -f {reference_genome} {bam_file} > {output_vcf}"
    os.system(command)
    print(f"✅ Variant calling complete. Output: {output_vcf}")
