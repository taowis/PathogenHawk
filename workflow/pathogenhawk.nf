#!/usr/bin/env nextflow

nextflow.enable.dsl=2

params.config = "configs/Cauris.yaml"

workflow {
    config_ch = Channel.fromPath(params.config)
    config_info = load_config(config_ch)

    preprocess_reads(config_info)
    call_variants(config_info)
    extract_features(config_info)
}

process load_config {
    input:
    path config_file

    output:
    tuple val("${config_file.getBaseName()}"), path(config_file)

    script:
    """
    echo "Loaded config: ${config_file}"
    """
}

process preprocess_reads {
    input:
    tuple val(name), path(config_file)

    output:
    path "data/processed/${name}_trimmed.fastq"

    script:
    """
    fastp -i data/raw/${name}.fastq -o data/processed/${name}_trimmed.fastq --detect_adapter_for_pe -w 4
    """
}

process call_variants {
    input:
    tuple val(name), path(config_file)

    output:
    path "data/processed/${name}.vcf"

    script:
    """
    bwa index ${name}.fa
    bwa mem ${name}.fa data/processed/${name}_trimmed.fastq > ${name}.sam
    samtools view -bS ${name}.sam > ${name}.bam
    samtools sort ${name}.bam -o ${name}.sorted.bam
    samtools index ${name}.sorted.bam
    freebayes -f ${name}.fa ${name}.sorted.bam > data/processed/${name}.vcf
    """
}

process extract_features {
    input:
    tuple val(name), path(config_file)

    output:
    path "data/processed/${name}_features.tsv"

    script:
    """
    python feature_engineering/build_features.py --config ${config_file} --output data/processed/${name}_features.tsv
    """
}
