#!/bin/bash -ue
fastp -i data/raw/Cauris.fastq -o data/processed/Cauris_trimmed.fastq --detect_adapter_for_pe -w 4
