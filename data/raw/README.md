# Data source

## For the Candida auris reference genome FASTA 
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/002/759/435/GCA_002759435.3_Cand_auris_B8441_V3/GCA_002759435.3_Cand_auris_B8441_V3_genomic.fna.gz
gunzip GCA_002759435.3_Cand_auris_B8441_V3_genomic.fna.gz
mv GCA_002759435.3_Cand_auris_B8441_V3_genomic.fna C_auris_B8441V3_ref.fa
```

## For the Escherichia coli reference genome FASTA (E_coli_ref.fa)
```bash
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna E_coli_ASM584V2_ref.fa
```