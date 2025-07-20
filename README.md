# 🧠 PathogenHawk

**PathogenHawk** is a modular machine learning toolkit for predicting antimicrobial resistance (AMR) from genome sequences of fungal and bacterial pathogens.

It supports multiple species including *Candida auris*, *E. coli*, and *Aspergillus fumigatus*, using genomic features (e.g. SNPs, resistance genes) to train interpretable ML models.

---

## 📄 Project Links
- 📂 [Source Code](https://github.com/biosciences/PathogenHawk): Explore the full repository
- 🔗 [Live Report](https://biosciences.github.io/PathogenHawk/demo_cauris.html): View the interactive HTML output

## 📚 Rationale

C. auris is an emerging pathogen with significant resistance concerns Novel antifungals and treatment approaches to tackle resistance and improve outcomes of invasive fungal disease.

---

## 📁 Project Structure

```
PathogenHawk/
├── configs/                 # YAML configs for each pathogen
├── data/                    # Raw and processed data
├── genome_processing/       # Preprocessing and variant calling
├── feature_engineering/     # Genomic feature extraction
├── ml_model/                # Model training, evaluation, interpretation
├── genes/                   # Resistance gene annotations
├── metadata/                # MIC values and resistance phenotype labels
├── ref/                     # Reference genomes (e.g. Cauris_CDC317.fa)
├── notebooks/               # Jupyter demo notebooks
└── scripts/                 # Utility scripts
```

---

## 🚀 Quick Start

1. Clone the repository and install dependencies:

```bash
conda env create -f environment.yml
conda activate pathogenhawk
```

2. Create a configuration file under `configs/`, e.g. `configs/Cauris.yaml`:
```yaml
pathogen: candida_auris
reference_genome: data/raw/C_auris_B8441V3_ref.fa
resistance_genes: data/metadata/Cauris_amr_genes.tsv
phenotype_file: data/metadata/Cauris_MIC.csv
feature_file: data/metadata/Cauris_features.tsv
labels: fluconazole_resistance
alignment_tool: bwa
variant_caller: freebayes
features: ["snp", "resistance_genes"]
ml_model: xgboost
output_dir: results/Candida_auris
```

3. Run your pipeline:
```bash
export PYTHONPATH=$(pwd)
python ml_model/train_model.py --config configs/Cauris.yaml
```

#### Run Nextflow

```bash
# For Candida auris
nextflow run pathogenhawk.nf --config configs/Cauris.yaml

# For Escherichia coli
nextflow run pathogenhawk.nf --config configs/Ecoli.yaml
```

#### Configuration
`workflow/nextflow.config` includes basic container settings. You can customize profiles for HPC/cloud usage.

---

## 📦 Example Files

### For Candida auris

- [`data/metadata/auris_amr_genes.tsv`](data/metadata/Cauris_amr_genes.tsv): Known resistance genes
- [`data/metadata/cauris_MIC.csv`](data/metadata/Cauris_MIC.csv): MIC values and resistance phenotypes

### For Escherichia coli

- [`data/metadata/Ecoli_res_genes.tsv`](data/metadata/Ecoli_res_genes.tsv): Known resistance genes
- [`data/metadata/Ecoli_MICs.csv`](data/metadata/Ecoli_MICs.csv): MIC values and resistance phenotypes

---

## 📚 Citation

If you use **PathogenHawk** in your research, please cite the associated JOSS paper (under review):

> Lai, K. (2025). *PathogenHawk: A Pathogen Machine Learning Toolkit for Predicting Antimicrobial Resistance from Genomic Features*. Journal of Open Source Software (under review). https://github.com/biosciences/PathogenHawk


---

## 👩‍💻 Author

Developed by Kaitao Lai

## 🪪 License

MIT © 2025 Kaitao Lai
