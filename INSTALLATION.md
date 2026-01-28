# Installation Guide

This document provides detailed installation options for the Virulence Factors Classifier (VF_classifier). For most users, the recommended approach using the environment.yml file is the simplest.

A functioning installation is composed of three main steps:

1. Clone the repository (or download the zip file)
2. Install the dependencies (using mamba or conda)
3. Download the VFDB fasta file and create the BLAST database

## Quick Start (Recommended)

For most users, follow these simple steps:

```bash
# Clone the repository
git clone https://github.com/brunogoncalves/VF_classifier.git
cd VF_classifier

# Make script executable
chmod +x *.py *.R

# Install mamba (if not already installed)
conda install -c conda-forge mamba

# Create environment from provided environment.yml file
mamba env create -f environment.yml

# Activate the environment
conda activate VF_classifier
```

That's it! You're ready to use the tool. The environment.yml file includes all necessary dependencies.

---

## Alternative Installation Methods

### Option 1: Complete Installation with Mamba

Mamba is a faster alternative to conda and can install all dependencies, including system tools:

```bash
# Create a new conda environment with mamba
conda create -n VF_classifier -c conda-forge mamba
conda activate VF_classifier

# Install all dependencies, including BLAST and Krona
mamba install -c bioconda -c conda-forge blast krona pandas biopython matplotlib seaborn
```

### Option 2: Complete Installation with Conda

If you prefer to use conda instead of mamba, be warned that it will be **much slower**:

```bash
# Create a new conda environment
conda create -n VF_classifier -c conda-forge
conda activate VF_classifier

# Install all dependencies, including BLAST and Krona
conda install -c bioconda -c conda-forge blast krona pandas biopython matplotlib seaborn
```

### Option 3: Mixed Installation (System BLAST + Python packages)

If you prefer to install BLAST+ system-wide and use pip for Python packages:

#### Install BLAST+ system-wide

```bash
# macOS (using Homebrew):
brew install blast

# Ubuntu/Debian:
sudo apt-get update
sudo apt-get install ncbi-blast+

# CentOS/RHEL:
sudo yum install ncbi-blast+

# From Source:
# Download from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
```

#### Install Python dependencies

```bash
# Create virtual environment (recommended)
python -m venv vf-env
source vf-env/bin/activate  # On Windows: vf-env\Scripts\activate

# Install required dependencies
pip install pandas biopython matplotlib seaborn
```

#### Install Krona Tools (optional)

##### Option A: Ubuntu/Debian System Package (Simplest)

```bash
# Install Krona via system package manager
sudo apt update && sudo apt install radiant
```

##### Option B: Conda Package

```bash
# Using conda (if not installed with system packages)
conda install -c bioconda krona
```

##### Option C: Manual Installation

```bash
# Manual installation
wget https://github.com/marbl/Krona/releases/download/v2.7/KronaTools-2.7.tar
tar -xvf KronaTools-2.7.tar
cd KronaTools-2.7
./install.pl
```

**Note**: All three options provide the same `ktImportText` command used by the script. The `radiant` package is the Ubuntu/Debian system package for Krona tools.

### Option 4: Minimal Installation (Core functionality only)

For basic functionality without visualizations:

```bash
# Install BLAST+ (see Option 3 for system-specific commands)
# Install only required Python packages
pip install pandas biopython
```

### Option 5: Advanced Heatmaps (R + ComplexHeatmap)

To generate advanced comparative heatmaps using `plot_vf_heatmap.R`, you must have R installed. The script will automatically attempt to install required packages (ComplexHeatmap, dplyr, etc.) on its first run.

#### Install R system-wide (in case you don't have it yet)

```bash
# macOS (using Homebrew):
brew install r

# Ubuntu/Debian:
sudo apt-get update
sudo apt-get install r-base

# CentOS/RHEL:
sudo yum install R
```

#### Install within Conda (Recommended if using Conda)

If you are using the provided environment, R and all dependencies are already included. Otherwise, you can add them manually:

```bash
conda install -c conda-forge r-base r-dplyr r-tidyr r-pheatmap r-rcolorbrewer
conda install -c bioconda bioconductor-complexheatmap
```

**Note**: The script is optimized for R 4.0 or later.

---

## System Requirements

### Core Dependencies (Required)

- **NCBI BLAST+** (version 2.12.0 or later)
- **Python 3.7+**
  - Pandas
  - Biopython

### Optional Dependencies

- **Git** (optional, for cloning/downloading the repository)
- **Mamba** (optional but highly recommended, for faster conda package installation)
- **Matplotlib** (optional, for visualization)
- **Seaborn** (optional, for visualization)
- **Krona Tools** (optional, for interactive charts) - Install from [Krona GitHub](https://github.com/marbl/Krona)
- **VFDB Database** - Download from [VFDB Official Site](https://www.mgc.ac.cn/VFs/)
- **R** (version 4.0+, optional, for advanced heatmap generation)
  - **ComplexHeatmap** (Bioconductor package, handles automated install if R is present)

#### Important Notes

- **Krona is optional** - the script will work without it, but won't generate interactive HTML charts
- If Krona is not installed, the script will gracefully skip HTML generation with a warning
- **Matplotlib/Seaborn are optional** - the script will work without them, but won't generate summary plots

---

## Database Setup

**Important information:**

The VFDB data is freely available under the Creative Commons Attribution-NonCommercial (CC BY-NC) license version 4.0 for personal and public non-commercial, research or academic use by individuals at academic, government or non-profit institutions. Users intending to use VFDB data for commercial purposes should contact the authors directly.

For more information, visit: <https://www.mgc.ac.cn/VFs/>

### Option 1: Automated VFDB Download Script (Recommended)

Use the enhanced download script with intelligent database management; note that the VFDB webpage may change, so the script may need to be updated.

```bash
# enter the directory
cd VF_classifier

# Download and create BLAST database with version checking
./VF_classifier.v0.1.0.py --setup

# To force an update even if up to date:
./VF_classifier.v0.1.0.py --setup --force
```

**Features of the enhanced script:**

- **Version comparison** between existing and new databases
- **Confirmation prompts** to prevent accidental overwrites
- **Comprehensive metadata** file creation

### Option 2: Manual Download

Download the VFDB dataset from the official repository (Approx. 10 MB). The database is updated regularly, so it's recommended to download the latest version.

```bash
# Download core dataset (Set B - nucleotide sequences)
wget https://www.mgc.ac.cn/VFs/Down/VFDB_setB_nt.fas.gz
gunzip VFDB_setB_nt.fas.gz
```

#### Create BLAST Database

Once you have the VFDB FASTA file, create a BLAST database:

```bash
makeblastdb -in VFDB_setB_nt.fas -dbtype nucl -out VFDB_db
```

```bash
# Move database files to a dedicated directory
mkdir -p database
mv VFDB_db* database/
mv VFDB_setB_nt.fas database/
```

**Important**: Keep the `VFDB_setB_nt.fas` file in the same directory as your BLAST database files. The script will automatically detect it, or you can specify its location with the `-v` parameter. You can put the database in any path you prefer, but make sure to update the path in the script accordingly.

---

## Verification

After installation, verify everything is working:

```bash
# Check BLAST installation
blastn -version

# Check Python packages
python -c "import pandas, biopython; print('Python packages OK')"

# Check optional packages
python -c "import matplotlib, seaborn; print('Visualization packages OK')"

# Check Krona (if installed)
ktImportText --help

# Check script
python VF_classifier.v0.1.0.py -h

# Check R and ComplexHeatmap (if R is installed)
Rscript -e "if (!require('ComplexHeatmap', quietly=TRUE)) { print('ComplexHeatmap not found, will be installed on first run') } else { print('ComplexHeatmap OK') }"
```

## Troubleshooting

### R Package Version Mismatches

If you see an error stating a package was "installed before R 4.0.0", it means R is finding an old, incompatible binary in your system or cache.

**Solution: The "Clean Start"**
If the standard installation fails with library errors, follow these steps to clear the cache and start fresh:

```bash
# Deactivate and remove the environment
conda deactivate
conda env remove -n VF_classifier

# CLEAR CACHE (The Nuclear Option)
# This ensures no broken binaries are reused
conda clean --all

# Re-create environment
mamba env create -f environment.yml
```

### Dependency Conflicts

If `conda` or `mamba` fails to solve the environment, try updating the solver first:

```bash
conda update -n base conda
# or if using mamba
mamba update -n base mamba
```

### Shadowing Libraries

Sometimes R finds packages in your personal folder (`~/R/...`) instead of the Conda environment. To see which paths R is searching, run:

```bash
conda activate VF_classifier
Rscript -e ".libPaths()"
```

If you see paths outside of your conda environment appearing first, you may need to temporarily rename those folders during analysis.
