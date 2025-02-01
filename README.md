# Quadrupia
Quadrupia is a database of G-Quadruplex sequences from 108,449 organismal genomes across the Tree of Life. This repository contains scripts used for analysis of the G-Quadruplex data.

**Paper**: https://www.biorxiv.org/content/10.1101/2024.07.09.602008v1.full

## Directories

- **raw_data/** (Input data): Contains extracted G4 sequences. Please contact us for access.
- **metadata/**: Contains text files and tabular data with metadata such as taxonomic groups, ftp paths for download and other information associated with reference genomes.
- **scripts/**: Contains scripts used for downloading, parsing and pre-processing data.
- **slurm/**: Contains scripts used for running jobs on the cluster.
    - **batch/**: Bash scripts for submitting batch jobs to slurm.
    - **files/**: Text files containing command line arguments for slurm scripts in the **slurm/batch/** directory.
    - **prep/**: Python scripts for generating text files in the **slurm/files/** directory.
- **data/**: Contains extracted G4 sequences and GFF data in bed format for further processing, and is generated using scripts in the **scripts/** and the **slurm/** directories.
- **results/**: Contains results from processing of data in **data/** using scripts in the **scripts/** directory.
- **notebooks/**: Contains Jupyter notebooks that use the processed data in the **results/** directory for further analysis and data visualization.
- **supplementary_data/**: Contains supplementary data for the paper.

## Getting Started

To get started with the project, clone the repository and install the required dependencies. Optionally, create a new conda environment for the project.

```bash
# create a new conda environment (optional)
conda create -n quadrupia python=3.9
conda activate quadrupia

# clone repository and install dependencies
git clone https://github.com/Georgakopoulos-Soares-lab/Quadrupia.git
cd Quadrupia
pip install -r requirements.txt
```

## Contact
For questions, or to report bugs, contact -
```
nmc6088@psu.edu (Nikol Chantzi)
abn5461@psu.edu (Akshatha Nayak)
izg5139@psu.edu (Ilias Georgakopoulos-Soares)
```
