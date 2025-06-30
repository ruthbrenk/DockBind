# DockBind: A High-Quality Dataset for Training Machine Learning Scoring Funtions

**DockBind** is a curated dataset designed to support the development and benchmarking of machine learning models for protein–ligand binding prediction. 
It includes high-quality protein-ligand complexes modeled through molecular docking annotated with binding affinity data from ChEMBL and BindingDB.

---

## Overview


DockBind addresses key limitations of existing datasets by:
- Merging structural data (Binding MOAD and the PDB) with curated bioactivity data (ChEMBL, BindingDB)
- Providing a simple and robust pipeline for docking pose selection

---

## Dataset Contents

The dataset includes:

- **Protein–ligand complexes** modeled with molecular docking
- **Ligands with binding data** (K<sub>d</sub>, K<sub>d</sub>, and IC<sub>50</sub>) sourced from BindingDB and ChEMBL

## Download DockBind scripts

```bash
# Clone the repository
git clone https://github.com/ruthbrenk/DockBind.git
cd DockBind

# Install dependencies
conda env create -f environment.yml
```

## Data availability
DockBing data can be downloaded from: [10.5281/zenodo.15768136](https://doi.org/10.5281/zenodo.15768136)

- DockBind v1 protein-ligand complexes
- Docking validation poses
- A complementary approach for template-ligand pairs selection based on a ML approach: ML\_for\_template\_prediction\_data.tar.gz
