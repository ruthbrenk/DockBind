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

The available scripts include the data selection and filtering, and Jupyter notebooks used for the results analysis


```bash
# Clone the repository
git clone https://github.com/ruthbrenk/DockBind.git
cd DockBind

# Install dependencies
conda env create -f environment.yml
```

## Download DockBind
DockBind data can be downloaded from: [10.5281/zenodo.15768136](https://doi.org/10.5281/zenodo.15768136)

- **DockBind v1 protein-ligand complexes**: protein-ligand complexes included in the DockBind dataset, and the csv file annotated with binding affinity data (DockBind\_complexes\_affinity.tar.gz)
- **DockBind raw files**: all the generated poses with molecular docking, pdb files, reference ligands, and the selected PDB structures released after 2020 (DockBind\_raw\_data.tar.gz)
- **Docking validation poses**: docking poses and datasets used to validate the DockBind workflow (docking\_validation\_data.tar.gz)
- **A complementary approach for template-ligand pairs selection based on a ML approach**: datasets used to test the machine learning classifiers for template-ligand pair selection (ML\_for\_template\_prediction\_data.tar.gz)
