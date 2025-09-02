# InteractionFragHub: Fragment Libraries from Protein-Ligand Interactions: Uncovering New Opportunities for Targeted Drug Design

![logo](static/picture/logo.png)
<br>This repository contains the code for setting up and implementing the InterFragHub local platform, a web-based tool for visualizing and interacting with fragment-based molecular design.The platform provides an interactive interface to facilitate the assembly, optimization, and analysis of molecular fragments.

Key features of the InterFragHub platform include:
    (1)Break Moleucles;
    (2)Split PDB Complex;
    (3)Generate Interaction Fragments;
    (4)Search Fragments;
    (5)Optimize molecules based on Fragments.

# How to install
To run InterFragHub, you need to:
* Clone this repository
* Setup Python and third-party dependencies.

## Cloning this repository
See "Clone" button on this page for further information

## Setting up Python and third-party dependencies
Our package has been developed and tested using Python 3.8 and the following
versions of the third-party packages:
* python==3.8.20
* gradio==4.44.1
* rdkit==2024.3.5
* biopython==1.83
* biopandas==0.5.1
* scipy==1.10.1

Since InterFragHub involves multiple parts, it is necessary to creat and configure the other two independent virtual environment.
* ESP-DNN:
In "Generate EC Fragments":
    Use https://github.com/AstexUK/ESP_DNN for generate pqr files for ligand and protein.
    Use https://github.com/AstexUK/esp-surface-generator for generate tmesh files for ligand and protein.
* Link-INVENT:
    See https://github.com/MolecularAI/Reinvent for more information
## Install ECscore conda environment and start up the local website
    * cd InteractionFragHub
    * conda env create -f environment.yml
    * conda activate ECscore
    * python localWebsite_interface.py
    
    
    
