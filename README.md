# BioPyView: Biopython-based Platform for Sequence Analysis

## Table of Contents:
- [Introduction](#introduction)
- [Features](#features)
- [License](#license)
- [Requirements](#requirements)
- [Installation](#installation)
- [Running the Application](#running-the-application)

## Introduction:
**What is BioPyView?**

BioPyView is a sequence analysis platform that focuses on working with unaligned and aligned biological sequences, supporting large files, and designed to be easy to use. The application primarily benefits from [Biopython](https://github.com/biopython/biopython).

This program was developed as part of my master's thesis at the Department of Computer Engineering, Boğaziçi University.

## Features:
- Supports common unaligned and aligned sequence file formats:
  - **Unaligned Formats:** FASTA, GenBank, EMBL, UniProt XML, Swiss-Prot, FASTQ, PIR, Nexus and more. 
  - **Aligned Formats:** Phylip, Clustal, Stockholm, MSF, Mase, emboss, maf, mauve and more.
  - **Full lists available at:** 
    - [Biopython SeqIO (unaligned formats)](https://biopython.org/wiki/SeqIO)
    - [Biopython AlignIO (aligned formats)](https://biopython.org/wiki/AlignIO)

- **Index Mode:** Open extremely large datasets beyond RAM limits.
- **Automatic Alignment:**
  - Muscle, ClustalOmega, Prank, Mafft, Probcons, MSAProbs
  - Needleman-Wunsch (Global) and Smith-Waterman (Local)
- **Manual Alignment:** Open, edit alignment files manually, or create new alignments from scratch.
- **Annotation Support:** View, edit, and save annotation-rich formats (GenBank, EMBL, UniProt-XML) and convert annotations between supported formats.
- **Motif Search:**
  - Exact and Approximate Search (PSSM-based)
  - Supports JASPAR, TRANSFAC, Pfm, and other common pattern files.
- **Phylogenetic Tree Features:**
  - Open tree files in FigTree
  - Tree Inference Methods: FastTree, RAxML
  - Tree file format conversion.
- **Biological Operations:** Concatenation, (Reverse) Complement, Transcription, Translation.
- **Editing Capabilities:** Undo-Redo, Copy-Paste, Delete/Insert Letters.
- **Memory-efficient handling of large files.**

## License:
- This application is licensed under **GPL-v3** and is **Biopython-based**.
- License file: [GPL-3.0](https://github.com/omerfarukcavass/biopyview?tab=GPL-3.0-1-ov-file)
- Biopython license file: [Biopython License](https://github.com/biopython/biopython?tab=License-2-ov-file)

## Requirements:
- **Python 3.10.11 or later** (lower versions may work but 3.10.11(or above) is recommended). Download from [Python official site](https://www.python.org)
- **Optional External Tools** (install any you need and configure the executable path in settings):

### Alignment Tools:
- [Muscle](https://drive5.com/muscle5/)
- [ClustalOmega](http://www.clustal.org/omega/)
- [Prank](https://ariloytynoja.github.io/prank-msa/)
- [Mafft](https://mafft.cbrc.jp/alignment/software/)
- [Probcons](http://probcons.stanford.edu/)
- [MSAProbs](https://msaprobs.sourceforge.net/homepage.htm#latest)
- [Needleman-Wunsch & Smith-Waterman](https://emboss.sourceforge.net/download/)

### Tree Tools:
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [FigTree](http://tree.bio.ed.ac.uk/software/Figtree/)


## Installation:

### General Instructions:
1. Download the content of this repository.
2. Unpack the downloaded zip file to a folder of your choice.

### Linux-based or MacOS:
1. Open a terminal.
2. Navigate to the project folder (replace `/path/to/project` with the actual path):
   ```
   cd /path/to/project
   ```
3. Grant execution permission to the installation script:
   ```
   chmod +x install.sh
   ```
4. Run the installation script:
   ```
   ./install.sh
   ```

### Windows:
1. Open File Explorer and navigate to the extracted folder.
2. Double-click `install.bat` to run it.
   - Alternatively, open Command Prompt (cmd), navigate to the project folder, and run:
     ```
     cd C:\path\to\project  # Replace with actual path
     install.bat
     ```

## Running BioPyView:

### Linux-based or MacOS:
1. Open a terminal.
2. Navigate to the project folder (replace `/path/to/project` with the actual path):
   ```
   cd /path/to/project
   ```
3. Grant execution permission to the run script (only required once):
   ```
   chmod +x run.sh
   ```
4. Run the program:
   ```
   ./run.sh
   ```

### Windows:
1. Open File Explorer and navigate to the extracted folder.
2. Double-click `run.bat` to start the application.
   - Alternatively, open Command Prompt (cmd), navigate to the project folder, and run:
     ```
     cd C:\path\to\project  # Replace with actual path
     run.bat
     ```

