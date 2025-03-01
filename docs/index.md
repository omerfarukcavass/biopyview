# User Manual

## Table of Contents
- [File](#file)  
  - [Open Sequence File](#open-sequence-file)  
  - [Open Sequence File as Index](#open-sequence-file-as-index)  
  - [Open Alignment File](#open-alignment-file)  
  - [Save File](#save-file)  
  - [Convert Sequence File Format](#convert-sequence-file-format)  
  - [Convert Alignment File Format](#convert-alignment-file-format)  
  - [Settings](#settings)   
- [Edit](#edit)  
  - [Editor Mode (Keyboard)](#editor-mode-keyboard)
  - [Undo-Redo](#undo-redo)
  - [Copy-Paste](#copy-paste)
  - [Delete-Insert](#delete-insert)
- [Select](#select)  
- [View](#view)  
  - [Color Mode](#color-mode)
  - [Font Size](#font-size)
  - [Navigate to Position](#navigate-to-position)
- [Sequence](#sequence)
  - [Add New Sequence](#add-new-sequence)
  - [Concatenate Selected](#concatenate-selected)
  - [Complement Selected](#complement-selected)
  - [Reverse Complement Selected](#reverse-complement-selected)
  - [Transcribe Selected](#transcribe-selected)
  - [Translate Selected](#translate-selected)
- [Alignment](#alignment)
  - [Switch Alignment](#switch-alignment)
  - [Multiple Sequence Alignment (MSA)](#multiple-sequence-alignment-msa)
  - [Progress Display](#progress-display)
  - [Pairwise Alignment](#pairwise-alignment)
  - [Manual Alignment](#manual-alignment)
- [Annotation](#annotation)
  - [Supported Annotation Formats](#supported-annotation-formats)
  - [Annotation Window](#annotation-window)
  - [Saving Changes](#saving-changes)
- [Phylogenetic Tree](#phylogenetic-tree)
  - [Open Tree](#open-tree)
  - [Infer Tree](#infer-tree)
  - [Convert File Format](#convert-file-format)
  - [Tree Tools](#tree-tools)
- [Motif](#motif)
  - [Create Motif](#create-motif)
  - [Read Motif](#read-motif)
  - [Motif Analysis](#motif-analysis)
  - [Save Motif](#save-motif)
  - [Motif Search](#motif-search)
  - [Search Options](#search-options)
- [Other Tools](#other-tools)


---

## File
### Open Sequence File

This feature allows users to open a sequence file. The application identifies the file as a **Sequence**, meaning it contains unaligned sequences rather than an alignment file. Each sequence is treated independently and can be of any type (DNA, RNA, or Protein) with varying lengths.

The tool supports 34 file formats, including:  
`abi, abi-trim, ace, fasta, fasta-2line, ig, embl, embl-cds, gb, gck, genbank, genbank-cds, imgt, nib, cif-seqres, cif-atom, pdb-atom, pdb-seqres, phd, pir, fastq, fastq-sanger, fastq-solexa, fastq-illumina, qual, seqxml, sff, snapgene, sff-trim, swiss, tab, twobit, uniprot-xml, xdna`.  

Additionally, compressed formats (`.gz`, `.bz2`) are supported.

#### Sequence File Dialog Window
<img src="images/seq-file-dialog.png" width="50%" alt="Sequence file dialog window">


Once opened, the file is displayed in the application.  

- **Left Panel:** Shows sequence names or IDs.  
- **Right Panel:** Displays sequence content, with a ruler marking every 10th position.  
- **Bottom Panel:** Indicates file type (`Sequence` or `Alignment`). Some editing options are specific to each type. The right side shows real-time row and column values.

#### Sequence Viewer Window  
<img src="images/seq-opened.png" width="80%" alt="Sequence file dialog window">

### Open Sequence File as Index

Index mode allows the user to open large datasets that exceed RAM limitations. However, editing capabilities are restricted in this mode (As you can see, some menu items are disabled). 

The supported formats include 19 file types:  
`ace, embl, fasta, fastq, fastq-sanger, fastq-solexa, fastq-illumina, genbank, gb, ig, imgt, phd, pir, sff, sff-trim, swiss, tab, qual, uniprot-xml`.

### Open Alignment File

This feature allows the user to open an alignment file and sets the file type to **Alignment**. Alignment files require sequences to have a fixed size, as they are part of an alignment. The alignment sequences are displayed the same way as sequence files.

No zip formats are supported for alignment files. The supported formats are:  
`clustal, emboss, fasta, fasta-m10, ig, maf, mauve, msf, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm`.

Typically, a file contains a single alignment. However, if the file contains multiple alignments, a selector is provided to choose an alignment based on its order in the file.

#### Alignment Selector  
<img src="images/select-align.png" width="50%" alt="Alignment selector">


### Save File

This option allows the user to save the currently worked file. During the save process, the user is prompted to select which annotations to keep. Annotations that are not supported by the selected format will be omitted.

#### Supported Formats for Sequence Files:  
`fasta, fasta-2line, gb, genbank, embl, imgt, nib, phd, pir, fastq, fastq-sanger, fastq-solexa, fastq-illumina, qual, seqxml, sff, tab, xdna`.

#### Supported Formats for Alignment Files:  
`clustal, maf, mauve, nexus, phylip, phylip-sequential, phylip-relaxed, stockholm`.

#### Save Window  
<img src="images/save-dialog.png" width="50%" alt="Save window">

### Convert Sequence File Format

This option allows the user to convert the format of a sequence file without the need to open it. The process is memory-efficient, as records are read and written one by one, without loading the entire file into memory. This ensures that large files can be converted without overwhelming system resources.

The input formats for this option are the same as those in the **Open Sequence File** option. The output formats available for conversion are the same as those supported in the **Save File** option for sequence files.

#### Convert Sequence Window  
<img src="images/convert-seq.png" width="50%" alt="Save window">

### Convert Alignment File Format

This option allows the user to convert an alignment file from one format to another without opening it.
The input formats for this option are the same as those in the **Open Alignment File** option. The output formats available for conversion are the same as those supported in the **Save File** option for alignment files.

#### Convert Alignment Window  
<img src="images/convert-align.png" width="50%" alt="Save window">

### Settings

This menu item opens a configuration interface with several tabs: **General**, **Alignment Tools**, **Tree Tools**, and **Other Tools**.

- **General**: This tab provides general application settings. The available options are explained throughout the document as related sections appear.

<img src="images/set-general.png" width="70%" alt="General settings">

- **Alignment Tools**: In this tab, the user can configure external alignment tool parameters. Currently, the application supports interfaces to these tools for multiple sequence alignment: **Muscle**, **ClustalOmega**, **Prank**, **Mafft**, **Probcons**, **MSAProbs**. For pairwise alignment, **Needleman-Wunsch** and **Smith-Waterman** from the EMBOSS package are available. Users can specify input, output, and output file format for MSA, as well as gap open and gap extension parameters for pairwise alignments.

<img src="images/set-align.png" width="70%" alt="Alignment tools settings">

- **Tree Tools**: This tab allows users to configure external tree tool parameters. The application supports **RAxML** and **FastTree** for phylogenetic tree inference, with **FigTree** for visualization. Users can specify input and output file parameters, along with additional tool-specific options.

<img src="images/set-tree.png" width="70%" alt="Tree tools settings">

- **Other Tools**: In this tab, users can add custom external tools by entering their respective commands. This provides immediate access to the user's preferred bioinformatics programs, which can be executed later via the **Other Tools** menu.
 
<img src="images/set-other.png" width="70%" alt="Other tools settings">


## Edit

### Editor Mode (Keyboard)

The application allows users to toggle editing mode. When editing mode is enabled, users can directly replace letters in the sequence using keystrokes. This includes the ability to insert characters from the DNA, RNA, or protein alphabets, as well as the indel-gap character (-).

### Undo-Redo

The application allows users to undo and redo any changes made during sequence editing. All actions are tracked, but the number of actions retained can be restricted to conserve memory. This can be configured in the settings: **Maximum number of operations retained on the undo stack**. By default, there is no limit on the number of actions that can be undone. Users can also **save the current state** of the editor, which locks in the changes made up to that point, preventing any further undos before the saved state.

### Copy-Paste

The application offers the following options for copying and pasting sequence data:

- **Copy Selection as FASTA**: Copies the selected sequences or regions in FASTA format, including sequence headers.
- **Copy Selection as Text**: Copies only the sequence letters (no headers) in plain text format.
- **Paste FASTA**: Pastes sequences from the clipboard in FASTA format directly into the application, enabling easy integration of external sequence data.

### Delete-Insert

The application provides the following options for deleting and inserting sequence data:

- **Delete Selected**: Deletes the entire selected sequence.
- **Delete Selected Letters**: Deletes only the selected letters within the sequence.
- **Insert Unknown Letter**: Inserts an unknown letter (`?`) at the selected position in the sequence.

## Select

Selection can be made using the mouse across three distinct panels: the **left panel** (sequence names), the **right panel** (sequence content), and the **top panel** (ruler).

- **Left Panel (Sequence Names):** Dragging across sequence names selects entire sequences, including all their letters.  
- **Right Panel (Sequence Content):** Dragging across this panel selects only a specific segment of the sequences, allowing for precise control over the selection range.  
- **Top Panel (Ruler):** Dragging across the ruler selects entire columns (sites), enabling the selection of specific alignment positions across all sequences.  

An example selection is shown below:

<img src="images/selection.png" width="70%" alt="Example selection">

The following options are also available for selection:

- **De-select**: Removes the current selection.
- **Select All**: Selects all sequences.
- **Expand Selection Right**: Expands the selection to the right.
- **Expand Selection Left**: Expands the selection to the left.
- **Expand Selection Up**: Expands the selection upward.
- **Expand Selection Down**: Expands the selection downward.
- **Move Selected Sequences Up**: Moves the selected sequences upward.
- **Move Selected Sequences Down**: Moves the selected sequences downward.
- **Move Selected Sequences Top**: Moves the selected sequences to the top.
- **Move Selected Sequences Bottom**: Moves the selected sequences to the bottom.

## View

### Color Mode
Switch between color modes: no color, foreground, and background. Coloring is performed based on the letter, with each letter having a unique color.

<img src="images/colored.png" width="70%" alt="Colored sequence view">

### Font Size
Allows adjustment of font size in the range of 8-32, using the Courier font.

### Navigate to Position
Allows the user to go to an arbitrary position, specified by the line and column number.


## Sequence

The Sequence menu allows users to manage and manipulate nucleotide and protein sequences. You can add new sequences, concatenate, complement, reverse complement, transcribe, and translate sequences with various configuration options.

### Add New Sequence
Users can add a new sequence by manually entering its sequence ID and content, or by opening a file to append sequences. This provides flexibility for handling new data.

### Concatenate Selected
Combines selected sequences into a single sequence.

### Complement Selected
Computes the complement of selected DNA or RNA sequences by replacing each nucleotide with its complement.

### Reverse Complement Selected
Generates the reverse complement of selected DNA or RNA sequences by reversing the sequence and replacing each nucleotide with its complement.

### Transcribe Selected
Converts nucleotide sequences into RNA, either from the coding or template strand.

### Translate Selected
Translates nucleotide sequences into protein sequences. Users can configure the following options:
  
<img src="images/translate.png" width="50%" alt="Translate window">

  - **Genetic Code Table**: Select the translation table based on NCBI's genetic code tables (supported by Biopython).
  - **Stop at First In-Frame Stop Codon**: Stops translation at the first in-frame stop codon.
  - **Complete Coding Sequence (CDS)**: If the nucleotide sequence is a complete CDS, select this option for validation.


## Alignment

This menu offers both automated and manual sequence alignment features.

### Switch Alignment: 
Switch between different alignments in a file (e.g., Stockholm format) if the file contains multiple MSAs.

### Multiple Sequence Alignment (MSA)
The app supports external MSA tools like Muscle, ClustalOmega, Prank, Mafft, Probcons, and MSAProbs. You can configure tool parameters and access official websites for tool installation. 

  - [Muscle](https://drive5.com/muscle5/)
  - [ClustalOmega](http://www.clustal.org/omega/)
  - [Prank](https://ariloytynoja.github.io/prank-msa/)
  - [Mafft](https://mafft.cbrc.jp/alignment/software/)
  - [Probcons](http://probcons.stanford.edu/)
  - [MSAProbs](https://msaprobs.sourceforge.net/homepage.htm#latest)
  - [Needleman-Wunsch & Smith-Waterman](https://emboss.sourceforge.net/download/)

The window for the **Muscle** tool:

<img src="images/muscle.png" width="50%" alt="Muscle tool window">

### Progress Display 
- The alignment progress is shown in real-time, capturing the tool's output via pipes.

<img src="images/progress.png" width="50%" alt="Muscle tool progress window">

### Pairwise Alignment 
EMBOSS tools like **Needleman-Wunsch** (Global) and **Smith-Waterman** (Local) are available for pairwise alignments. Gap settings can be adjusted for fine-tuning. A progress window is also shown during execution.

### Manual Alignment

Refine alignments by directly editing sequences. The following operations are available:
  - **Move Selected Position Right**: Shifts selected letters to the right if space is available.
  - **Move Selected Position Left**: Shifts selected letters to the left if space is available.
  - **Insert Gap Move Right**: Inserts a gap and shifts sequences to the right.
  - **Insert Gap Move Left**: Inserts a gap and shifts sequences to the left.
  - **Replace Selected With Gap(s)**: Replaces selected letters with gaps.
  - **Delete Gap-Only Columns**: Deletes columns consisting solely of gaps.


## Annotation

The annotation menu of the app allows users to view and modify sequence annotations. Annotations in sequence files provide **additional information about the sequence, such as functional features, structural elements, and metadata**, which enhance understanding and analysis of biological data. The platform supports annotations from various formats and features a dedicated **Annotation Window**, enabling users to edit general information, key-value pairs, feature tables, letter annotations, and cross-references to external databases.

### Supported Annotation Formats

Annotation support varies across file formats. Simple formats like **FASTA** have limited annotations, while formats like **GenBank**, **EMBL**, or **UniProt-XML** contain rich annotations. **FASTQ** and **Stockholm** files support letter annotations.

An example FASTA file:

<img src="images/fasta.jpg" width="70%" alt="An example FASTA file">

An example GenBank file:

<img src="images/genbank.jpg" width="70%" alt="An example GenBank file">


### Annotation Window

The application provides a dedicated **Annotation Window** that allows users to inspect and modify annotations of the selected sequence.

The **Annotation Window** consists of five tabs:

#### **Information Tab**
Displays general information about the sequence, allowing users to edit fields using menu items.

- **Id:** The primary identifier for a sequence, usually an accession number.
- **Name:** An alternative identifier, such as a clone name or alias.
- **Description:** A detailed description providing context about the sequence.

<img src="images/info-tab.png" width="70%" alt="Information tab">

#### **Annotations Tab**
A two-column table storing additional metadata about the sequence. Users can add, delete, or modify key-value pairs.

<img src="images/annotation-tab.png" width="70%" alt="Annotation tab">

Annotations shown in the GenBank file:

<img src="images/annot-all-c.png" width="70%" alt="Annotations in a GenBank file">

The annotations displayed in the app’s Annotation Window:

<img src="images/app-annot-all-c.png" width="70%" alt="Annotations displayed in the app">

#### **Features Tab**
Displays and edits feature table elements, which annotate specific regions of a sequence.

<img src="images/app-features-tab.png" width="70%" alt="Features tab in the app">

Each feature represents a biologically significant function. Users can edit feature properties, including **CDS (coding sequences)**.

CDS feature in the GenBank file:

<img src="images/features-cds-c.png" width="70%" alt="CDS feature in the GenBank file">

CDS feature in the app:

<img src="images/app-features-cds-c.png" width="70%" alt="CDS feature in the app">

Feature locations can be edited through a dedicated location editor.

<img src="images/app-features-tab-location.png" width="70%" alt="Editing location in features tab">

#### **Letter Annotation Tab**
Used for per-letter annotations, such as **quality scores** (from FASTQ files) or **secondary structure data** (from Stockholm files).

<img src="images/letter-annot-tab.png" width="70%" alt="Letter annotation tab">


#### **Database Cross-References Tab**
Lists cross-references to external databases, allowing users to link sequences to external resources.

<img src="images/database-tab.png" width="70%" alt="Database cross-references tab">


### Saving Changes
After modifying annotations, users must save their changes. Unlike edits in the **Sequence Viewer**, changes made in the **Annotation Window** cannot be undone after saving, so users should review edits carefully before confirming them.


## Phylogenetic Tree

This menu allows for parsing, visualizing, and inferring phylogenetic trees from sequence data. Users can open existing trees, infer trees from alignments, and convert between various tree formats. It integrates with external tools like **FigTree**, **FastTree**, and **RAxML**, providing a user-friendly interface. The real-time progress of tree tools is displayed in the application, similar to alignment tools.

### Open Tree: 
Visualize phylogenetic trees using **FigTree**.
  
<img src="images/figtree.png" width="50%" alt="FigTree window">

### Infer Tree 
Perform tree inference from alignments using **FastTree** or **RAxML**.

<img src="images/fasttree.png" width="50%" alt="FastTree window">

### Convert File Format: 
Converts tree files between formats such as `newick`, `nexus`, `phyloxml`, and `nexml`. 
  
<img src="images/tree-convert.png" width="50%" alt="Tree convert window">

### Tree Tools:
- [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/index.html)
- [FastTree](http://www.microbesonline.org/fasttree/)
- [FigTree](http://tree.bio.ed.ac.uk/software/Figtree/)


## Motif

Motif search identifies biologically significant patterns in sequences. This component allows users to **create, read, edit, save, and search motifs**. It supports both **instance-based and count matrix-based** motif creation, offers detailed motif analysis through Position-Weight Matrices (PWM) and **Position-Specific Scoring Matrices (PSSM)**, and provides flexible **search options**.

### **Create Motif**
Users can create motifs using two methods:

- **Motif Instances**: Enter motif sequences manually.

<img src="images/motif-instance.png" width="70%" alt="Motif creation using instances">

- **Count Matrix**: Define motif probabilities using a matrix.
  
<img src="images/motif-counts.png" width="70%" alt="Motif creation using count matrix">

### **Read Motif**
Motifs can be loaded from various file formats, including **AlignAce, MEME, TRANSFAC, JASPAR, and PFM**. Loaded motifs can be viewed and edited using **instances** or **count matrices**.

### **Motif Analysis**
The **motif window** contains four key tabs for analysis:

- **Instances**: Displays motif sequences.
- **Counts**: Shows frequency of nucleotides at each position.

<img src="images/count-matrix.png" width="70%" alt="Count matrix">

- **PWM (Position-Weight Matrix)**: Normalized matrix representation.

<img src="images/pwm-matrix.png" width="70%" alt="PWM matrix">

- **PSSM (Position-Specific Scoring Matrix)**: Log-odds matrix used in searches.

<img src="images/pssm-matrix.png" width="70%" alt="PSSM matrix">

Additional options:
- **Pseudocounts**: Adjusts values in PWM to prevent zero probabilities.
- **Background Probabilities**: Allows custom nucleotide frequencies.

### **Save Motif**
Motifs can be saved in formats like **ClusterBuster, PFM, JASPAR, and TRANSFAC** (only for fixed-length motifs).

### **Motif Search**
Users can search motifs in sequence data using two methods:

- **Exact Search**: Finds perfect motif matches.
- **PSSM-Based Search**: Identifies approximate matches based on similarity scores.

Users can **navigate results** using "Next" and "Previous" buttons, and **minimize the search panel** for better sequence viewing.

<img src="images/motif-min.png" width="70%" alt="Minimized search window">

### **Search Options**
The search tool provides flexible options:
- Search within **all sequences** or **selected sequences**.
- Start from **the first sequence** or **the currently selected sequence**.

<img src="images/search-motif.png" width="70%" alt="Motif search options">


## Other Tools
Tools added in the settings window are listed under this menu with the given names, providing **easy access to frequently used external programs**.

## Troubleshooting
If you encounter issues, open an [issue](https://github.com/omerfarukcavass/biopyview/issues).
