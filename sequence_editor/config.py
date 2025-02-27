COLOR_MAPPING = {  # check https://www.bioinformatics.org/sms/iupac.html
    # DNA
    'A': 'green',
    'T': 'blue',
    'C': 'orange',
    'G': 'red',
    'N': 'lightgrey',  # Any nucleotide
    'Y': 'lightblue',  # Pyrimidine nucleotide (C or T)
    'K': 'lightgreen',  # Keto nucleotide (T or G)
    'M': 'lightcoral',  # Amino nucleotide (A or C)
    'W': 'lightskyblue',  # Weak nucleotide (A or T)
    'S': 'lightpink',  # Strong nucleotide (C or G)
    'H': 'lightseagreen',  # A, C, or T
    'B': 'lightyellow',  # C, G, or T
    'V': 'lightsalmon',  # A, C, or G
    'D': 'lightcyan',  # A, G, or T
    'X': 'lightsteelblue',  # Unknown nucleotide

    # RNA
    'U': 'purple',

    # Proteins (amino acids)
    'R': 'magenta',
    'E': 'olive',
    'Q': 'grey',
    'I': 'violet',
    'L': 'darkorange',
    'F': 'lightgreen',
    'P': 'lightpink',
    '*': 'black',  # Stop codon
    '-': 'grey',  # Gap
    '.': 'grey',  # Gap or indel
    '_': 'grey',  # Underscore (commonly used as a placeholder)
    '?': 'grey',  # Unknown symbol
}

# Add lowercase equivalents
COLOR_MAPPING.update({k.lower(): v for k, v in COLOR_MAPPING.items()})
