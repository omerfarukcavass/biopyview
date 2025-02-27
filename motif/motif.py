from Bio.motifs import Motif, Instances
# https://biopython.org/docs/1.83/api/Bio.motifs.html?highlight=motifs#module-Bio.motifs

DNA_MOTIF = "DNA_MOTIF"
RNA_MOTIF = "RNA_MOTIF"
AA_MOTIF = "AA_MOTIF"


class MotifSearch:
    def __init__(self, motifs: list, search_in_selected=False, start_from_selected=False, pssm_search=False, pssm_threshold = None):
        self.motifs = motifs
        self.search_in_selected = search_in_selected
        self.start_from_selected = start_from_selected
        self.pssm_search = pssm_search
        self.pssm_threshold = pssm_threshold

    def start_search(self):
        pass

    def stop_search(self):
        pass

    def next_item(self):
        pass

    def prev_item(self):
        pass


class FixedLengthMotif(Motif):
    def __init__(self, alphabet, instances=None, counts=None):  # put either instances or counts
        if instances is not None:
            super().__init__(alphabet=alphabet, counts=counts, instances=Instances(instances, alphabet=alphabet))
        elif counts is not None:
            super().__init__(alphabet=alphabet, counts=counts)


class VariableLengthMotif:
    def __init__(self, alphabet, instances):
        self.alphabet = alphabet
        self.instances = instances  # a list of Seq/str objects.

# def is_motif_nucleotide(motif):
#     if isinstance(motif, Motif):
#         return len(motif.alphabet) == 4  # A', 'C', 'G', 'T' (for DNA), or 'A', 'C', 'G', 'U' (for RNA)
#     elif isinstance(motif, VariableLengthMotif):
#         return motif.is_nucleotide

