from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord

from sequence_editor.sequence_writer import filter_annots


class AlignmentWriter:
    def __init__(self, data, file_path, file_format, sel_annots):
        self.data = data
        self.file_path = file_path
        self.file_format = file_format
        self.sel_annots = sel_annots
        self.aligns = []

        self.create_align_data()

    def create_align_data(self):
        curr_table = self.data.which_table()
        no_of_table = self.data.number_of_tables

        for i in range(no_of_table):
            self.data.switch_table(i)
            records = [filter_annots(seq_record, self.sel_annots) for seq_record in self.data.get_values()]
            align = MultipleSeqAlignment(records)  # sequences are all the same length.
            self.aligns.append(align)

        self.data.switch_table(curr_table)  # switch back to initial table

    def write(self):
        """
        Create MultipleSeqAlignment object for each table and write them to file.
        :return:
        """
        with open(self.file_path, "w") as output_handle:
            AlignIO.write(self.aligns, output_handle, self.file_format)

        print(f"Alignment file written to {self.file_path}")

