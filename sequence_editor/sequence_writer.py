from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from common.data import Data, IndexData
from sequence_editor.sequence_reader import binary_file_formats


class SequenceWriter:
    def __init__(self, data: Data, file_path, file_format, sel_annots):
        self.data = data
        self.file_path = file_path
        self.file_format = file_format
        self.sel_annots = sel_annots

    def write(self):
        # binary or text
        if self.file_format in binary_file_formats:  # binary files
            save_mode = "wb"
        else:
            save_mode = "wt"

        # Save SeqRecords to the chosen file location
        with open(self.file_path, save_mode) as output_handle:
            data = filter_annots_generator(self.data.get_values(), self.sel_annots)
            SeqIO.write(data, output_handle, self.file_format)
        print(f"Sequence file written to {self.file_path}")


class SequenceIndexWriter:
    def __init__(self, data: IndexData, file_path, file_format):
        self.data = data
        self.file_path = file_path
        self.file_format = file_format

    def write(self):
        # binary or text
        if self.file_format in binary_file_formats:  # binary files
            save_mode = "wb"
        else:
            save_mode = "wt"

        # Save SeqRecords to the chosen file location
        with open(self.file_path, save_mode) as output_handle:
            SeqIO.write(self.data.get_iter(), output_handle, self.file_format)

        print(f"Sequence index file written to {self.file_path}")


def filter_annots_generator(records, sel_annots):
    for record in records:
        if not isinstance(record, SeqRecord):
            record = record.to_seq_record()
        yield filter_annots(record, sel_annots)


def filter_annots(seq_record: SeqRecord, sel_annots):
    if not sel_annots["Name"]:
        seq_record.name = ""
    if not sel_annots["Desc"]:
        seq_record.description = ""
    if not sel_annots["Annot"]:
        seq_record.annotations = {}
    if not sel_annots["Feature"]:
        seq_record.features = []
    if not sel_annots["Letter Annot"]:
        seq_record.letter_annotations = {}
    if not sel_annots["Dbxrefs"]:
        seq_record.dbxrefs = []
    return seq_record
