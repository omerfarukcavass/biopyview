from Bio import SeqIO
import time
from Bio.SeqRecord import SeqRecord

from common.data import Data, IndexData
import gzip, bz2


binary_file_formats = ["abi", "abi-trim", "gck", "nib", "seqxml", "sff", "snapgene", "sff-trim", "twobit", "uniprot-xml", "xdna"]


class MySeqRecord:
    __slots__ = ['seq', 'description']

    def __init__(self, seq, description):
        self.seq = seq
        self.description = description

    def get_id(self):
        return self.description.split(" ")[0]

    def to_seq_record(self):
        return SeqRecord(seq=self.seq, id=self.get_id(), name=self.get_id(), description=self.description,
                         dbxrefs=[], features=[], annotations={},
                         letter_annotations={})


class SequenceReader:
    def __init__(self, file_path, file_format, zip_format):
        self.file_path = file_path
        self.file_format = file_format
        self.zip_format = zip_format

        # to store sequences and names
        self._seq_records = dict()  # key: integer id starts from 1, value: seq record

    def read_seq_file(self):
        # Read the file and store sequences in the list
        data = self.read_sequences(self.file_path)
        return data

    def read_sequences(self, file_path):
        start_time = time.time()

        # Read file and store sequences in the dictionary
        try:
            # zip format
            if self.zip_format == "-":
                open_fun = open
            elif self.zip_format == "gz":
                open_fun = gzip.open
            elif self.zip_format == "bz2":
                open_fun = bz2.open
            else:
                raise ValueError("Undefined zip format !")

            # binary or text
            if self.file_format in binary_file_formats:  # binary files
                open_mode = "rb"
            else:
                open_mode = "rt"

            # fasta/fasta-2line or other format
            if self.file_format == "fasta":
                read_fun = self.read_fasta
            elif self.file_format == "fasta-2line":
                read_fun = self.read_fasta_2line
            else:
                read_fun = self.read

            with open_fun(file_path, open_mode) as file_handle:
                read_fun(file_handle)

        except Exception as e:
            print(e)
            raise  # raise e

        data = Data([self._seq_records], 0, False)

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for reading {len(self._seq_records)} seqs.")

        return data

    def read(self, file):
        for index, record in enumerate(SeqIO.parse(file, self.file_format)):
            self._seq_records[index + 1] = record

    def read_fasta(self, file):
        from Bio.SeqIO.FastaIO import SimpleFastaParser
        for index, (title, seq) in enumerate(SimpleFastaParser(file)):
            record = MySeqRecord(seq, title)
            self._seq_records[index + 1] = record

    def read_fasta_2line(self, file):
        from Bio.SeqIO.FastaIO import FastaTwoLineParser
        for index, (title, seq) in enumerate(FastaTwoLineParser(file)):
            record = MySeqRecord(seq, title)
            self._seq_records[index + 1] = record


class SequenceIndexReader:
    def __init__(self, file_path, file_format, cache_limit):
        self.file_path = file_path
        self.file_format = file_format
        self.cache_limit = cache_limit

        # to store sequences and names
        self._seq_records = dict()  # key: integer id starts from 1, value: seq record
        self.max_length = 0

    def read_sequences(self):
        start_time = time.time()

        data = IndexData(self.cache_limit)
        try:
            # As of Biopython 1.69, this will preserve the ordering of the records in
            #     file when iterating over the entries.
            records = SeqIO.index(self.file_path, self.file_format, key_function=data.id_to_key)  # returns SeqRecord dict
            data.set_records(records)
        except Exception as e:
            print(e)
            raise  # raise e

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for indexing {len(data)} seqs.")

        return data

