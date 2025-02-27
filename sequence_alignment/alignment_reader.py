from Bio import AlignIO
import time
from common.data import Data


class AlignmentReader:
    def __init__(self, file_path, file_format):
        self.file_path = file_path
        self.file_format = file_format

        # to store sequences and names
        self._aligns = list()  # key: int id starts from 1, value: seq record

    def read_align_file(self):
        # Removes all elements from the list
        self._aligns.clear()

        # Read the file and store sequences in the list
        data = self.read_aligns(self.file_path)

        return data

    def read_aligns(self, file_path):
        start_time = time.time()

        try:
            with open(file_path, "r") as file:
                for alignment in AlignIO.parse(file, self.file_format):
                    align_dict = dict()

                    for index, record in enumerate(alignment):
                        align_dict[index + 1] = record

                    self._aligns.append(align_dict)
        except Exception as e:
            print(e)
            raise  # raise e

        data = Data(self._aligns, 0, True)  # user select index.

        end_time = time.time()
        runtime = end_time - start_time
        print(f"Runtime: {runtime} seconds for reading {len(self._aligns)} aligns.")

        return data
