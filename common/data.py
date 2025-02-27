from collections import deque


class Data:
    def __init__(self, table_list, table_index, is_alignment: bool):
        """
        Table: a dictionary of SeqRecord objects.
               key: auto-increment id (printed alongside seq_id)
               value: SeqRecord object.
        :param table_list: for sequence files it is a list of one item. for alignment file, multiple tables.
        :param table_index:
        :param is_alignment
        :return:
        """
        self._table_list = table_list  # list of seq_record_dict
        self._table_index = table_index  # index of current seq_record_dict
        self._max_ids = [len(table) for table in self._table_list]  # auto-increment ids
        self.is_alignment = is_alignment
        self.number_of_tables = len(self._table_list)

    def __iter__(self):
        # Return an iterator for the dictionary items (in current table)
        table = self._table_list[self._table_index]
        return iter(table.items())

    def __reversed__(self):
        # Return a reversed iterator for the dictionary items in the current table
        table = self._table_list[self._table_index]
        return iter(reversed(table.items()))  # Create a reversed iterator

    def __len__(self):
        table = self._table_list[self._table_index]  # len of current table
        return len(table)

    def get_keys(self):
        # Return an iterator for the dictionary values (in current table)
        table = self._table_list[self._table_index]
        return iter(table.keys())

    def get_values(self):
        # Return an iterator for the dictionary values (in current table)
        table = self._table_list[self._table_index]
        return iter(table.values())

    def get_items(self):
        table = self._table_list[self._table_index]
        return table

    def iter_chunks(self, chunk_size=100):
        my_dict = self._table_list[self._table_index]
        it = iter(my_dict.items())
        for i in range(0, len(my_dict), chunk_size):
            yield {item[0]: item[1] for index, item in zip(range(chunk_size), it)}

    def get_new_id(self):
        # increment max_id and return
        max_id = self._max_ids[self._table_index]
        new_id = max_id + 1
        self._max_ids[self._table_index] = new_id
        return new_id

    def add_update_record(self, key, record):
        """
        updates if exist, add elsewhere.
        :param key:
        :param record:
        :return:
        """
        table = self._table_list[self._table_index]
        table[key] = record

    def get_record(self, key):
        table = self._table_list[self._table_index]
        return table.get(key)  # Using .get() for safe access

    def remove_record(self, key):
        table = self._table_list[self._table_index]

        if key in table:
            del table[key]
        else:
            print("no such key")

    def reorder(self, new_order):
        """
        Reorder dictionary based on seq viwer order before saving.
        Note: Only the sequences present in new_order are taken.
        New_order includes whatever is in the left panel. Deleted sequences will not appear in new_order.

        :param new_order:
        :return:
        """
        reordered_table = {key: self._table_list[self._table_index][key] for key in new_order}
        self._table_list[self._table_index] = reordered_table

    def switch_table(self, new_index):
        """
        Zero indexed.
        :return:
        """
        self._table_index = new_index

    def which_table(self):
        """
        Zero indexed.
        :return:
        """
        return self._table_index


class IndexData:
    def __init__(self, cache_limit):
        self.records = None
        self.query_id = None
        self.count = 0
        self.edited_dict = dict()
        self.cache_dict = dict()   # Stores cached sequences
        self.cache_queue = deque()  # Provides the limit for the cache.
        self.cache_limit = cache_limit

    def __len__(self):
        return len(self.records)

    def set_records(self, records):
        self.records = records  # index file

    def get_iter(self):
        for key in range(1, len(self) + 1):
            yield self.get_record(key)

    def get_record(self, key):
        # check edited
        record = self.edited_dict.get(key)
        if record is None:
            # check cache
            record = self.cache_dict.get(key)
            if record is None:
                # get from file
                self.query_id = key
                record = self.records[key]  # reads from file
        return record

    def id_to_key(self, identifier):
        """
        :param identifier: not used, we give id based on seq order.
        :return:
        """
        if self.query_id is not None:  # for retrieving (after creating dict).
            return self.query_id

        self.count += 1  # for filling dict
        return self.count

    def add_to_edited(self, key, seq_record):
        # remove if in cache
        if key in self.cache_dict:
            del self.cache_dict[key]
            self.cache_queue.remove(key)  # O(n) complexity (only called after an edit is performed, np)

        # add to edited
        self.edited_dict[key] = seq_record

        # print("d,q:",len(self.cache_dict), len(self.cache_queue))

    def add_to_cache(self, key, seq_record):
        # add if not in edited and not already added
        if self.cache_limit > 0 and key not in self.edited_dict and key not in self.cache_dict:
            # delete oldest if hit limit
            if len(self.cache_dict) >= self.cache_limit:
                oldest_key = self.cache_queue.popleft()
                # print(oldest_key)
                del self.cache_dict[oldest_key]

            # add record
            self.cache_queue.append(key)
            self.cache_dict[key] = seq_record

        # print("d,q:",len(self.cache_dict), len(self.cache_queue))
        # print(self.cache_queue)
