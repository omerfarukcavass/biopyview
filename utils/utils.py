import os

from common.configs import MAX_FONT_SIZE, MIN_FONT_SIZE
from sequence_editor.config import COLOR_MAPPING

DNA_BASES = "ACGT"
RNA_BASES = "ACGU"
AA_BASES = 'ACDEFGHIKLMNPQRSTVWY'


def toggle_menu_on_off(menu, menu_items, on = False):
    """If app is closed when child window is open, it tries to access menu of parent for enabling but menu is gone  """
    try:
        state = "normal" if on else "disabled"
        for name in menu_items:
            menu.entryconfig(name, state=state)
    except Exception as e:
        # print("EXCEPTION: toggle_menu_on_off_all ", e)
        return


def toggle_menu_on_off_all(menu, on = False):
    """If app is closed when child window is open, it tries to access menu of parent for enabling but menu is gone  """
    try:
        menu_items = menu.index("end")
        state = "normal" if on else "disabled"

        for i in range(menu_items + 1):
            if menu.type(i) == "separator":
                continue  # Skip separators
            menu.entryconfig(i, state=state)
    except Exception as e:
        # print("EXCEPTION: toggle_menu_on_off_all ", e)
        return


def print_seq_record(seq_record):
    print("==== SEQ RECORD ====")
    print(f"ID: {seq_record.id}")
    print(f"Name: {seq_record.name}")
    print(f"Description: {seq_record.description}")

    print("\nAnnotations:")
    for key, value in seq_record.annotations.items():
        print(f"  {key}: {value}")

    print("\nLetter Annotations:")
    for index, (seq_char, *annotations) in enumerate(
            zip(seq_record.seq, *seq_record.letter_annotations.values()), 1):
        print(f"  Index {index}: {seq_char}, Annotations: {annotations}")

    print("\nDatabase cross-references:")
    for dbxref in seq_record.dbxrefs:
        print(f"  {dbxref}")
    print("==== ==== ==== ====")


def print_all_table_items(tree, count=3):
    """
    Print items in given tree/table.
    """
    print("Items in table:")
    all_items = tree.get_children()
    for item_id in all_items[:count]:
        item_values = tree.item(item_id, 'values')
        print(f"Item {item_id}: {item_values}")

def seperator(title):
    print(f"======= {title} =======")


def is_valid_letter(letter:str):
    return letter.upper() in COLOR_MAPPING.keys()


def is_dna(seq:str):
    seq_set = set(seq.upper())

    # Check if all characters are DNA bases
    if seq_set.issubset(set(DNA_BASES)):
        return True
    return False


def is_rna(seq:str):
    seq_set = set(seq.upper())

    # Check if all characters are RNA bases
    if seq_set.issubset(set(RNA_BASES)):
        return True
    return False


def is_nucleotide(seq:str):
    # Check if all characters are DNA or RNA bases
    if is_dna(seq):
        return True
    elif is_rna(seq):
        return True
    return False


def is_amino_acid(seq:str):
    seq_set = set(seq.upper())

    # Check if any characters are not DNA/RNA bases but are amino acids
    if seq_set.issubset(set(AA_BASES)):
        return True
    return False


def safe_index(lst: list, element):
    try:
        return lst.index(element)
    except ValueError:
        return None


def check_font_size_range(font_size: int):
    if font_size > MAX_FONT_SIZE:
        font_size = MAX_FONT_SIZE
    elif font_size < MIN_FONT_SIZE:
        font_size = MIN_FONT_SIZE
    return font_size


def is_valid_path(path):
    """
    Check if the given path is valid for opening a file.
    """
    # Check if the path exists and is accessible
    if os.path.exists(path) and os.access(path, os.R_OK):
        return True
    return False
