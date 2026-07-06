from typing import Dict
def extract_condition(meta_file:str) -> Dict:
    """
    Reads a metadata file and returns a dictionary mapping sample IDs to their corresponding conditions.

    Args:
        meta_file (str): Path to the metadata file.

    Returns:
        dict: A dictionary mapping sample IDs to their corresponding conditions.
    """
    sample_condition_dict = {}
    with open(meta_file, "r") as f:
        header = f.readline().strip().split("\t")
        if "sample_id" not in header or "design" not in header:
            raise ValueError("Metadata file must contain a 'sample_id' and 'design' column.")
        for line in f:
            values = line.strip().split("\t")
            sample_id = values[header.index("sample_id")]
            condition = values[header.index("design")]
            sample_condition_dict[sample_id] = condition
    return sample_condition_dict