from blazeqc.helpers import _seq_to_hash
from collections import Dict


# TODO: Move those to a config file
##############
fn hash_list() -> List[UInt64]:
    li = List[UInt64](capacity=6)
    li.append(_seq_to_hash("AGATCGGAAGAG"))
    li.append(_seq_to_hash("TGGAATTCTCGG"))
    li.append(_seq_to_hash("GATCGTCGGACT"))
    li.append(_seq_to_hash("CTGTCTCTTATA"))
    li.append(_seq_to_hash("AAAAAAAAAAAA"))
    li.append(_seq_to_hash("GGGGGGGGGGGG"))
    return li^


# TODO: Check how to unpack this variadic
def hash_names() -> List[StringSlice[StaticConstantOrigin]]:
    var names = [
        "Illumina Universal Adapter",
        "Illumina Small RNA 3' Adapter",
        "Illumina Small RNA 5' Adapter",
        "Nextera Transposase Sequence",
        "PolyA",
        "PolyG",
    ]

    return names^


def get_hashes() -> (
    Dict[StringSlice[StaticConstantOrigin], StringSlice[StaticConstantOrigin]]
):
    var hashes = Dict[
        StringSlice[StaticConstantOrigin], StringSlice[StaticConstantOrigin]
    ]()
    hashes["Illumina Universal Adapter"] = "AGATCGGAAGAG"
    hashes["Illumina Small RNA 3' Adapter"] = "TGGAATTCTCGG"
    hashes["Illumina Small RNA 5' Adapter"] = "GATCGTCGGACT"
    hashes["Nextera Transposase Sequence"] = "CTGTCTCTTATA"
    hashes["PolyA"] = "AAAAAAAAAAAA"
    hashes["PolyG"] = "GGGGGGGGGGGG"
    return hashes^
