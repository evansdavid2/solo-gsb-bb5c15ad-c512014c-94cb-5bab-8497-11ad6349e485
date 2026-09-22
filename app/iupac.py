"""IUPAC 简并碱基编码表与匹配工具。"""

# 每个 IUPAC 编码对应可能代表的碱基集合（DNA）
IUPAC_SETS = {
    "A": frozenset("A"),
    "C": frozenset("C"),
    "G": frozenset("G"),
    "T": frozenset("T"),
    "R": frozenset("AG"),   # puRine
    "Y": frozenset("CT"),   # pYrimidine
    "S": frozenset("GC"),   # Strong
    "W": frozenset("AT"),   # Weak
    "K": frozenset("GT"),   # Keto
    "M": frozenset("AC"),   # aMino
    "B": frozenset("CGT"),  # not A
    "D": frozenset("AGT"),  # not C
    "H": frozenset("ACT"),  # not G
    "V": frozenset("ACG"),  # not T
    "N": frozenset("ACGT"), # aNy
}

# 互补映射（用于求反向互补序列）
COMPLEMENT = {
    "A": "T", "T": "A", "C": "G", "G": "C",
    "R": "Y", "Y": "R", "S": "S", "W": "W",
    "K": "M", "M": "K",
    "B": "V", "V": "B", "D": "H", "H": "D",
    "N": "N",
}

VALID_BASES = frozenset(IUPAC_SETS)


def is_valid_sequence(seq: str) -> bool:
    """序列是否仅由 IUPAC 碱基组成。"""
    return bool(seq) and all(c in VALID_BASES for c in seq)


def reverse_complement(seq: str) -> str:
    """返回 IUPAC 序列的反向互补序列。"""
    return "".join(COMPLEMENT[c] for c in reversed(seq))


def bases_match(pattern_base: str, seq_base: str) -> bool:
    """简并匹配：两者可能代表的碱基集合有交集即视为可匹配。

    例如模式中的 N 可匹配序列中的任意碱基；序列中的 R(A/G)
    可被模式中的 A 匹配（该位置"可能"被酶识别）。
    """
    return bool(IUPAC_SETS[pattern_base] & IUPAC_SETS[seq_base])
