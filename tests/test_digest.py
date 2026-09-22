"""酶切核心算法测试。"""

from app.digest import run_digest
from app.iupac import reverse_complement


class FakeEnzyme:
    def __init__(self, name, pattern, cut_top, cut_bottom):
        self.name = name
        self.pattern = pattern
        self.cut_top = cut_top
        self.cut_bottom = cut_bottom


# 常用酶：EcoRI G^AATTC（5' 黏端 AATT）
ECORI = FakeEnzyme("EcoRI", "GAATTC", 1, 5)
# EcoRV GAT^ATC（平端）
ECORV = FakeEnzyme("EcoRV", "GATATC", 3, 3)
# KpnI GGTAC^C（3' 黏端 GTAC）
KPNI = FakeEnzyme("KpnI", "GGTACC", 5, 1)


def test_reverse_complement_iupac():
    assert reverse_complement("GAATTC") == "GAATTC"  # 回文
    assert reverse_complement("GCGGCC") == "GGCCGC"
    assert reverse_complement("ACGTYR") == "YRACGT"


def test_ecori_linear_sticky_5prime():
    seq = "AAAAGAATTCTTTT"  # 位点起于 4
    result = run_digest(seq, "linear", [ECORI])
    assert len(result["cuts"]) == 1
    cut = result["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (5, 9)
    assert cut["end_type"] == "5prime"
    assert cut["overhang"] == "AATT"
    assert cut["enzymes"] == ["EcoRI"]
    assert cut["sites"][0]["start"] == 4
    assert cut["sites"][0]["strand"] == "+"

    frags = result["fragments"]
    assert len(frags) == 2
    assert (frags[0]["start"], frags[0]["end"], frags[0]["length"]) == (0, 5, 5)
    assert frags[0]["sequence"] == "AAAAG"
    assert frags[0]["left_end"]["type"] == "none"
    assert frags[0]["right_end"] == {"type": "5prime", "overhang": "AATT"}
    assert (frags[1]["start"], frags[1]["end"], frags[1]["length"]) == (5, 14, 9)
    assert frags[1]["sequence"] == "AATTCTTTT"
    assert frags[1]["left_end"] == {"type": "5prime", "overhang": "AATT"}
    assert frags[1]["right_end"]["type"] == "none"


def test_ecorv_blunt():
    seq = "CCGATATCGG"  # 位点起于 2
    result = run_digest(seq, "linear", [ECORV])
    cut = result["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (5, 5)
    assert cut["end_type"] == "blunt"
    assert cut["overhang"] == ""
    assert [f["length"] for f in result["fragments"]] == [5, 5]


def test_kpnI_sticky_3prime():
    seq = "AAAGGTACCTTT"  # 位点起于 3
    result = run_digest(seq, "linear", [KPNI])
    cut = result["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (8, 4)
    assert cut["end_type"] == "3prime"
    assert cut["overhang"] == "GTAC"


def test_circular_site_spanning_origin():
    # GAATTC 跨越首尾：位置 8,9,0,1,2,3
    seq = "ATTCGGGGGA"
    result = run_digest(seq, "circular", [ECORI])
    assert len(result["cuts"]) == 1
    cut = result["cuts"][0]
    assert cut["sites"][0]["start"] == 8
    assert (cut["top_cut"], cut["bottom_cut"]) == (9, 3)
    assert cut["end_type"] == "5prime"
    assert cut["overhang"] == "AATT"
    # 单切口开环：一条全长片段
    assert len(result["fragments"]) == 1
    frag = result["fragments"][0]
    assert frag["length"] == len(seq)
    assert frag["start"] == 9
    assert frag["sequence"] == seq[9:] + seq[:9]
    assert frag["left_end"]["type"] == "5prime"
    assert frag["right_end"]["type"] == "5prime"


def test_circular_two_cuts_two_fragments():
    # 两个 EcoRI 位点，均不跨原点
    seq = "GAATTC" + "CCCC" + "GAATTC" + "GGGG"  # 位点 0 与 10
    result = run_digest(seq, "circular", [ECORI])
    assert len(result["cuts"]) == 2
    frags = result["fragments"]
    assert len(frags) == 2
    assert sorted(f["length"] for f in frags) == [10, 10]
    assert all(f["left_end"]["type"] == "5prime" for f in frags)


def test_ambiguous_bases_in_pattern_and_sequence():
    # 模式含 N：ANNT 匹配 AGTT
    enz = FakeEnzyme("Amb1", "ANNT", 1, 1)
    result = run_digest("CAGTT", "linear", [enz])
    assert len(result["cuts"]) == 1
    assert result["cuts"][0]["sites"][0]["start"] == 1
    # 序列含简并碱基：模式 GCGC 可匹配 GCGN（N 与 C 有交集）
    enz2 = FakeEnzyme("Amb2", "GCGC", 2, 2)
    result2 = run_digest("GGCGCN", "linear", [enz2])
    assert len(result2["cuts"]) == 1
    assert result2["cuts"][0]["sites"][0]["start"] == 1


def test_overlapping_sites():
    # GATATC 与 TATCGA 重叠于同一区域
    enz1 = FakeEnzyme("E1", "GATATC", 3, 3)
    enz2 = FakeEnzyme("E2", "TATCGA", 2, 2)
    seq = "AAGATATCGAA"
    result = run_digest(seq, "linear", [enz1, enz2])
    assert len(result["cuts"]) == 2
    starts = sorted(s["start"] for c in result["cuts"] for s in c["sites"])
    assert starts == [2, 4]


def test_reverse_strand_site():
    # 非回文模式 GCGGCC，其反向互补 GGCCGC 出现在 top 链上
    enz = FakeEnzyme("Rev", "GCGGCC", 2, 4)
    seq = "AAGGCCGCTT"  # GGCCGC 起于 2
    result = run_digest(seq, "linear", [enz])
    assert len(result["cuts"]) == 1
    cut = result["cuts"][0]
    assert cut["sites"][0]["strand"] == "-"
    # 反向命中：top = 2 + 6 - 4 = 4，bottom = 2 + 6 - 2 = 6
    assert (cut["top_cut"], cut["bottom_cut"]) == (4, 6)
    assert cut["end_type"] == "5prime"
    assert cut["overhang"] == "CC"  # seq[4:6]


def test_physical_cut_deduplicated():
    # 两个酶产生完全相同的物理切口 -> 只保留一条切口记录
    other = FakeEnzyme("EcoRI-twin", "GAATTC", 1, 5)
    result = run_digest("AAAAGAATTCTTTT", "linear", [ECORI, other])
    assert len(result["cuts"]) == 1
    assert result["cuts"][0]["enzymes"] == ["EcoRI", "EcoRI-twin"]
    assert len(result["fragments"]) == 2


def test_no_cut_keeps_whole_molecule():
    seq = "AAAAAAAAAA"
    for topology in ("linear", "circular"):
        result = run_digest(seq, topology, [ECORI])
        assert result["cuts"] == []
        assert len(result["fragments"]) == 1
        frag = result["fragments"][0]
        assert frag["length"] == len(seq)
        assert frag["sequence"] == seq
        assert frag["left_end"]["type"] == "none"
        assert frag["right_end"]["type"] == "none"


def test_multiple_enzymes_combined():
    seq = "GAATTC" + "AAAA" + "GATATC" + "TTTT"  # EcoRI@0, EcoRV@10
    result = run_digest(seq, "linear", [ECORI, ECORV])
    assert len(result["cuts"]) == 2
    assert [f["length"] for f in result["fragments"]] == [1, 12, 7]
    # 切口按 top 位置排序
    assert [c["top_cut"] for c in result["cuts"]] == [1, 13]


def test_deterministic_output_order():
    seq = "GAATTCAAAAAGATATCTTTTGGTACCGGGG"
    enzymes = [KPNI, ECORI, ECORV]
    first = run_digest(seq, "circular", enzymes)
    second = run_digest(seq, "circular", list(reversed(enzymes)))
    assert first == second  # 与酶传入顺序无关，结果完全一致
