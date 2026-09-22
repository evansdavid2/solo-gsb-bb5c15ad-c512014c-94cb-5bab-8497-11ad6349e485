"""限制性酶切模拟核心算法。

坐标约定
--------
* 所有坐标均为 0-based，作用于 top 链（5'->3'，即用户提交的序列方向）。
* "切口位置 c" 表示断裂发生在第 c-1 与第 c 个碱基之间（c ∈ [0, L]）。
* top 链切口与 bottom 链切口位置相同即为平端；不同则产生黏端。
* 环形序列的切口位置对 L 取模归一化到 [0, L)；位置 0 即跨越原点的切口。

末端类型
--------
* top_cut == bottom_cut        -> blunt（平端）
* top_cut <  bottom_cut        -> 5prime（5' 突出端），突出序列 = top 链 [top, bottom) 区段
* top_cut >  bottom_cut        -> 3prime（3' 突出端），突出序列 = top 链 [bottom, top) 区段
* 线性分子的天然末端 / 未切割的完整分子 -> none
"""

from __future__ import annotations

from .iupac import bases_match, reverse_complement

BLUNT = "blunt"
FIVE_PRIME = "5prime"
THREE_PRIME = "3prime"
NONE = "none"


def find_sites(sequence: str, pattern: str, circular: bool):
    """在 top 链上查找 pattern 的所有匹配起点（允许重叠匹配）。

    环形序列在 seq + seq[:k-1] 的扩展序列上搜索，使跨越首尾的
    识别位点也能命中；起点范围仍为 [0, L)。识别模式长于序列时不匹配。
    """
    length = len(sequence)
    k = len(pattern)
    if k == 0 or k > length:
        return
    haystack = sequence + sequence[: k - 1] if circular else sequence
    last_start = length if circular else length - k + 1
    for i in range(last_start):
        for j, pbase in enumerate(pattern):
            if not bases_match(pbase, haystack[i + j]):
                break
        else:
            yield i


def _end_type_and_overhang(haystack: str, top_raw: int, bottom_raw: int):
    """由未归一化的切口位置计算末端类型与突出端序列。

    必须使用取模前的原始坐标：环形分子上切口跨原点时，归一化后的
    数值大小关系会被打乱，而原始坐标始终保持真实的先后关系。
    """
    if top_raw < bottom_raw:
        return FIVE_PRIME, haystack[top_raw:bottom_raw]
    if top_raw > bottom_raw:
        return THREE_PRIME, haystack[bottom_raw:top_raw]
    return BLUNT, ""


def compute_cuts(sequence: str, topology: str, enzymes) -> list[dict]:
    """计算所有酶的切割事件，并按物理切口 (top, bottom) 去重。

    每个酶同时在两个方向搜索：top 链上的 pattern（正向）与
    revcomp(pattern)（反向，即酶结合在 bottom 链上）。回文模式只搜一次。
    反向命中时切口坐标按酶自身坐标系换算回 top 链坐标：
        top_cut    = start + k - cut_bottom
        bottom_cut = start + k - cut_top
    同一物理切口（top、bottom 均相同）只保留一条记录，酶名与位点合并。
    返回列表按 (top_cut, bottom_cut) 排序，保证相同输入的输出顺序一致。
    """
    circular = topology == "circular"
    length = len(sequence)
    # 扩展序列用于提取跨原点的突出端序列
    max_k = max(len(e.pattern) for e in enzymes)
    haystack = sequence + sequence[:max_k] if circular else sequence

    merged: dict[tuple[int, int], dict] = {}
    for enz in sorted(enzymes, key=lambda e: (e.name, e.pattern, e.cut_top, e.cut_bottom)):
        k = len(enz.pattern)
        for oriented in sorted({enz.pattern, reverse_complement(enz.pattern)}):
            forward = oriented == enz.pattern
            for start in find_sites(sequence, oriented, circular):
                if forward:
                    top_raw = start + enz.cut_top
                    bottom_raw = start + enz.cut_bottom
                else:
                    top_raw = start + k - enz.cut_bottom
                    bottom_raw = start + k - enz.cut_top
                top = top_raw % length if circular else top_raw
                bottom = bottom_raw % length if circular else bottom_raw
                end_type, overhang = _end_type_and_overhang(haystack, top_raw, bottom_raw)

                key = (top, bottom)
                entry = merged.get(key)
                if entry is None:
                    entry = merged[key] = {
                        "top_cut": top,
                        "bottom_cut": bottom,
                        "end_type": end_type,
                        "overhang": overhang,
                        "enzymes": set(),
                        "sites": [],
                    }
                entry["enzymes"].add(enz.name)
                entry["sites"].append(
                    {
                        "enzyme": enz.name,
                        "start": start,
                        "strand": "+" if forward else "-",
                        "matched": haystack[start : start + k],
                    }
                )

    cuts = []
    for entry in merged.values():
        entry["enzymes"] = sorted(entry["enzymes"])
        entry["sites"].sort(key=lambda s: (s["enzyme"], s["start"], s["strand"]))
        cuts.append(entry)
    cuts.sort(key=lambda c: (c["top_cut"], c["bottom_cut"]))
    return cuts


def _end_dict(cut: dict | None) -> dict:
    if cut is None:
        return {"type": NONE, "overhang": ""}
    return {"type": cut["end_type"], "overhang": cut["overhang"]}


def build_fragments(sequence: str, topology: str, cuts: list[dict]) -> list[dict]:
    """按拓扑结构由切口生成片段。

    * 线性：n 个有效切口产生 n+1 个片段；位于分子端点 (0 或 L) 的
      切口不产生新边界。
    * 环形：n 个切口产生 n 个片段（仅 1 个切口时开环为 1 条全长片段）；
      无切口时保留整条环状分子。
    """
    length = len(sequence)
    fragments = []

    if topology == "circular":
        if not cuts:
            fragments.append(
                {
                    "start": 0,
                    "end": 0,
                    "length": length,
                    "sequence": sequence,
                    "wraps_origin": False,
                    "left_end": _end_dict(None),
                    "right_end": _end_dict(None),
                }
            )
        else:
            n = len(cuts)
            for idx, cut in enumerate(cuts):
                nxt = cuts[(idx + 1) % n]
                start = cut["top_cut"]
                if n == 1:
                    frag_len = length  # 单切口开环，得到整条分子的线性化序列
                else:
                    frag_len = (nxt["top_cut"] - start) % length
                    if frag_len == 0:
                        continue  # 同一切点的退化零长片段
                end = (start + frag_len) % length
                if start + frag_len >= length:
                    frag_seq = sequence[start:] + sequence[: (start + frag_len) % length]
                else:
                    frag_seq = sequence[start : start + frag_len]
                fragments.append(
                    {
                        "start": start,
                        "end": end,
                        "length": frag_len,
                        "sequence": frag_seq,
                        "wraps_origin": start + frag_len > length,
                        "left_end": _end_dict(cut),
                        "right_end": _end_dict(nxt),
                    }
                )
    else:
        events_by_pos: dict[int, dict] = {}
        for cut in cuts:  # cuts 已按 (top, bottom) 排序，同位置取第一条，保证确定性
            pos = cut["top_cut"]
            if 0 < pos < length and pos not in events_by_pos:
                events_by_pos[pos] = cut
        points = [0] + sorted(events_by_pos) + [length]
        for i in range(len(points) - 1):
            s, e = points[i], points[i + 1]
            fragments.append(
                {
                    "start": s,
                    "end": e,
                    "length": e - s,
                    "sequence": sequence[s:e],
                    "wraps_origin": False,
                    "left_end": _end_dict(events_by_pos.get(s)),
                    "right_end": _end_dict(events_by_pos.get(e)),
                }
            )

    return [
        {"index": i, **frag} for i, frag in enumerate(fragments)
    ]


def run_digest(sequence: str, topology: str, enzymes) -> dict:
    """执行一次完整酶切，返回 {"cuts": [...], "fragments": [...]}。"""
    cuts = compute_cuts(sequence, topology, enzymes)
    fragments = build_fragments(sequence, topology, cuts)
    return {"cuts": cuts, "fragments": fragments}
