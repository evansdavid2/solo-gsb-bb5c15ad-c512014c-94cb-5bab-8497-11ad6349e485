"""甲基化阻断分析：只读评估甲基化对已保存酶切结果的影响。

分析口径
--------
* 请求为每种酶给出识别模式内的**敏感偏移**集合（0-based，相对于模式 5' 端）。
* 对已保存切口上的每个识别位点，把敏感偏移换算为 top 链物理坐标：
  - 正向位点（strand == "+"）：``start + offset``；
  - 反向位点（strand == "-"）：``start + k - 1 - offset``（k 为模式长度）；
  - 环形序列对序列长度取模，使跨原点位点的偏移落到真实碱基坐标。
* 位点的任一敏感物理坐标出现在 ``methylated_positions`` 中，即视为该位点
  **被甲基化阻断**。
* 按物理切口 ``(top_cut, bottom_cut)`` 汇总：
  - 仍存在未阻断位点 -> 保留该切口，并重算 ``enzymes`` / ``sites``
    （酶名去重、位点按酶名/起点/链排序，与完整酶切口径一致）；
  - 全部位点均被阻断 -> 移除该切口。
* 片段由保留下来的切口复用 ``digest.build_fragments`` 重建，线性/环形的
  全部边界规则与完整酶切一致。

整个分析只在内存中进行，不读写数据库，原记录保持不变。
"""

from __future__ import annotations

from .digest import build_fragments


def _sensitive_positions(site: dict, offsets, circular: bool, length: int) -> set[int]:
    """把一个位点上的敏感偏移集合换算为 top 链物理坐标集合。"""
    start = site["start"]
    k = len(site["matched"])
    if site["strand"] == "+":
        raw = (start + offset for offset in offsets)
    else:
        raw = (start + k - 1 - offset for offset in offsets)
    if circular:
        return {position % length for position in raw}
    return set(raw)


def analyze_methylation(
    sequence: str,
    topology: str,
    cuts: list[dict],
    rules,
    methylated_positions,
) -> dict:
    """评估甲基化对已保存切口的影响，调用方须已完成全部校验。

    ``rules`` 为 ``(酶名, 敏感偏移可迭代对象)`` 序列；酶名允许重复出现
    （调用方的请求模型保证不重复）。``methylated_positions`` 为 top 链
    0-based 物理坐标集合。

    返回 {"blocked_sites": [...], "retained_cuts": [...],
    "removed_cuts": [...], "fragments": [...]}，各列表均稳定排序。
    """
    circular = topology == "circular"
    length = len(sequence)
    methylated = set(methylated_positions)
    rules_by_enzyme: dict[str, set[int]] = {}
    for enzyme, offsets in rules:
        rules_by_enzyme.setdefault(enzyme, set()).update(offsets)

    blocked_sites: list[dict] = []
    retained_cuts: list[dict] = []
    removed_cuts: list[dict] = []

    for cut in cuts:
        surviving_sites: list[dict] = []
        for site in cut["sites"]:
            offsets = rules_by_enzyme.get(site["enzyme"], set())
            hits = _sensitive_positions(site, offsets, circular, length) & methylated
            if hits:
                blocked_sites.append(
                    {
                        "top_cut": cut["top_cut"],
                        "bottom_cut": cut["bottom_cut"],
                        "site": site,
                        "methylated_positions": sorted(hits),
                    }
                )
            else:
                surviving_sites.append(site)

        if surviving_sites:
            retained_cuts.append(
                {
                    "top_cut": cut["top_cut"],
                    "bottom_cut": cut["bottom_cut"],
                    "end_type": cut["end_type"],
                    "overhang": cut["overhang"],
                    "enzymes": sorted({site["enzyme"] for site in surviving_sites}),
                    "sites": sorted(
                        surviving_sites,
                        key=lambda s: (s["enzyme"], s["start"], s["strand"]),
                    ),
                }
            )
        else:
            # 切口整体消失：回显其被阻断前的完整对象（酶名与位点均保留）
            removed_cuts.append(
                {
                    "top_cut": cut["top_cut"],
                    "bottom_cut": cut["bottom_cut"],
                    "end_type": cut["end_type"],
                    "overhang": cut["overhang"],
                    "enzymes": list(cut["enzymes"]),
                    "sites": list(cut["sites"]),
                }
            )

    retained_cuts.sort(key=lambda c: (c["top_cut"], c["bottom_cut"]))
    removed_cuts.sort(key=lambda c: (c["top_cut"], c["bottom_cut"]))
    blocked_sites.sort(
        key=lambda b: (
            b["top_cut"],
            b["bottom_cut"],
            b["site"]["enzyme"],
            b["site"]["start"],
            b["site"]["strand"],
        )
    )

    fragments = build_fragments(sequence, topology, retained_cuts)

    return {
        "blocked_sites": blocked_sites,
        "retained_cuts": retained_cuts,
        "removed_cuts": removed_cuts,
        "fragments": fragments,
    }
