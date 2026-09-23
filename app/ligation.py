"""跨记录片段末端连接（ligation）兼容性分析（只读）。

分析口径
--------
* 每条已保存记录取其完整酶切结果的全部片段，每个片段的左/右末端
  （``left_end`` / ``right_end``）各作为一个候选末端；
* 末端对只枚举**跨记录**组合：同一记录内的末端互不配对，同一末端
  更不会与自身配对（请求中的重复 id 在模型校验阶段即被 422 拒绝）；
* 兼容判定：
  - 黏端（sticky）：一端为 5' 突出端、另一端为 3' 突出端，且 5' 端
    序列**恰好等于** 3' 端序列的反向互补（IUPAC 逐字符反向互补后的
    字符串相等，不做简并匹配）；
  - 平端（blunt）：仅与平端兼容，且仅当请求 ``include_blunt`` 为真时
    才纳入结果；
  - 其余组合（含 ``none`` 天然末端、5'-5'、3'-3'、黏端-平端）均不兼容；
* 连接后突出端 ``ligated_overhang``：黏端对取该对中 5' 突出端的序列
  （即接头处退火的突出端），平端对为空字符串；
* 全部末端先按（记录请求顺序、片段索引、端方向 left<right）排成全序，
  再按下标 i < j 枚举配对：每对的 first 恒为全序中靠前的一端（必然
  属于请求顺序靠前的记录），返回顺序即该全序，相同请求的结果必然一致；
* 整个分析只在内存中进行，不读写数据库。
"""

from __future__ import annotations

from .iupac import reverse_complement

STICKY = "sticky"
BLUNT = "blunt"

_END_DIRECTIONS = ("left", "right")


def _fragment_ends(record_position: int, digest_id: int, fragments) -> list[dict]:
    """把一条记录的全部片段末端展开为带定位信息的末端列表。

    返回顺序即（片段索引、端方向）升序；``record_position`` 是记录在
    请求 digest_ids 中的下标，仅用于跨记录判定与排序，不对外返回。
    """
    ends = []
    for fragment in fragments:
        for direction in _END_DIRECTIONS:
            end = fragment[f"{direction}_end"]
            ends.append(
                {
                    "record_position": record_position,
                    "digest_id": digest_id,
                    "fragment_index": fragment["index"],
                    "end": direction,
                    "end_type": end["type"],
                    "overhang": end["overhang"],
                }
            )
    return ends


def match_ends(end_a: dict, end_b: dict, include_blunt: bool):
    """判定两个末端是否可连接，返回 ``(兼容类型, 连接后突出端)`` 或 None。

    末端字典至少含 ``end_type`` / ``overhang``（与片段末端结构一致）。
    本函数是末端配对枚举（:func:`find_compatible_pairs`）与片段拼接链
    规划（``app.assembly``）共用的唯一兼容判定入口，两处口径完全一致：
    黏端要求一方 5' 突出端序列等于另一方 3' 突出端序列的反向互补，
    平端仅在 ``include_blunt`` 为真时与平端兼容。判定对入参对称
    （``match_ends(a, b)`` 与 ``match_ends(b, a)`` 结论一致）。
    """
    type_a = end_a["end_type"]
    type_b = end_b["end_type"]
    if type_a == BLUNT and type_b == BLUNT:
        return (BLUNT, "") if include_blunt else None
    if type_a == "5prime" and type_b == "3prime":
        if end_a["overhang"] == reverse_complement(end_b["overhang"]):
            return (STICKY, end_a["overhang"])
        return None
    if type_a == "3prime" and type_b == "5prime":
        if end_b["overhang"] == reverse_complement(end_a["overhang"]):
            return (STICKY, end_b["overhang"])
        return None
    return None


def find_compatible_pairs(records, include_blunt: bool) -> list[dict]:
    """枚举跨记录的可连接末端对，调用方须已完成全部校验。

    ``records`` 为 ``{"digest_id": int, "fragments": [...]}`` 序列，
    顺序即请求中 digest_ids 的顺序。返回按（记录顺序、片段索引、
    端方向）稳定排列的配对列表，每项：
    ``{"first": end, "second": end, "compatibility": "sticky"|"blunt",
    "ligated_overhang": str}``；first/second 为末端的对外定位信息。
    """
    # 全部末端的全序列表：记录请求顺序 -> 片段索引 -> 端方向（left < right）
    all_ends: list[dict] = []
    for position, record in enumerate(records):
        all_ends.extend(
            _fragment_ends(position, record["digest_id"], record["fragments"])
        )

    pairs = []
    for i, end_a in enumerate(all_ends):
        for end_b in all_ends[i + 1 :]:
            # 只枚举跨记录组合；i < j 保证同一末端不会与自身配对
            if end_b["record_position"] == end_a["record_position"]:
                continue
            match = match_ends(end_a, end_b, include_blunt)
            if match is None:
                continue
            compatibility, ligated_overhang = match
            pairs.append(
                {
                    "first": _public_end(end_a),
                    "second": _public_end(end_b),
                    "compatibility": compatibility,
                    "ligated_overhang": ligated_overhang,
                }
            )
    return pairs


def _public_end(end: dict) -> dict:
    """末端的对外表示：去掉内部排序用的 record_position。"""
    return {
        "digest_id": end["digest_id"],
        "fragment_index": end["fragment_index"],
        "end": end["end"],
        "end_type": end["end_type"],
        "overhang": end["overhang"],
    }
