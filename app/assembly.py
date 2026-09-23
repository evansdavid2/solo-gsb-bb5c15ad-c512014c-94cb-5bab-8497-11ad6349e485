"""片段拼接链规划（只读）。

分析口径
--------
* 每条已保存记录取其完整酶切结果的全部片段，每个片段的左/右末端
  （``left_end`` / ``right_end``）各作为一个候选末端；
* 以末端为顶点、兼容末端对为边构建有向图：兼容判定复用
  ``ligation.match_ends``——黏端要求一方 5' 突出端序列等于另一方 3'
  突出端序列的反向互补，平端仅与平端兼容（``include_blunt`` 可关闭）；
* 链中每个节点使用一条记录的一个片段，**每条记录最多出现一次**
  （请求中的重复 id 在模型校验阶段即被 422 拒绝）；相邻节点必须使用
  一对兼容末端连接；
* 链为线性：中间节点以 left/right 两端分别连接前驱与后继，首/尾节点
  各只使用一个端，未参与连接的自由端可以是任意类型（含 ``none``
  天然末端）；
* 链长 = 节点数，上限为 ``max_chain_length``；先以（自由端签名,
  已用记录掩码）为状态做记忆化搜索求出可达的最长链长，再按
  （记录请求顺序、片段索引、端方向）逐节点字典序贪心选出键最小的
  方案——首节点比较其连接后继的端方向，其余节点比较其连接前驱的
  端方向（left < right），相同请求必然得到相同方案；
* 整个规划只在内存中进行，不读写数据库。
"""

from __future__ import annotations

from functools import lru_cache

from .ligation import match_ends

ALL_RECORDS_USED = "all_records_used"
MAX_CHAIN_LENGTH_REACHED = "max_chain_length_reached"
NO_COMPATIBLE_END = "no_compatible_end"

_DIRECTIONS = ("left", "right")
_OTHER_DIRECTION = {"left": "right", "right": "left"}


def _signature(end: dict) -> tuple:
    """末端的兼容签名：连接判定只取决于末端类型与突出端序列。"""
    return (end["type"], end["overhang"])


def plan_assembly_chain(records, include_blunt: bool, max_chain_length: int) -> dict:
    """搜索不重复使用记录的最长线性拼接链，调用方须已完成全部校验。

    ``records`` 为 ``{"digest_id": int, "fragments": [...]}`` 序列，
    顺序即请求中 digest_ids 的顺序。返回：
    ``{"chain_length", "junction_count", "stop_reason", "nodes", "junctions"}``。
    """
    # 片段按（记录请求顺序、片段索引）全序展开；末端签名 = (end_type, overhang)
    fragments = []  # (rec_pos, frag_idx, digest_id, {"left": sig, "right": sig})
    for rec_pos, record in enumerate(records):
        for fragment in record["fragments"]:
            fragments.append(
                (
                    rec_pos,
                    fragment["index"],
                    record["digest_id"],
                    {
                        direction: _signature(fragment[f"{direction}_end"])
                        for direction in _DIRECTIONS
                    },
                )
            )

    # 候选末端全序列表：（记录请求顺序、片段索引、端方向 left < right）
    candidate_ends = []  # (rec_pos, frag_idx, direction, sig)
    frag_ends = {}  # (rec_pos, frag_idx) -> {"left": sig, "right": sig}
    digest_id_of = {}  # (rec_pos, frag_idx) -> digest_id
    for rec_pos, frag_idx, digest_id, ends in fragments:
        frag_ends[(rec_pos, frag_idx)] = ends
        digest_id_of[(rec_pos, frag_idx)] = digest_id
        for direction in _DIRECTIONS:
            candidate_ends.append((rec_pos, frag_idx, direction, ends[direction]))

    def _compatible(sig_a, sig_b) -> bool:
        return (
            match_ends(
                {"end_type": sig_a[0], "overhang": sig_a[1]},
                {"end_type": sig_b[0], "overhang": sig_b[1]},
                include_blunt,
            )
            is not None
        )

    @lru_cache(maxsize=None)
    def max_remaining(end_type, overhang, mask):
        """当前自由端签名 + 已用记录掩码下，还能接入链尾的最大节点数。"""
        best = 0
        for rec_pos, frag_idx, direction, sig in candidate_ends:
            bit = 1 << rec_pos
            if mask & bit or not _compatible((end_type, overhang), sig):
                continue
            other = frag_ends[(rec_pos, frag_idx)][_OTHER_DIRECTION[direction]]
            best = max(best, 1 + max_remaining(other[0], other[1], mask | bit))
        return best

    # 首节点选项：（记录请求顺序、片段索引、连接后继的端方向）全序
    start_options = [
        (rec_pos, frag_idx, direction)
        for rec_pos, frag_idx, _digest_id, _ends in fragments
        for direction in _DIRECTIONS
    ]

    # 最长链长：任一单节点即长度 1，故至少为 1；再受 max_chain_length 截断
    best_length = min(
        max_chain_length,
        max(
            1
            + max_remaining(
                *frag_ends[(rec_pos, frag_idx)][direction], 1 << rec_pos
            )
            for rec_pos, frag_idx, direction in start_options
        ),
    )

    # 链节点：(rec_pos, frag_idx, in_dir|None, out_dir|None)
    chain = []
    if best_length == 1:
        # 无连接可用（或被限长为 1）：取全序最小的片段，不使用任何端
        chain.append((fragments[0][0], fragments[0][1], None, None))
    else:
        # 首节点：全序中第一个能延伸达到 best_length 的选项
        for rec_pos, frag_idx, out_dir in start_options:
            sig = frag_ends[(rec_pos, frag_idx)][out_dir]
            if 1 + max_remaining(*sig, 1 << rec_pos) >= best_length:
                chain.append((rec_pos, frag_idx, None, out_dir))
                break
        mask = 1 << chain[0][0]
        free_sig = frag_ends[(chain[0][0], chain[0][1])][chain[0][3]]
        # 逐位贪心：每次取全序最小且仍可达成总长度的候选端
        for position in range(2, best_length + 1):
            for rec_pos, frag_idx, in_dir, sig in candidate_ends:
                bit = 1 << rec_pos
                if mask & bit or not _compatible(free_sig, sig):
                    continue
                new_mask = mask | bit
                other = frag_ends[(rec_pos, frag_idx)][_OTHER_DIRECTION[in_dir]]
                if max_remaining(*other, new_mask) >= best_length - position:
                    out_dir = (
                        None
                        if position == best_length
                        else _OTHER_DIRECTION[in_dir]
                    )
                    chain.append((rec_pos, frag_idx, in_dir, out_dir))
                    mask = new_mask
                    free_sig = other
                    break
    assert len(chain) == best_length

    nodes = [
        {
            "digest_id": digest_id_of[(rec_pos, frag_idx)],
            "fragment_index": frag_idx,
            "left_end": {
                "type": frag_ends[(rec_pos, frag_idx)]["left"][0],
                "overhang": frag_ends[(rec_pos, frag_idx)]["left"][1],
            },
            "right_end": {
                "type": frag_ends[(rec_pos, frag_idx)]["right"][0],
                "overhang": frag_ends[(rec_pos, frag_idx)]["right"][1],
            },
        }
        for rec_pos, frag_idx, _in_dir, _out_dir in chain
    ]

    junctions = []
    for prev, nxt in zip(chain, chain[1:]):
        from_end = _junction_end(prev[0], prev[1], prev[3], frag_ends, digest_id_of)
        to_end = _junction_end(nxt[0], nxt[1], nxt[2], frag_ends, digest_id_of)
        compatibility, ligated_overhang = match_ends(
            {"end_type": from_end["end_type"], "overhang": from_end["overhang"]},
            {"end_type": to_end["end_type"], "overhang": to_end["overhang"]},
            include_blunt,
        )
        junctions.append(
            {
                "from_end": from_end,
                "to_end": to_end,
                "compatibility": compatibility,
                "ligated_overhang": ligated_overhang,
            }
        )

    used = {rec_pos for rec_pos, _frag_idx, _i, _o in chain}
    if len(used) == len(records):
        stop_reason = ALL_RECORDS_USED
    elif best_length == max_chain_length:
        stop_reason = MAX_CHAIN_LENGTH_REACHED
    else:
        stop_reason = NO_COMPATIBLE_END

    return {
        "chain_length": len(chain),
        "junction_count": len(junctions),
        "stop_reason": stop_reason,
        "nodes": nodes,
        "junctions": junctions,
    }


def _junction_end(rec_pos, frag_idx, direction, frag_ends, digest_id_of) -> dict:
    """链中一处连接使用的末端的对外表示。"""
    sig = frag_ends[(rec_pos, frag_idx)][direction]
    return {
        "digest_id": digest_id_of[(rec_pos, frag_idx)],
        "fragment_index": frag_idx,
        "end": direction,
        "end_type": sig[0],
        "overhang": sig[1],
    }
