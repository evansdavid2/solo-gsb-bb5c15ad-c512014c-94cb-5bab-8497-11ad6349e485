"""跨记录片段拼接链规划（ligation assembly chain planning，只读）。

规划口径
--------
* 每条已保存记录贡献其完整酶切结果中的**一个**片段参与拼接（一条记录
  在链中最多出现一次）；片段是无方向元件，可按原始方向（``left→right``）
  或翻转方向（``right→left``）摆入链中；
* **有向图**：节点 = （记录请求位置, 片段索引, 摆放方向），共 2 × 片段数个；
  每个节点有一个朝链右的"出端"（正向取片段右端、反向取片段左端）与一个
  朝链左的"入端"（正向取片段左端、反向取片段右端）；边 u→v 存在当且仅当
  u、v 来自**不同记录**且 u 的出端与 v 的入端通过
  :func:`app.ligation.match_ends` 判定兼容（黏端反向互补 / 平端，
  ``include_blunt`` 可关闭平端）；
* 在图上搜索**节点数最多的简单路径**（记录不重复使用——同记录节点间
  本就不设边，搜索中再以未用记录集合保证不重复），节点数受请求
  ``max_chain_length`` 限制；
* **同长度稳定择优**：长度相同的方案按节点序列
  （记录请求位置, 片段索引, 端方向 left<right——即链左侧所用的原片段端）
  逐键字典序取最小；与边的枚举顺序无关，相同请求结果必然一致；
* 搜索以记忆化（状态 = 当前出端 × 未用记录集合 × 剩余槽位）完成，并在
  同一记录的候选节点内带上界剪枝：该记录首个达到理论最长后缀的候选即其
  本记录最优；不同记录之间必须完整比较以保证跨记录字典序择优。枚举规模
  有界：8 条记录、总片段数不超过 :data:`MAX_ASSEMBLY_FRAGMENTS`；
* 整个规划只在内存中进行，不读写数据库。
"""

from __future__ import annotations

from functools import lru_cache

from .ligation import match_ends

# 参与规划的全部记录片段总数上限（每个可连接片段展开为 fwd/rev 两个图节点）：
# 记忆化状态数 = 末端数 × 2^(记录数-1) × 链长槽位，超出该上限无法在保证
# 耗时的前提下完成全枚举，路由以 422 拒绝。
MAX_ASSEMBLY_FRAGMENTS = 64

# 摆放方向：链左端所用的原片段端（left < right 即稳定择优的端方向次序）
_ORIENTATIONS = (
    # orientation, 入端(链左)取原片段, 出端(链右)取原片段, 端方向次序
    ("forward", "left", "right", 0),
    ("reverse", "right", "left", 1),
)

# 未能延伸 / 无方案的原因代码
REASON_MAX_LENGTH = "max_chain_length_reached"
REASON_ALL_USED = "all_records_used"
REASON_NO_EXTENSION = "no_compatible_extension"
REASON_NO_ENDS = "no_compatible_ends"

_STOP_DETAILS = {
    REASON_MAX_LENGTH: (
        "chain reached max_chain_length; compatible ends on unused records "
        "remain but extending further is disallowed by the limit"
    ),
    REASON_ALL_USED: (
        "every submitted record is used in the chain; there is no remaining "
        "record to extend with"
    ),
    REASON_NO_EXTENSION: (
        "the rightmost chain end is not compatible with any candidate end of "
        "the unused records"
    ),
    REASON_NO_ENDS: (
        "no pair of compatible fragment ends exists across the submitted "
        "records (a linear chain needs at least one ligation between two "
        "fragments); no chain can be assembled"
    ),
}


def _candidate_nodes(record_position: int, digest_id: int, fragments) -> list[dict]:
    """把一条记录的片段展开为有向图节点。

    两端均为天然末端（``none``）的片段不可能参与任何连接，直接排除；
    其余片段展开为正向/反向两个节点，按择优键
    （片段索引, 端方向 left<right）排序返回。
    """
    nodes = []
    for fragment in fragments:
        left = fragment["left_end"]
        right = fragment["right_end"]
        if left["type"] == "none" and right["type"] == "none":
            continue
        for orientation, in_end_name, out_end_name, direction_rank in _ORIENTATIONS:
            in_end = fragment[f"{in_end_name}_end"]
            out_end = fragment[f"{out_end_name}_end"]
            nodes.append(
                {
                    "record_position": record_position,
                    "digest_id": digest_id,
                    "fragment_index": fragment["index"],
                    "orientation": orientation,
                    "direction_rank": direction_rank,
                    "in_end": {
                        "fragment_end": in_end_name,
                        "end_type": in_end["type"],
                        "overhang": in_end["overhang"],
                    },
                    "out_end": {
                        "fragment_end": out_end_name,
                        "end_type": out_end["type"],
                        "overhang": out_end["overhang"],
                    },
                }
            )
    nodes.sort(
        key=lambda n: (n["record_position"], n["fragment_index"], n["direction_rank"])
    )
    return nodes


def _node_key(node: dict) -> tuple:
    """节点的稳定择优键：记录请求位置、片段索引、链左端方向（left<right）。"""
    return (node["record_position"], node["fragment_index"], node["direction_rank"])


def _end_key(end: dict, node: dict) -> tuple:
    """末端身份：（记录位置, 片段索引, 原片段端方向 left<right）。"""
    return (
        node["record_position"],
        node["fragment_index"],
        0 if end["fragment_end"] == "left" else 1,
    )


def plan_assembly_chain(records, include_blunt: bool, max_chain_length: int) -> dict:
    """在跨记录有向图上搜索最长线性拼接链，调用方须已完成全部校验。

    ``records`` 为 ``{"digest_id": int, "fragments": [...]}`` 序列，顺序即
    请求中 digest_ids 的顺序。返回：
    ``{"nodes": [节点...], "connections": [连接...], "stop_reason": str,
    "stop_reason_detail": str, "used_positions": frozenset}``；不存在任何
    可连接的片段对（连接数为 0）时 ``nodes`` / ``connections`` 为空，
    ``stop_reason`` 为 ``no_compatible_ends``。
    """
    record_count = len(records)

    # 全部图节点（按记录请求顺序分组，组内已按择优键排序）
    nodes_by_record: list[list[dict]] = []
    all_nodes: list[dict] = []
    for position, record in enumerate(records):
        nodes = _candidate_nodes(position, record["digest_id"], record["fragments"])
        nodes_by_record.append(nodes)
        all_nodes.extend(nodes)

    all_positions = frozenset(range(record_count))

    # 末端身份 -> 末端内容；记忆化状态只传身份键，避免字典不可哈希
    ends: dict[tuple, dict] = {}
    for node in all_nodes:
        ends[_end_key(node["in_end"], node)] = node["in_end"]
        ends[_end_key(node["out_end"], node)] = node["out_end"]

    def compatible(out_key: tuple, candidate_in_end: dict) -> bool:
        return match_ends(ends[out_key], candidate_in_end, include_blunt) is not None

    @lru_cache(maxsize=None)
    def best_suffix(out_key: tuple, unused_positions: frozenset, slots: int) -> tuple:
        """从当前出端出发可接出的最优后缀，返回节点择优键元组（不含当前节点）。

        先最大化后缀节点数；等长时取字典序最小的节点键序列。``slots`` 为
        还允许接入的节点数（链长上限减去已用节点数）。后缀长度的理论上界
        为 min(slots, 未用记录数)。剪枝只在**同一记录内部**进行：该记录
        节点按择优键排列，第一个达到上界的候选即本记录最优；不同记录之间
        必须完整比较，不能因已达长度上界而跳过请求顺序更靠后的记录
        （其首节点键可能更小，影响等长择优）。
        """
        if slots == 0:
            return ()
        best: tuple = ()
        upper_bound = min(slots, len(unused_positions))
        for position in sorted(unused_positions):  # 记录请求顺序枚举
            best_in_record: tuple = ()
            for node in nodes_by_record[position]:
                if not compatible(out_key, node["in_end"]):
                    continue
                key = _node_key(node)
                suffix = best_suffix(
                    _end_key(node["out_end"], node),
                    unused_positions - {position},
                    slots - 1,
                )
                candidate = (key,) + suffix
                if len(candidate) > len(best_in_record) or (
                    len(candidate) == len(best_in_record)
                    and candidate < best_in_record
                ):
                    best_in_record = candidate
                if len(best_in_record) == upper_bound:
                    break  # 同记录内已找到理论最长后缀
            if len(best_in_record) > len(best) or (
                len(best_in_record) == len(best) and best_in_record < best
            ):
                best = best_in_record
        return best

    # 枚举起点（任意记录的任意节点均可作链首），按同一择优口径选全局最优。
    # 链长上界 = min(max_chain_length, 记录数)。同枚举后缀：必须比较全部
    # 记录的全部候选起点，不能在达到长度上界后跳过请求顺序更靠后的记录
    # （其起点节点键可能更小，影响等长择优）。
    best_path: tuple = ()
    for position in range(record_count):
        for node in nodes_by_record[position]:
            suffix = best_suffix(
                _end_key(node["out_end"], node),
                all_positions - {position},
                max_chain_length - 1,
            )
            candidate = (_node_key(node),) + suffix
            if len(candidate) > len(best_path) or (
                len(candidate) == len(best_path) and candidate < best_path
            ):
                best_path = candidate

    # 拼接链至少要有一处连接（两个节点）；不足即视为不存在可行方案
    if len(best_path) < 2:
        return {
            "nodes": [],
            "connections": [],
            "stop_reason": REASON_NO_ENDS,
            "stop_reason_detail": _STOP_DETAILS[REASON_NO_ENDS],
            "used_positions": frozenset(),
        }

    # 择优键 -> 节点对象（同一键唯一标识一个图节点）
    node_by_key = {_node_key(node): node for node in all_nodes}
    chain_nodes = [node_by_key[key] for key in best_path]

    connections = []
    for index in range(len(chain_nodes) - 1):
        left_node = chain_nodes[index]
        right_node = chain_nodes[index + 1]
        match = match_ends(
            left_node["out_end"], right_node["in_end"], include_blunt
        )
        # best_suffix 只放行兼容边；不匹配属于内部错误，尽早暴露
        assert match is not None
        compatibility, ligated_overhang = match
        connections.append(
            {
                "index": index,
                "from_end": _public_end(left_node, left_node["out_end"]),
                "to_end": _public_end(right_node, right_node["in_end"]),
                "compatibility": compatibility,
                "ligated_overhang": ligated_overhang,
            }
        )

    used_positions = frozenset(n["record_position"] for n in chain_nodes)
    stop_reason = _stop_reason(
        chain_nodes, nodes_by_record, used_positions, all_positions,
        max_chain_length, include_blunt,
    )
    return {
        "nodes": chain_nodes,
        "connections": connections,
        "stop_reason": stop_reason,
        "stop_reason_detail": _STOP_DETAILS[stop_reason],
        "used_positions": used_positions,
    }


def _stop_reason(
    chain_nodes, nodes_by_record, used_positions, all_positions,
    max_chain_length, include_blunt,
) -> str:
    """判定最优链无法继续延伸的原因（仅用于至少含一处连接的方案）。"""
    unused_positions = all_positions - used_positions
    if not unused_positions:
        return REASON_ALL_USED

    # 仍有未用记录：检查链右端是否还能与其中任一候选入端兼容
    out_end = chain_nodes[-1]["out_end"]
    extendable = any(
        match_ends(out_end, node["in_end"], include_blunt) is not None
        for position in unused_positions
        for node in nodes_by_record[position]
    )
    if extendable and len(chain_nodes) >= max_chain_length:
        return REASON_MAX_LENGTH
    return REASON_NO_EXTENSION


def _public_end(node: dict, end: dict) -> dict:
    """节点端点的对外表示：记录 id、片段索引、原片段端方向与末端类型/突出端。"""
    return {
        "digest_id": node["digest_id"],
        "fragment_index": node["fragment_index"],
        "end": end["fragment_end"],
        "end_type": end["end_type"],
        "overhang": end["overhang"],
    }


def public_node(node: dict) -> dict:
    """图节点的对外表示：链左/右端点均带其在原片段上的端方向。"""
    return {
        "digest_id": node["digest_id"],
        "fragment_index": node["fragment_index"],
        "orientation": node["orientation"],
        "left_end": _public_end(node, node["in_end"]),
        "right_end": _public_end(node, node["out_end"]),
    }
