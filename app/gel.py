"""凝胶条带比较：多条已保存酶切记录的片段长度聚类（只读）。

比较口径
--------
* 每条已保存记录视为一个**泳道**，取其完整酶切结果的全部片段长度；
* 所有泳道的片段长度汇总后按**降序**贪心聚类为条带：
  - 每个条带以其最大长度（簇内首个元素）为锚点，新长度 x 并入当前
    条带当且仅当 ``(max - x) / max <= resolution_ratio``；
  - 判定始终相对条带的**最大值**，而非相邻的上一个元素，因此不会
    出现"相邻两两可合并、首尾却超差"的链式归并；
  - 长度降序遍历时当前条带是已有条带中 max 最小者，元素若进不了
    当前条带，也不可能进任何更早的条带，贪心结果即确定；
* 各条带的 max 必然互不相同（相同长度必然并入同一条带），条带按
  max_length 降序返回的次序因此完全确定；
* 仅出现在一个泳道的条带标记为 diagnostic（诊断条带）；
* 整个比较只在内存中进行，不读写数据库。
"""

from __future__ import annotations


def compare_gel_lanes(lanes: list[dict], resolution_ratio: float) -> list[dict]:
    """对若干泳道的片段长度聚类，返回按 max_length 降序的条带列表。

    ``lanes`` 为 ``{"digest_id": int, "fragment_lengths": [int, ...]}``
    序列，泳道顺序即请求中 digest_ids 的顺序；返回的
    ``lane_fragment_counts`` 与该顺序一一对应。调用方须已完成全部校验
    （泳道数 2~8、resolution_ratio ∈ [0, 0.2]、记录均存在）。

    每个条带：{"min_length", "max_length", "lane_fragment_counts",
    "fragment_count", "total_bases", "diagnostic"}。
    """
    lane_count = len(lanes)
    # 汇总全部泳道的片段长度并降序排列；酶切片段长度恒 >= 1，
    # 因此条带 max 不会为 0，聚类判定不存在除零
    pooled: list[tuple[int, int]] = [
        (length, lane_index)
        for lane_index, lane in enumerate(lanes)
        for length in lane["fragment_lengths"]
    ]
    pooled.sort(key=lambda item: -item[0])

    bands: list[dict] = []
    for length, lane_index in pooled:
        if not bands or _spread(bands[-1]["max_length"], length) > resolution_ratio:
            bands.append(
                {
                    "min_length": length,
                    "max_length": length,
                    "lane_fragment_counts": [0] * lane_count,
                    "fragment_count": 0,
                    "total_bases": 0,
                }
            )
        band = bands[-1]
        band["min_length"] = length  # 降序遍历，最新并入的元素即簇内最小值
        band["lane_fragment_counts"][lane_index] += 1
        band["fragment_count"] += 1
        band["total_bases"] += length

    for band in bands:
        occupied = sum(1 for count in band["lane_fragment_counts"] if count > 0)
        band["diagnostic"] = occupied == 1
    return bands


def _spread(max_length: int, length: int) -> float:
    """长度并入条带后的相对跨度 (max - length) / max。"""
    return (max_length - length) / max_length
