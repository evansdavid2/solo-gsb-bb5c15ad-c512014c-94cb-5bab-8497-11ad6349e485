"""部分酶切（partial digest）概率分布。

对一次已保存的完整酶切结果，其每个物理切口都只以给定概率独立"成功"
切割；其余情形视为该切口未切。枚举所有成功/失败组合并重建片段：

* 至多枚举 12 个切口（2^12 = 4096 种状态），保证耗时有界；
* 片段重建复用 ``digest.build_fragments``，线性/环形拓扑的全部边界
  规则（端点切口、单切口开环、同点退化零长片段等）与完整酶切一致；
* 片段长度多重集相同的状态（切割位置不同但产物长度分布相同）合并，
  概率用 ``math.fsum`` 精确累加；
* 结果按概率降序稳定返回，概率相同的次序由片段长度多重集字典序确定，
  因此相同请求（无论条目顺序）结果必然一致；
* 整个分析只在内存中进行，不读写数据库。

各切口相互独立时，全部 2^n 个状态的概率之和恒为 1
（数学上）；``math.fsum`` 累加使浮点误差远小于 1e-9。
"""

from __future__ import annotations

import math

from .digest import build_fragments

# 最多枚举的物理切口数：2^12 = 4096 种成功/失败状态
MAX_PARTIAL_CUTS = 12


def partial_digest_distribution(sequence: str, topology: str, cuts, probabilities):
    """枚举切口的成功/失败状态，返回按概率降序合并的片段长度分布。

    ``cuts`` 为已保存结果中的物理切口（dict，含 top_cut/bottom_cut 等），
    ``probabilities`` 与之等长、一一对应，给出每个切口独立的切割概率。
    调用方须已完成全部校验（数量、坐标匹配、概率范围、切口上限）。

    返回 {"total_probability": float, "outcomes": [
        {"fragment_lengths": [int...], "fragment_count": int, "probability": float}
    ]}；``fragment_lengths`` 为排序后的片段长度多重集。
    """
    n = len(cuts)
    # 分组键 = 排序后的片段长度多重集；同一组保留所有状态概率，
    # 最后统一用 fsum 精确累加，避免逐项相加带来的误差累积。
    groups: dict[tuple[int, ...], list[float]] = {}

    for mask in range(1 << n):
        probability = 1.0
        active = []
        for i, cut in enumerate(cuts):
            if mask >> i & 1:
                active.append(cut)
                probability *= probabilities[i]
            else:
                probability *= 1.0 - probabilities[i]
        fragments = build_fragments(sequence, topology, active)
        key = tuple(sorted(f["length"] for f in fragments))
        groups.setdefault(key, []).append(probability)

    outcomes = []
    for key, values in groups.items():
        # 概率为 0 的切口使整条状态链不可能发生（如 p=0/1）；
        # fsum 后仍为 0 的分组是不可能结果，直接剔除
        total = math.fsum(values)
        if total > 0.0:
            outcomes.append(
                {
                    "fragment_lengths": list(key),
                    "fragment_count": len(key),
                    "probability": total,
                }
            )
    # 稳定的全序：先按概率降序，概率相同再按长度多重集字典序
    outcomes.sort(key=lambda o: (-o["probability"], o["fragment_lengths"]))

    return {
        "total_probability": math.fsum(o["probability"] for o in outcomes),
        "outcomes": outcomes,
    }
