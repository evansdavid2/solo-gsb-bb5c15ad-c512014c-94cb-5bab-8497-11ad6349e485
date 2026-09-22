"""酶组合面板规划：枚举候选酶组合，寻找达到目标片段数的最少酶集合。

规划规则
--------
* 复用酶切核心算法 ``run_digest``，对候选酶的每个子集计算片段数；
* 优先选择酶数量最少的组合（空子集也参与枚举：线性/环形分子不切割
  均为 1 条片段，因此目标片段数为 1 时可能返回空的酶集合）；
* 酶数相同时，按酶名序列（排序后的名字元组）字典序最小者胜出；
* 完全相同（名称/模式/偏移均相同）的候选先去重，避免无意义的重复组合；
* 整个规划只在内存中枚举计算，不读写数据库，原记录保持不变。
"""

from __future__ import annotations

from itertools import combinations

from .digest import run_digest


def _fragment_count_and_cuts(sequence: str, topology: str, enzymes) -> tuple[int, list]:
    """复用酶切计算，返回 (片段数, 切口列表)。

    空酶列表不切割：无论线性还是环形，均为 1 条全长片段、无切口
    （run_digest 不支持空酶列表，在此单独处理）。
    """
    if not enzymes:
        return 1, []
    result = run_digest(sequence, topology, enzymes)
    return len(result["fragments"]), result["cuts"]


def plan_enzyme_panel(sequence: str, topology: str, candidates, target_fragment_count: int):
    """枚举候选酶子集，返回达到目标片段数的最少酶组合，无解时返回 None。

    返回 {"enzymes": [...], "cuts": [...], "fragment_count": int}；
    enzymes 按 (name, pattern, cut_top, cut_bottom) 排序，保证相同输入的
    输出与候选在请求中的传入顺序无关。
    """
    # 完全相同的候选去重，并规范化排序，使枚举顺序确定
    unique = sorted(
        {(e.name, e.pattern, e.cut_top, e.cut_bottom): e for e in candidates}.values(),
        key=lambda e: (e.name, e.pattern, e.cut_top, e.cut_bottom),
    )

    for size in range(0, len(unique) + 1):
        best = None  # (酶名序列, 组合, 切口列表)
        for combo in combinations(unique, size):
            fragment_count, cuts = _fragment_count_and_cuts(sequence, topology, combo)
            if fragment_count != target_fragment_count:
                continue
            names = tuple(e.name for e in combo)  # combo 已按名称排序
            if best is None or names < best[0]:
                best = (names, combo, cuts)
        if best is not None:
            # 当前酶数已存在解，即为最少酶集合，不再枚举更大的组合
            _, combo, cuts = best
            return {
                "enzymes": list(combo),
                "cuts": cuts,
                "fragment_count": target_fragment_count,
            }
    return None
