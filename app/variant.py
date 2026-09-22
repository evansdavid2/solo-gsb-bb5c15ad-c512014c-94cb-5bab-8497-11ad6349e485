"""变异影响分析：在不改动已保存记录的前提下重跑酶切并对比差异。

对比口径
--------
* 切口以物理坐标 ``(top_cut, bottom_cut)`` 作为身份：
  - 仅出现在变异后结果中的切口 -> added（新增）；
  - 仅出现在原结果中的切口     -> removed（消失）；
  - 坐标相同但末端（end_type/overhang）或归属酶/位点变化 -> changed；
  - 切口"移动位置"（坐标改变）表现为一条 removed + 一条 added。
* 片段以其在骨架上的起点 ``start`` 配对（替换不改变序列长度，起点坐标
  在变异前后含义一致）：起点仅一侧存在即为 added/removed；起点相同而
  长度不同则记入 changed 并给出 length_delta。
"""

from __future__ import annotations

# 切口条目上参与"末端/归属变化"判定的字段
_CUT_CHANGE_FIELDS = ("end_type", "overhang", "enzymes", "sites")


def apply_substitutions(sequence: str, variants) -> str:
    """按 0-based 单碱基替换生成变异序列（调用方须已完成全部校验）。

    variants 为含 ``position`` / ``alt`` 属性的对象列表，按位置升序施加，
    结果与传入顺序无关。
    """
    chars = list(sequence)
    for v in sorted(variants, key=lambda v: v.position):
        chars[v.position] = v.alt
    return "".join(chars)


def _cut_key(cut: dict) -> tuple[int, int]:
    return cut["top_cut"], cut["bottom_cut"]


def diff_cuts(before: list[dict], after: list[dict]) -> dict:
    """对比两次酶切的切口列表，返回稳定排序后的 added/removed/changed。"""
    before_map = {_cut_key(c): c for c in before}
    after_map = {_cut_key(c): c for c in after}
    before_keys, after_keys = set(before_map), set(after_map)

    added = [after_map[k] for k in sorted(after_keys - before_keys)]
    removed = [before_map[k] for k in sorted(before_keys - after_keys)]

    changed = []
    for key in sorted(before_keys & after_keys):
        old, new = before_map[key], after_map[key]
        fields = [f for f in _CUT_CHANGE_FIELDS if old[f] != new[f]]
        if fields:
            changed.append({"before": old, "after": new, "changes": fields})
    return {"added": added, "removed": removed, "changed": changed}


def diff_fragments(before: list[dict], after: list[dict]) -> dict:
    """按片段起点配对，返回新增/消失片段与长度差异。

    before_lengths/after_lengths 保持片段在各自结果中的物理顺序；
    added/removed/changed 均按起点坐标排序，保证输出确定。
    """
    before_map = {f["start"]: f for f in before}
    after_map = {f["start"]: f for f in after}

    added = [after_map[s] for s in sorted(set(after_map) - set(before_map))]
    removed = [before_map[s] for s in sorted(set(before_map) - set(after_map))]

    changed = []
    for start in sorted(set(before_map) & set(after_map)):
        old_len = before_map[start]["length"]
        new_len = after_map[start]["length"]
        if old_len != new_len:
            changed.append(
                {
                    "start": start,
                    "before_length": old_len,
                    "after_length": new_len,
                    "length_delta": new_len - old_len,
                }
            )

    return {
        "before_lengths": [f["length"] for f in before],
        "after_lengths": [f["length"] for f in after],
        "added": added,
        "removed": removed,
        "changed": changed,
    }


def variant_summary(sequence: str, variants) -> list[dict]:
    """构造回显用的变异列表（按位置稳定排序，附变异前碱基 ref）。"""
    return [
        {"position": v.position, "ref": sequence[v.position], "alt": v.alt}
        for v in sorted(variants, key=lambda v: v.position)
    ]
