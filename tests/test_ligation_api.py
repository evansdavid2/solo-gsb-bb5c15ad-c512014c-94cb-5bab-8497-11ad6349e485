"""末端连接兼容性接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


# 5' 突出端酶：识别 GAAGTC，切口 1/5 -> 5prime，突出端 AAGT（非回文）
FIVE_PRIME_ENZYME = {"name": "F5", "pattern": "GAAGTC", "cut_top": 1, "cut_bottom": 5}
# 3' 突出端酶：识别 GACTTC，切口 5/1 -> 3prime，突出端 ACTT = 反向互补(AAGT)
THREE_PRIME_ENZYME = {"name": "F3", "pattern": "GACTTC", "cut_top": 5, "cut_bottom": 1}
# 3' 突出端对照酶：识别 GAAGTC，切口 5/1 -> 3prime，突出端 AAGT
# （反向互补(AAGT)=ACTT ≠ AAGT，不应与 5' AAGT 端兼容）
THREE_PRIME_CTRL_ENZYME = {
    "name": "F3C", "pattern": "GAAGTC", "cut_top": 5, "cut_bottom": 1
}
# 平端酶：识别 GATATC，切口 3/3 -> blunt
BLUNT_ENZYME = {"name": "Blunt", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}
# 识别 GAATTC，在全 A 序列上无位点，用于构造"不切割"的记录
ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}


def _create(client, sequence, enzymes, topology="linear"):
    resp = client.post(
        "/digests",
        json={"sequence": sequence, "topology": topology, "enzymes": enzymes},
    )
    assert resp.status_code == 201
    return resp.json()


def _create_five_prime(client, topology="linear"):
    """5prime 记录：序列 AAAAAGAAGTCAAAA，单切口 (5,9)，突出端 AAGT。

    线性片段 [5, 9]：f0 右端 / f1 左端为 5prime(AAGT)；
    环形单切口开环为一条全长片段，左、右两端均为 5prime(AAGT)。
    """
    return _create(client, "AAAAGAAGTCAAAA", [FIVE_PRIME_ENZYME], topology)


def _create_three_prime(client):
    """3prime 记录：单切口 (9,5)（按 top_cut=9 切分），突出端 ACTT。

    片段 [9, 5]：f0 右端 / f1 左端为 3prime(ACTT)。
    """
    return _create(client, "AAAAGACTTCAAAA", [THREE_PRIME_ENZYME])


def _create_blunt(client):
    """blunt 记录：序列 AAAGATATCAAA，单切口 (6,6)，f0 右端 / f1 左端为平端。"""
    return _create(client, "AAAGATATCAAA", [BLUNT_ENZYME])


def _create_uncut(client):
    """无切口记录：唯一片段两端均为 none。"""
    return _create(client, "AAAAAAAAAA", [ECORI])


def _check(client, digest_ids, include_blunt=None):
    payload = {"digest_ids": digest_ids}
    if include_blunt is not None:
        payload["include_blunt"] = include_blunt
    return client.post("/ligation-compatibility", json=payload)


def _end(digest_id, fragment_index, end, end_type, overhang):
    return {
        "digest_id": digest_id,
        "fragment_index": fragment_index,
        "end": end,
        "end_type": end_type,
        "overhang": overhang,
    }


# ---------- 黏端兼容判定 ----------


def test_sticky_pairs_follow_reverse_complement_rule(client):
    rec5 = _create_five_prime(client)  # 5prime AAGT：f0.right、f1.left
    rec3 = _create_three_prime(client)  # 3prime ACTT：f0.right、f1.left
    resp = _check(client, [rec5["id"], rec3["id"]])
    assert resp.status_code == 200
    data = resp.json()
    assert data["digest_ids"] == [rec5["id"], rec3["id"]]
    assert data["include_blunt"] is True
    assert data["record_count"] == 2
    assert data["pair_count"] == 4
    a, b = rec5["id"], rec3["id"]
    assert data["pairs"] == [
        {
            "first": _end(a, 0, "right", "5prime", "AAGT"),
            "second": _end(b, 0, "right", "3prime", "ACTT"),
            "compatibility": "sticky",
            "ligated_overhang": "AAGT",
        },
        {
            "first": _end(a, 0, "right", "5prime", "AAGT"),
            "second": _end(b, 1, "left", "3prime", "ACTT"),
            "compatibility": "sticky",
            "ligated_overhang": "AAGT",
        },
        {
            "first": _end(a, 1, "left", "5prime", "AAGT"),
            "second": _end(b, 0, "right", "3prime", "ACTT"),
            "compatibility": "sticky",
            "ligated_overhang": "AAGT",
        },
        {
            "first": _end(a, 1, "left", "5prime", "AAGT"),
            "second": _end(b, 1, "left", "3prime", "ACTT"),
            "compatibility": "sticky",
            "ligated_overhang": "AAGT",
        },
    ]


def test_sticky_requires_exact_reverse_complement(client):
    rec5 = _create_five_prime(client)  # 5prime AAGT
    ctrl = _create(
        client, "AAAAGAAGTCAAAA", [THREE_PRIME_CTRL_ENZYME]
    )  # 3prime AAGT
    # AAGT != 反向互补(AAGT)=ACTT：序列不互补则不兼容
    assert _check(client, [rec5["id"], ctrl["id"]]).json()["pair_count"] == 0


def test_same_direction_overhangs_never_pair(client):
    rec5a = _create_five_prime(client)
    rec5b = _create_five_prime(client)
    rec3a = _create_three_prime(client)
    rec3b = _create_three_prime(client)
    # 5'-5' 与 3'-3' 均不兼容，即使突出端序列完全相同
    assert _check(client, [rec5a["id"], rec5b["id"]]).json()["pair_count"] == 0
    assert _check(client, [rec3a["id"], rec3b["id"]]).json()["pair_count"] == 0


def test_first_end_comes_from_earlier_requested_record(client):
    rec5 = _create_five_prime(client)
    rec3 = _create_three_prime(client)
    # 反转请求顺序：first 恒属于请求顺序靠前的记录（此处为 3prime 记录）
    data = _check(client, [rec3["id"], rec5["id"]]).json()
    assert data["pair_count"] == 4
    assert all(p["first"]["digest_id"] == rec3["id"] for p in data["pairs"])
    assert all(p["second"]["digest_id"] == rec5["id"] for p in data["pairs"])
    # 连接后突出端始终取该对中 5' 突出端序列，与其位于 first/second 无关
    assert all(p["ligated_overhang"] == "AAGT" for p in data["pairs"])
    assert data["pairs"][0]["first"] == _end(
        rec3["id"], 0, "right", "3prime", "ACTT"
    )
    assert data["pairs"][0]["second"] == _end(
        rec5["id"], 0, "right", "5prime", "AAGT"
    )


def test_pairs_sorted_across_three_records(client):
    rec5 = _create_five_prime(client)
    rec3a = _create_three_prime(client)
    rec3b = _create_three_prime(client)
    data = _check(client, [rec5["id"], rec3a["id"], rec3b["id"]]).json()
    assert data["record_count"] == 3
    # rec3a 与 rec3b 同为 3prime 不兼容，因此 8 对全部是 rec5 × 两条 3prime 记录
    assert data["pair_count"] == 8
    order = [
        (
            p["first"]["digest_id"],
            p["first"]["fragment_index"],
            p["first"]["end"],
            p["second"]["digest_id"],
            p["second"]["fragment_index"],
            p["second"]["end"],
        )
        for p in data["pairs"]
    ]
    a, b, c = rec5["id"], rec3a["id"], rec3b["id"]
    # 先按 first（记录、片段、端方向），再按 second 同一键序：
    # 同一 first 与 rec3a、rec3b 的配对按 second 记录顺序交错排列
    assert order == [
        (a, 0, "right", b, 0, "right"),
        (a, 0, "right", b, 1, "left"),
        (a, 0, "right", c, 0, "right"),
        (a, 0, "right", c, 1, "left"),
        (a, 1, "left", b, 0, "right"),
        (a, 1, "left", b, 1, "left"),
        (a, 1, "left", c, 0, "right"),
        (a, 1, "left", c, 1, "left"),
    ]


def test_circular_record_ends_participate(client):
    circ = _create_five_prime(client, topology="circular")
    # 环形单切口开环：一条全长片段，两端均为 5prime(AAGT)
    assert len(circ["fragments"]) == 1
    rec3 = _create_three_prime(client)
    data = _check(client, [circ["id"], rec3["id"]]).json()
    assert data["pair_count"] == 4
    assert {p["first"]["fragment_index"] for p in data["pairs"]} == {0}
    assert {p["first"]["end"] for p in data["pairs"]} == {"left", "right"}


# ---------- 平端与天然末端 ----------


def test_blunt_pairs_included_by_default(client):
    rec_a = _create_blunt(client)
    rec_b = _create_blunt(client)
    data = _check(client, [rec_a["id"], rec_b["id"]]).json()  # 省略 include_blunt
    assert data["include_blunt"] is True
    assert data["pair_count"] == 4
    assert all(p["compatibility"] == "blunt" for p in data["pairs"])
    assert all(p["ligated_overhang"] == "" for p in data["pairs"])
    a, b = rec_a["id"], rec_b["id"]
    assert data["pairs"][0] == {
        "first": _end(a, 0, "right", "blunt", ""),
        "second": _end(b, 0, "right", "blunt", ""),
        "compatibility": "blunt",
        "ligated_overhang": "",
    }
    assert data["pairs"][-1] == {
        "first": _end(a, 1, "left", "blunt", ""),
        "second": _end(b, 1, "left", "blunt", ""),
        "compatibility": "blunt",
        "ligated_overhang": "",
    }


def test_blunt_pairs_excluded_when_include_blunt_false(client):
    rec_a = _create_blunt(client)
    rec_b = _create_blunt(client)
    data = _check(
        client, [rec_a["id"], rec_b["id"]], include_blunt=False
    ).json()
    assert data["include_blunt"] is False
    assert data["pair_count"] == 0
    assert data["pairs"] == []


def test_include_blunt_false_keeps_sticky_pairs(client):
    rec5 = _create_five_prime(client)
    rec3 = _create_three_prime(client)
    data = _check(
        client, [rec5["id"], rec3["id"]], include_blunt=False
    ).json()
    assert data["pair_count"] == 4
    assert all(p["compatibility"] == "sticky" for p in data["pairs"])


def test_blunt_not_compatible_with_sticky_or_none(client):
    rec_blunt = _create_blunt(client)
    rec5 = _create_five_prime(client)
    rec3 = _create_three_prime(client)
    uncut = _create_uncut(client)
    assert _check(client, [rec_blunt["id"], rec5["id"]]).json()["pair_count"] == 0
    assert _check(client, [rec_blunt["id"], rec3["id"]]).json()["pair_count"] == 0
    assert _check(client, [rec_blunt["id"], uncut["id"]]).json()["pair_count"] == 0


def test_natural_none_ends_never_pair(client):
    uncut_a = _create_uncut(client)
    uncut_b = _create_uncut(client)
    rec5 = _create_five_prime(client)
    assert _check(client, [uncut_a["id"], rec5["id"]]).json()["pair_count"] == 0
    assert _check(client, [uncut_a["id"], uncut_b["id"]]).json()["pair_count"] == 0


def test_ends_within_same_record_are_not_paired(client):
    # 单条记录内同时存在可互补的 5'/3' 突出端：
    # 切口 (5,9) 5'AAGT、(9,5) 3'AAGT、(15,19) 5'ACTT、(19,15) 3'ACTT
    both = _create(
        client,
        "AAAAGAAGTCAAAAGACTTCAAAA",
        [FIVE_PRIME_ENZYME, THREE_PRIME_ENZYME],
    )
    assert len(both["fragments"]) == 5
    uncut = _create_uncut(client)
    # 配对只枚举跨记录组合：与无切口记录一起请求时不存在任何兼容对
    assert _check(client, [both["id"], uncut["id"]]).json()["pair_count"] == 0
    # 该记录自身的 5' 端与另一记录的 3' 端仍可正常跨记录配对
    rec3 = _create_three_prime(client)
    assert _check(client, [both["id"], rec3["id"]]).json()["pair_count"] > 0


# ---------- 记录数上限 ----------


def test_eight_records_supported(client):
    records = [_create_blunt(client) for _ in range(8)]
    data = _check(client, [r["id"] for r in records]).json()
    assert data["record_count"] == 8
    # 每对记录产生 2 × 2 = 4 个平端配对：C(8,2) × 4 = 112
    assert data["pair_count"] == 112


# ---------- 校验与状态码 ----------


def test_duplicate_digest_ids_return_422(client):
    rec = _create_blunt(client)
    resp = _check(client, [rec["id"], rec["id"]])
    assert resp.status_code == 422
    assert "duplicate digest ids" in resp.text


def test_digest_id_count_limits_return_422(client):
    rec = _create_blunt(client)
    # 少于 2 条
    assert _check(client, [rec["id"]]).status_code == 422
    # 多于 8 条（数量校验先于记录存在性检查）
    assert _check(client, list(range(1, 10))).status_code == 422


def test_invalid_argument_types_return_422(client):
    rec = _create_blunt(client)
    assert _check(client, ["x", "y"]).status_code == 422
    # include_blunt 非法布尔值：请求体校验先于记录存在性（rec+1 不存在也不返回 404）
    resp = _check(client, [rec["id"], rec["id"] + 1], include_blunt="notabool")
    assert resp.status_code == 422


def test_missing_record_returns_404(client):
    rec = _create_blunt(client)
    resp = _check(client, [rec["id"], 999])
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_first_missing_id_in_request_order_is_reported(client):
    resp = _check(client, [999, 998])
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_validation_precedes_existence_check(client):
    # 重复 id 直接 422，即使 id 不存在也不报 404
    assert _check(client, [999, 999]).status_code == 422


# ---------- 只读 ----------


def test_analysis_is_read_only(client):
    rec5 = _create_five_prime(client)
    rec3 = _create_three_prime(client)
    resp = _check(client, [rec5["id"], rec3["id"]])
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{rec5['id']}").json() == rec5
    assert client.get(f"/digests/{rec3['id']}").json() == rec3
    assert len(client.get("/digests").json()) == 2
