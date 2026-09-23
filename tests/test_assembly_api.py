"""片段拼接链规划接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


# 5' 突出端酶：识别回文 GAATTC，切口 1/5 -> 5prime，突出端 AATT（回文）
FIVE_PRIME_ENZYME = {"name": "P5", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
# 3' 突出端酶：同一回文模式，切口 5/1 -> 3prime，突出端 AATT
THREE_PRIME_ENZYME = {"name": "P3", "pattern": "GAATTC", "cut_top": 5, "cut_bottom": 1}
# 平端酶：识别 GATATC，切口 3/3 -> blunt
BLUNT_ENZYME = {"name": "Blunt", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}


def _create(client, sequence, enzymes, topology="linear"):
    resp = client.post(
        "/digests",
        json={"sequence": sequence, "topology": topology, "enzymes": enzymes},
    )
    assert resp.status_code == 201
    return resp.json()


def _create_five_prime(client, topology="linear"):
    """5prime 记录：单切口 (5,9)，f0 右端 / f1 左端为 5prime(AATT)。

    环形单切口开环为一条全长片段，左、右两端均为 5prime(AATT)。
    """
    return _create(client, "AAAAGAATTCTTTT", [FIVE_PRIME_ENZYME], topology)


def _create_three_prime(client):
    """3prime 记录：单切口 (9,5)，f0 右端 / f1 左端为 3prime(AATT)。"""
    return _create(client, "AAAAGAATTCTTTT", [THREE_PRIME_ENZYME])


def _create_bridge(client):
    """桥接记录：切口 (5,9) 与 (9,5)，中间片段 f1 左端 5prime、右端 3prime。"""
    return _create(
        client, "AAAAGAATTCTTTT", [FIVE_PRIME_ENZYME, THREE_PRIME_ENZYME]
    )


def _create_blunt(client, sequence="AAAGATATCAAA"):
    """blunt 记录：单切口 (6,6)，f0 右端 / f1 左端为平端。"""
    return _create(client, sequence, [BLUNT_ENZYME])


def _create_uncut(client):
    """无切口记录：唯一片段两端均为 none。"""
    return _create(client, "AAAAAAAAAA", [FIVE_PRIME_ENZYME])


def _plan(client, digest_ids, include_blunt=None, max_chain_length=None):
    payload = {"digest_ids": digest_ids}
    if include_blunt is not None:
        payload["include_blunt"] = include_blunt
    if max_chain_length is not None:
        payload["max_chain_length"] = max_chain_length
    return client.post("/assembly-chains", json=payload)


def _end(digest_id, fragment_index, end, end_type, overhang):
    return {
        "digest_id": digest_id,
        "fragment_index": fragment_index,
        "end": end,
        "end_type": end_type,
        "overhang": overhang,
    }


# ---------- 最长链搜索与择优 ----------


def test_three_record_sticky_chain(client):
    rec5 = _create_five_prime(client)  # 5prime AATT
    bridge = _create_bridge(client)  # f1：左 5prime / 右 3prime
    rec3 = _create_three_prime(client)  # 3prime AATT
    resp = _plan(client, [rec5["id"], bridge["id"], rec3["id"]])
    assert resp.status_code == 200
    data = resp.json()
    assert data["digest_ids"] == [rec5["id"], bridge["id"], rec3["id"]]
    assert data["include_blunt"] is True
    assert data["max_chain_length"] == 8
    assert data["record_count"] == 3
    assert data["chain_length"] == 3
    assert data["junction_count"] == 2
    assert data["stop_reason"] == "all_records_used"
    a, b, c = rec5["id"], bridge["id"], rec3["id"]
    assert data["nodes"] == [
        {
            "digest_id": a,
            "fragment_index": 0,
            "left_end": {"type": "none", "overhang": ""},
            "right_end": {"type": "5prime", "overhang": "AATT"},
        },
        {
            "digest_id": b,
            "fragment_index": 1,
            "left_end": {"type": "5prime", "overhang": "AATT"},
            "right_end": {"type": "3prime", "overhang": "AATT"},
        },
        {
            "digest_id": c,
            "fragment_index": 0,
            "left_end": {"type": "none", "overhang": ""},
            "right_end": {"type": "3prime", "overhang": "AATT"},
        },
    ]
    assert data["junctions"] == [
        {
            "from_end": _end(a, 0, "right", "5prime", "AATT"),
            "to_end": _end(b, 1, "right", "3prime", "AATT"),
            "compatibility": "sticky",
            "ligated_overhang": "AATT",
        },
        {
            "from_end": _end(b, 1, "left", "5prime", "AATT"),
            "to_end": _end(c, 0, "right", "3prime", "AATT"),
            "compatibility": "sticky",
            "ligated_overhang": "AATT",
        },
    ]


def test_request_order_determines_chain_direction(client):
    rec5 = _create_five_prime(client)
    bridge = _create_bridge(client)
    rec3 = _create_three_prime(client)
    # 反转请求顺序：首节点变为请求顺序最靠前的 3prime 记录，链方向随之反转
    data = _plan(client, [rec3["id"], bridge["id"], rec5["id"]]).json()
    assert data["chain_length"] == 3
    assert data["stop_reason"] == "all_records_used"
    a, b, c = rec3["id"], bridge["id"], rec5["id"]
    assert [n["digest_id"] for n in data["nodes"]] == [a, b, c]
    assert data["junctions"] == [
        {
            "from_end": _end(a, 0, "right", "3prime", "AATT"),
            "to_end": _end(b, 1, "left", "5prime", "AATT"),
            "compatibility": "sticky",
            "ligated_overhang": "AATT",
        },
        {
            "from_end": _end(b, 1, "right", "3prime", "AATT"),
            "to_end": _end(c, 0, "right", "5prime", "AATT"),
            "compatibility": "sticky",
            "ligated_overhang": "AATT",
        },
    ]


def test_tie_break_prefers_earlier_requested_record(client):
    rec5a = _create_five_prime(client)
    rec5b = _create_five_prime(client)
    rec3 = _create_three_prime(client)
    # 两条 5prime 记录等价：同长度时选请求顺序靠前者
    data = _plan(client, [rec5a["id"], rec5b["id"], rec3["id"]]).json()
    assert data["chain_length"] == 2
    assert data["nodes"][0]["digest_id"] == rec5a["id"]
    assert data["nodes"][1]["digest_id"] == rec3["id"]
    # 反转请求顺序后改选 rec5b
    data = _plan(client, [rec5b["id"], rec5a["id"], rec3["id"]]).json()
    assert data["nodes"][0]["digest_id"] == rec5b["id"]


def test_tie_break_prefers_earlier_fragment_and_left_end(client):
    # 环形单切口开环：唯一片段两端均为 5prime(AATT)，left 先于 right 入选
    circ = _create_five_prime(client, topology="circular")
    rec3 = _create_three_prime(client)
    data = _plan(client, [circ["id"], rec3["id"]]).json()
    assert data["chain_length"] == 2
    assert data["junctions"][0]["from_end"] == _end(
        circ["id"], 0, "left", "5prime", "AATT"
    )
    assert data["junctions"][0]["to_end"] == _end(
        rec3["id"], 0, "right", "3prime", "AATT"
    )


def test_record_used_at_most_once(client):
    circ = _create_five_prime(client, topology="circular")  # 两端均 5prime
    rec3 = _create_three_prime(client)
    # 若允许重复使用记录，可得到 rec3 -> circ -> rec3 的 3 节点链
    data = _plan(client, [circ["id"], rec3["id"]]).json()
    assert data["chain_length"] == 2
    assert data["junction_count"] == 1
    assert data["stop_reason"] == "all_records_used"


def test_circular_fragment_bridges_chain(client):
    circ = _create_five_prime(client, topology="circular")
    rec3a = _create_three_prime(client)
    rec3b = _create_three_prime(client)
    data = _plan(client, [circ["id"], rec3a["id"], rec3b["id"]]).json()
    assert data["chain_length"] == 3
    assert data["junction_count"] == 2
    assert data["stop_reason"] == "all_records_used"
    # 环形开环片段作为中间节点：left 接前驱、right 接后继
    assert data["nodes"][1] == {
        "digest_id": circ["id"],
        "fragment_index": 0,
        "left_end": {"type": "5prime", "overhang": "AATT"},
        "right_end": {"type": "5prime", "overhang": "AATT"},
    }
    assert data["junctions"][0]["to_end"]["end"] == "left"
    assert data["junctions"][1]["from_end"]["end"] == "right"


# ---------- 最长链长上限与停止原因 ----------


def test_max_chain_length_truncates_chain(client):
    rec5 = _create_five_prime(client)
    bridge = _create_bridge(client)
    rec3 = _create_three_prime(client)
    data = _plan(
        client,
        [rec5["id"], bridge["id"], rec3["id"]],
        max_chain_length=2,
    ).json()
    assert data["max_chain_length"] == 2
    assert data["chain_length"] == 2
    assert data["junction_count"] == 1
    assert data["stop_reason"] == "max_chain_length_reached"
    # 长度 2 的方案独立择优：第二节点取全序最小的兼容端（bridge f1 右端）
    assert [n["digest_id"] for n in data["nodes"]] == [rec5["id"], bridge["id"]]
    assert data["junctions"][0]["to_end"] == _end(
        bridge["id"], 1, "right", "3prime", "AATT"
    )


def test_max_chain_length_one_yields_single_node(client):
    rec5 = _create_five_prime(client)
    bridge = _create_bridge(client)
    rec3 = _create_three_prime(client)
    data = _plan(
        client,
        [rec5["id"], bridge["id"], rec3["id"]],
        max_chain_length=1,
    ).json()
    assert data["chain_length"] == 1
    assert data["junction_count"] == 0
    assert data["junctions"] == []
    assert data["stop_reason"] == "max_chain_length_reached"
    # 单节点取全序最小片段：请求顺序最靠前记录的第一个片段
    assert data["nodes"] == [
        {
            "digest_id": rec5["id"],
            "fragment_index": 0,
            "left_end": {"type": "none", "overhang": ""},
            "right_end": {"type": "5prime", "overhang": "AATT"},
        }
    ]


def test_no_compatible_end_returns_single_node(client):
    uncut = _create_uncut(client)  # 两端均为 none，无法参与任何连接
    rec5 = _create_five_prime(client)
    data = _plan(client, [uncut["id"], rec5["id"]]).json()
    assert data["chain_length"] == 1
    assert data["junction_count"] == 0
    assert data["junctions"] == []
    assert data["stop_reason"] == "no_compatible_end"
    assert data["nodes"][0]["digest_id"] == uncut["id"]
    assert data["nodes"][0]["left_end"] == {"type": "none", "overhang": ""}
    assert data["nodes"][0]["right_end"] == {"type": "none", "overhang": ""}


# ---------- 平端连接 ----------


def test_blunt_chain_included_by_default(client):
    rec_a = _create_blunt(client)
    rec_b = _create_blunt(client)
    data = _plan(client, [rec_a["id"], rec_b["id"]]).json()  # 省略 include_blunt
    assert data["include_blunt"] is True
    assert data["chain_length"] == 2
    assert data["junction_count"] == 1
    assert data["stop_reason"] == "all_records_used"
    assert data["junctions"] == [
        {
            "from_end": _end(rec_a["id"], 0, "right", "blunt", ""),
            "to_end": _end(rec_b["id"], 0, "right", "blunt", ""),
            "compatibility": "blunt",
            "ligated_overhang": "",
        }
    ]


def test_blunt_chain_excluded_when_include_blunt_false(client):
    rec_a = _create_blunt(client)
    rec_b = _create_blunt(client)
    data = _plan(
        client, [rec_a["id"], rec_b["id"]], include_blunt=False
    ).json()
    assert data["include_blunt"] is False
    assert data["chain_length"] == 1
    assert data["junction_count"] == 0
    assert data["stop_reason"] == "no_compatible_end"


def test_eight_record_chain(client):
    # 双位点序列：中间片段 f1 两端均为平端，可任意串联
    records = [
        _create_blunt(client, sequence="AAAGATATCAAAGATATCAAA") for _ in range(8)
    ]
    data = _plan(client, [r["id"] for r in records]).json()
    assert data["record_count"] == 8
    assert data["chain_length"] == 8
    assert data["junction_count"] == 7
    assert data["stop_reason"] == "all_records_used"
    assert all(j["compatibility"] == "blunt" for j in data["junctions"])
    # 每条记录恰好出现一次
    assert sorted(n["digest_id"] for n in data["nodes"]) == sorted(
        r["id"] for r in records
    )


# ---------- 校验与状态码 ----------


def test_duplicate_digest_ids_return_422(client):
    rec = _create_blunt(client)
    resp = _plan(client, [rec["id"], rec["id"]])
    assert resp.status_code == 422
    assert "duplicate digest ids" in resp.text


def test_digest_id_count_limits_return_422(client):
    rec = _create_blunt(client)
    # 少于 2 条
    assert _plan(client, [rec["id"]]).status_code == 422
    # 多于 8 条（数量校验先于记录存在性检查）
    assert _plan(client, list(range(1, 10))).status_code == 422


def test_invalid_max_chain_length_returns_422(client):
    rec_a = _create_blunt(client)
    rec_b = _create_blunt(client)
    ids = [rec_a["id"], rec_b["id"]]
    assert _plan(client, ids, max_chain_length=0).status_code == 422
    assert _plan(client, ids, max_chain_length=-1).status_code == 422
    assert _plan(client, ids, max_chain_length="x").status_code == 422


def test_invalid_argument_types_return_422(client):
    rec_a = _create_blunt(client)
    rec_b = _create_blunt(client)
    assert _plan(client, ["x", "y"]).status_code == 422
    # include_blunt 非法布尔值：请求体校验先于记录存在性
    resp = _plan(client, [rec_a["id"], rec_b["id"]], include_blunt="notabool")
    assert resp.status_code == 422


def test_missing_record_returns_404(client):
    rec = _create_blunt(client)
    resp = _plan(client, [rec["id"], 999])
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_first_missing_id_in_request_order_is_reported(client):
    resp = _plan(client, [999, 998])
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_validation_precedes_existence_check(client):
    # 重复 id 直接 422，即使 id 不存在也不报 404
    assert _plan(client, [999, 999]).status_code == 422


# ---------- 只读 ----------


def test_planning_is_read_only(client):
    rec5 = _create_five_prime(client)
    bridge = _create_bridge(client)
    rec3 = _create_three_prime(client)
    resp = _plan(client, [rec5["id"], bridge["id"], rec3["id"]])
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{rec5['id']}").json() == rec5
    assert client.get(f"/digests/{bridge['id']}").json() == bridge
    assert client.get(f"/digests/{rec3['id']}").json() == rec3
    assert len(client.get("/digests").json()) == 3
