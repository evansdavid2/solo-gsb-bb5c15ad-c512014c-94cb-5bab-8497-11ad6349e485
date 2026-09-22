"""部分酶切概率分布接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
# 回文模式 GC 的一位切口酶（cut 1/1 平端），只在 GC 位点产生切口，
# 不会像单碱基模式 A 那样因反向互补为 T 而在混合序列上到处切割
GC_ENZYME = {"name": "GC-cutter", "pattern": "GC", "cut_top": 1, "cut_bottom": 1}

# 长度 14：EcoRI 唯一物理切口 (5, 9)，完整酶切片段为 [5, 9]
ECO_SEQ = "AAAAGAATTCTTTT"
# 长度 9、恰有两个对称 GC 位点：物理切口 (2, 2) 与 (7, 7)
SYMMETRIC_SEQ = "TGCAAAGCT"
# 长度 10、三个等距 GC 位点：物理切口 (2,2),(5,5),(8,8)
THREE_CUT_SEQ = "TGCAGCAGCT"


def _create(client, sequence, topology="linear", enzymes=None):
    resp = client.post(
        "/digests",
        json={
            "sequence": sequence,
            "topology": topology,
            "enzymes": enzymes if enzymes is not None else [GC_ENZYME],
        },
    )
    assert resp.status_code == 201
    return resp.json()


def _dist(client, digest_id, cuts):
    return client.post(
        f"/digests/{digest_id}/partial-digest-distribution",
        json={"cuts": cuts},
    )


def _probs_for(body, prob=0.5):
    """为已保存的每个切口生成等概率请求。"""
    return [
        {"top_cut": c["top_cut"], "bottom_cut": c["bottom_cut"], "probability": prob}
        for c in body["cuts"]
    ]


def test_zero_probabilities_yields_uncut_molecule(client):
    body = _create(client, ECO_SEQ, enzymes=[ECORI])
    assert [(c["top_cut"], c["bottom_cut"]) for c in body["cuts"]] == [(5, 9)]
    resp = _dist(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": 0.0}])
    assert resp.status_code == 200
    data = resp.json()
    assert data["cut_count"] == 1
    assert data["outcome_count"] == 1
    assert data["outcomes"] == [
        {"fragment_lengths": [14], "fragment_count": 1, "probability": 1.0}
    ]
    assert data["total_probability"] == pytest.approx(1.0, abs=1e-12)


def test_probability_one_reproduces_full_digest(client):
    body = _create(client, ECO_SEQ, enzymes=[ECORI])
    resp = _dist(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": 1.0}])
    assert resp.status_code == 200
    lengths = sorted(f["length"] for f in body["fragments"])
    assert resp.json()["outcomes"] == [
        {"fragment_lengths": lengths, "fragment_count": 2, "probability": 1.0}
    ]


def test_single_cut_half_probability_distribution(client):
    # p=0.5：切 [5,9] 与不切 [14] 各 0.5；等概率时按长度多重集字典序
    body = _create(client, ECO_SEQ, enzymes=[ECORI])
    resp = _dist(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": 0.5}])
    assert resp.status_code == 200
    assert resp.json()["outcomes"] == [
        {"fragment_lengths": [5, 9], "fragment_count": 2, "probability": 0.5},
        {"fragment_lengths": [14], "fragment_count": 1, "probability": 0.5},
    ]


def test_unequal_probabilities_sorted_descending(client):
    body = _create(client, ECO_SEQ, enzymes=[ECORI])
    resp = _dist(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": 0.2}])
    outcomes = resp.json()["outcomes"]
    # 不切概率 0.8 在前，切概率 0.2 在后
    assert [o["probability"] for o in outcomes] == [0.8, 0.2]
    assert outcomes[0]["fragment_lengths"] == [14]
    assert outcomes[1]["fragment_lengths"] == [5, 9]


def test_symmetric_cuts_merge_equal_length_multisets_linear(client):
    # 切口 2 与 7（长度 9）：
    #   都不切 -> [9]        0.25
    #   只切 2 -> [2,7]      0.25  ┐ 长度多重集相同，合并为 0.5
    #   只切 7 -> [7,2]      0.25  ┘
    #   都切   -> [2,2,5]    0.25
    body = _create(client, SYMMETRIC_SEQ)
    assert [(c["top_cut"], c["bottom_cut"]) for c in body["cuts"]] == [(2, 2), (7, 7)]
    resp = _dist(client, body["id"], [
        {"top_cut": 2, "bottom_cut": 2, "probability": 0.5},
        {"top_cut": 7, "bottom_cut": 7, "probability": 0.5},
    ])
    assert resp.status_code == 200
    data = resp.json()
    assert data["outcome_count"] == 3  # 4 个状态合并为 3 种结果
    assert data["outcomes"] == [
        {"fragment_lengths": [2, 7], "fragment_count": 2, "probability": 0.5},
        {"fragment_lengths": [2, 2, 5], "fragment_count": 3, "probability": 0.25},
        {"fragment_lengths": [9], "fragment_count": 1, "probability": 0.25},
    ]


def test_total_probability_within_tolerance_for_twelve_cuts(client):
    # 线性 12 个 A：单碱基酶 A（纯 A 序列上反向等价）物理切口在 1..12
    # 共 12 个（12 为端点切口）；4096 个状态、各 1/3 概率，总概率误差须 < 1e-9
    a_enzyme = {"name": "A-cutter", "pattern": "A", "cut_top": 1, "cut_bottom": 1}
    body = _create(client, "A" * 12, enzymes=[a_enzyme])
    keys = [(c["top_cut"], c["bottom_cut"]) for c in body["cuts"]]
    assert keys == [(i, i) for i in range(1, 13)]
    resp = _dist(client, body["id"], [
        {"top_cut": t, "bottom_cut": b, "probability": 1.0 / 3.0} for t, b in keys
    ])
    assert resp.status_code == 200
    data = resp.json()
    # 内部切口 1..11 的每个子集对应 12 的一种有序拆分，
    # 合并后的长度多重集数 = 12 的整数分拆数 p(12) = 77
    assert data["outcome_count"] == 77
    assert abs(data["total_probability"] - 1.0) < 1e-9
    # 概率严格按降序返回
    probs = [o["probability"] for o in data["outcomes"]]
    assert probs == sorted(probs, reverse=True)


def test_circular_distribution(client):
    # 环形 "CGCAAAG"（长度 7）：GC 位点恰为两个物理切口 (0,0)（跨原点）
    # 与 (2,2)
    #   无切口        -> 整条环 [7]   0.25
    #   仅单切口（任） -> 开环 [7]    0.50（两种状态合并）
    #   双切口        -> 环形两段 [2,5] 0.25
    body = _create(client, "CGCAAAG", topology="circular")
    keys = {(c["top_cut"], c["bottom_cut"]) for c in body["cuts"]}
    assert keys == {(0, 0), (2, 2)}
    resp = _dist(client, body["id"], [
        {"top_cut": 0, "bottom_cut": 0, "probability": 0.5},
        {"top_cut": 2, "bottom_cut": 2, "probability": 0.5},
    ])
    assert resp.status_code == 200
    data = resp.json()
    by_lengths = {tuple(o["fragment_lengths"]): o["probability"] for o in data["outcomes"]}
    assert by_lengths[(7,)] == pytest.approx(0.75)
    assert by_lengths[(2, 5)] == pytest.approx(0.25)
    assert abs(data["total_probability"] - 1.0) < 1e-12


def test_result_independent_of_request_entry_order(client):
    body = _create(client, THREE_CUT_SEQ)
    assert [(c["top_cut"], c["bottom_cut"]) for c in body["cuts"]] == [
        (2, 2), (5, 5), (8, 8)
    ]
    entries = [
        {"top_cut": 2, "bottom_cut": 2, "probability": 0.8},
        {"top_cut": 5, "bottom_cut": 5, "probability": 0.3},
        {"top_cut": 8, "bottom_cut": 8, "probability": 0.6},
    ]
    r1 = _dist(client, body["id"], entries).json()
    r2 = _dist(client, body["id"], list(reversed(entries))).json()
    assert r1 == r2


def test_record_without_cuts_accepts_empty_list(client):
    # 序列上没有 EcoRI 位点 -> 无已保存切口，空请求得到唯一的全长结果
    body = _create(client, "CCCC", enzymes=[ECORI])
    assert body["cuts"] == []
    resp = _dist(client, body["id"], [])
    assert resp.status_code == 200
    data = resp.json()
    assert data["cut_count"] == 0
    assert data["outcomes"] == [
        {"fragment_lengths": [4], "fragment_count": 1, "probability": 1.0}
    ]
    assert abs(data["total_probability"] - 1.0) < 1e-12


def test_missing_record_returns_404(client):
    assert _dist(client, 999, []).status_code == 404


def test_missing_probability_returns_422(client):
    body = _create(client, SYMMETRIC_SEQ)
    # 已保存 2 个切口，只给了 1 个
    resp = _dist(client, body["id"], [
        {"top_cut": 2, "bottom_cut": 2, "probability": 0.5}
    ])
    assert resp.status_code == 422
    assert "(7, 7)" in resp.json()["detail"]


def test_extra_probability_returns_422(client):
    body = _create(client, SYMMETRIC_SEQ)
    resp = _dist(client, body["id"], _probs_for(body) + [
        {"top_cut": 8, "bottom_cut": 8, "probability": 0.5}
    ])
    assert resp.status_code == 422
    assert "(8, 8)" in resp.json()["detail"]


def test_missing_and_extra_together_return_422(client):
    body = _create(client, SYMMETRIC_SEQ)
    resp = _dist(client, body["id"], [
        {"top_cut": 7, "bottom_cut": 7, "probability": 0.5},
        {"top_cut": 8, "bottom_cut": 8, "probability": 0.5},
    ])
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "(2, 2)" in detail and "(8, 8)" in detail


def test_duplicate_coordinates_return_422(client):
    body = _create(client, ECO_SEQ, enzymes=[ECORI])
    resp = _dist(client, body["id"], [
        {"top_cut": 5, "bottom_cut": 9, "probability": 0.5},
        {"top_cut": 5, "bottom_cut": 9, "probability": 0.2},
    ])
    assert resp.status_code == 422


def test_invalid_probability_returns_422(client):
    body = _create(client, ECO_SEQ, enzymes=[ECORI])
    for bad in (-0.01, 1.01):
        resp = _dist(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": bad}])
        assert resp.status_code == 422, bad


def test_nan_probability_returns_422(client):
    body = _create(client, ECO_SEQ, enzymes=[ECORI])
    resp = _dist(client, body["id"], [{"top_cut": 5, "bottom_cut": 9, "probability": "NaN"}])
    assert resp.status_code == 422


def test_more_than_twelve_saved_cuts_returns_422(client):
    # 线性 14 个 A：单碱基酶在 1..14 共 14 个物理切口，超出 12 的枚举上限
    a_enzyme = {"name": "A-cutter", "pattern": "A", "cut_top": 1, "cut_bottom": 1}
    body = _create(client, "A" * 14, enzymes=[a_enzyme])
    assert len(body["cuts"]) == 14
    # 即使请求本身为空，也先因记录切口超限被拒绝
    assert _dist(client, body["id"], []).status_code == 422
    resp = _dist(client, body["id"], _probs_for(body))
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "14" in detail and "12" in detail


def test_exactly_twelve_saved_cuts_is_allowed(client):
    a_enzyme = {"name": "A-cutter", "pattern": "A", "cut_top": 1, "cut_bottom": 1}
    body = _create(client, "A" * 12, enzymes=[a_enzyme])
    assert len(body["cuts"]) == 12
    resp = _dist(client, body["id"], _probs_for(body))
    assert resp.status_code == 200
    assert abs(resp.json()["total_probability"] - 1.0) < 1e-9


def test_analysis_is_read_only(client):
    body = _create(client, SYMMETRIC_SEQ)
    resp = _dist(client, body["id"], _probs_for(body, 0.3))
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{body['id']}").json() == body
    assert len(client.get("/digests").json()) == 1
