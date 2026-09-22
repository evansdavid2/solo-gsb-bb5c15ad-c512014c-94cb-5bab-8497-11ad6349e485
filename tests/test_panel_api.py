"""酶组合面板规划接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


# 线性序列（长度 24）：AAAA GAATTC TTTT GGATCC GGGG
#   EcoRI (GAATTC, 1/5)：top 切口 5   -> 单酶 2 片段
#   BamHI (GGATCC, 1/5)：top 切口 15  -> 单酶 2 片段
#   两酶联用                           -> 3 片段
SEQ = "AAAAGAATTCTTTTGGATCCGGGG"
ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
BAMHI = {"name": "BamHI", "pattern": "GGATCC", "cut_top": 1, "cut_bottom": 5}
ECORV = {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3}  # SEQ 上无位点


def _create(client, sequence=SEQ, topology="linear"):
    resp = client.post(
        "/digests",
        json={"sequence": sequence, "topology": topology, "enzymes": [ECORI]},
    )
    assert resp.status_code == 201
    return resp.json()


def _plan(client, digest_id, target, candidates):
    return client.post(
        f"/digests/{digest_id}/enzyme-panel-plan",
        json={"target_fragment_count": target, "candidates": candidates},
    )


def test_plan_selects_all_enzymes_when_target_needs_both(client):
    body = _create(client)
    resp = _plan(client, body["id"], 3, [ECORI, BAMHI])
    assert resp.status_code == 200
    data = resp.json()

    assert data["id"] == body["id"]
    assert data["target_fragment_count"] == 3
    assert data["fragment_count"] == 3
    # 选中的酶按名称排序返回
    assert [e["name"] for e in data["selected_enzymes"]] == ["BamHI", "EcoRI"]
    # 切口坐标与直接酶切一致
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["cuts"]] == [(5, 9), (15, 19)]


def test_minimal_enzyme_count_wins(client):
    # EcoRI 在该序列上有两个位点（切口 1、11），单酶即可得 3 片段；
    # [EcoRI, EcoRV] 虽同样得到 3 片段，但酶数更多，不应被选中
    body = _create(client, sequence="GAATTCAAAAGAATTCGGATCC")
    resp = _plan(client, body["id"], 3, [ECORV, ECORI, BAMHI])
    assert resp.status_code == 200
    data = resp.json()
    assert [e["name"] for e in data["selected_enzymes"]] == ["EcoRI"]
    assert data["fragment_count"] == 3
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["cuts"]] == [(1, 5), (11, 15)]


def test_tie_break_by_enzyme_name_sequence(client):
    # EcoRI 与 BamHI 单酶均得 2 片段，按酶名序列择优：BamHI < EcoRI
    body = _create(client)
    resp = _plan(client, body["id"], 2, [ECORI, BAMHI])
    assert resp.status_code == 200
    data = resp.json()
    assert [e["name"] for e in data["selected_enzymes"]] == ["BamHI"]
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["cuts"]] == [(15, 19)]


def test_result_independent_of_candidate_order(client):
    body = _create(client)
    r1 = _plan(client, body["id"], 2, [ECORI, BAMHI, ECORV]).json()
    r2 = _plan(client, body["id"], 2, [ECORV, BAMHI, ECORI]).json()
    assert r1 == r2


def test_duplicate_candidates_are_deduped(client):
    # 完全相同的候选去重后参与枚举，不会拼出"同酶两次"的伪组合
    body = _create(client)
    resp = _plan(client, body["id"], 2, [ECORI, dict(ECORI)])
    assert resp.status_code == 200
    assert [e["name"] for e in resp.json()["selected_enzymes"]] == ["EcoRI"]


def test_empty_selection_when_target_is_one(client):
    # 不切割即 1 条片段：最少酶集合为空集
    body = _create(client)
    resp = _plan(client, body["id"], 1, [ECORI, BAMHI])
    assert resp.status_code == 200
    data = resp.json()
    assert data["selected_enzymes"] == []
    assert data["fragment_count"] == 1
    assert data["cuts"] == []


def test_circular_planning(client):
    # 环形序列：EcoRI 单切口开环为 1 片段，加 EcoRV 两个切口得 2 片段
    body = _create(client, sequence="GAATTCGGGGGATATC", topology="circular")
    resp = _plan(client, body["id"], 2, [ECORI, ECORV])
    assert resp.status_code == 200
    data = resp.json()
    assert [e["name"] for e in data["selected_enzymes"]] == ["EcoRI", "EcoRV"]
    assert data["fragment_count"] == 2
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["cuts"]] == [(1, 5), (13, 13)]


def test_no_solution_returns_422(client):
    body = _create(client)
    for impossible_target in (4, 5, 100):  # 两酶最多切出 3 片段
        resp = _plan(client, body["id"], impossible_target, [ECORI, BAMHI])
        assert resp.status_code == 422, (impossible_target, resp.text)


def test_missing_record_returns_404(client):
    resp = _plan(client, 999, 2, [ECORI])
    assert resp.status_code == 404


def test_plan_is_read_only(client):
    body = _create(client)
    resp = _plan(client, body["id"], 3, [ECORI, BAMHI])
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{body['id']}").json() == body
    assert len(client.get("/digests").json()) == 1


def test_too_many_candidates_returns_422(client):
    body = _create(client)
    candidates = [
        {"name": f"E{i}", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
        for i in range(9)
    ]
    assert _plan(client, body["id"], 2, candidates).status_code == 422


def test_empty_candidate_list_returns_422(client):
    body = _create(client)
    assert _plan(client, body["id"], 2, []).status_code == 422


def test_invalid_target_returns_422(client):
    body = _create(client)
    for bad_target in (0, -1):
        assert _plan(client, body["id"], bad_target, [ECORI]).status_code == 422


def test_invalid_candidate_enzyme_returns_422(client):
    body = _create(client)
    bad = {"name": "Bad", "pattern": "GAATTC", "cut_top": 7, "cut_bottom": 5}
    assert _plan(client, body["id"], 2, [ECORI, bad]).status_code == 422
    bad_pattern = {"name": "Bad", "pattern": "GAUZ", "cut_top": 1, "cut_bottom": 1}
    assert _plan(client, body["id"], 2, [bad_pattern]).status_code == 422
