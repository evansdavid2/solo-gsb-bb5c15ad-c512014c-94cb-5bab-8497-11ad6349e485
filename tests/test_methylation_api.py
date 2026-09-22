"""甲基化阻断分析接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
# 非回文酶（GGATC，cut 1/4）：反向命中 GATCC，用于验证反向位点坐标映射
FIVER = {"name": "Fiver", "pattern": "GGATC", "cut_top": 1, "cut_bottom": 4}
# 回文平端酶（GC，cut 1/1），用于环形片段重建
GC_ENZYME = {"name": "GC-cutter", "pattern": "GC", "cut_top": 1, "cut_bottom": 1}
# 两种同模式、同切口但名称不同的回文平端酶（GGCC，cut 2/2，切口 (4,4)）：
# 模式相同故位点完全重合，可通过不同敏感规则分别阻断，验证共享物理切口的保留/移除
ENZ_A = {"name": "EnzA", "pattern": "GGCC", "cut_top": 2, "cut_bottom": 2}
ENZ_B = {"name": "EnzB", "pattern": "GGCC", "cut_top": 2, "cut_bottom": 2}

ECO_SEQ = "AAAAGAATTCTTTT"  # 长度 14，EcoRI 位点起于 4，切口 (5, 9)
TWO_SITE_SEQ = "GAATTCAAAAGAATTC"  # 长度 16，EcoRI 位点起于 0、10
# AAGGCCAA：GGCC 起于 2，EnzA/EnzB 共享唯一物理切口 (4,4)
SHARED_CUT_SEQ = "AAGGCCAA"


def _create(client, sequence, topology="linear", enzymes=None):
    resp = client.post(
        "/digests",
        json={
            "sequence": sequence,
            "topology": topology,
            "enzymes": enzymes if enzymes is not None else [ECORI],
        },
    )
    assert resp.status_code == 201
    return resp.json()


def _analyze(client, digest_id, methylated_positions, enzyme_rules):
    return client.post(
        f"/digests/{digest_id}/methylation-analysis",
        json={
            "methylated_positions": methylated_positions,
            "enzyme_rules": enzyme_rules,
        },
    )


# ---------- 基本阻断行为 ----------


def test_methylation_on_sensitive_offset_blocks_site_and_cut(client):
    body = _create(client, ECO_SEQ)
    # EcoRI 模式 GAATTC 起于 4；偏移 0 -> 物理坐标 4（G）
    resp = _analyze(client, body["id"], [4], [{"enzyme": "EcoRI", "sensitive_offsets": [0]}])
    assert resp.status_code == 200
    data = resp.json()

    assert data["methylated_positions"] == [4]
    assert len(data["blocked_sites"]) == 1
    blocked = data["blocked_sites"][0]
    assert blocked["top_cut"] == 5 and blocked["bottom_cut"] == 9
    assert blocked["site"] == {
        "enzyme": "EcoRI", "start": 4, "strand": "+", "matched": "GAATTC"
    }
    assert blocked["methylated_positions"] == [4]

    assert data["retained_cuts"] == []
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["removed_cuts"]] == [(5, 9)]
    # 切口全部消失 -> 片段重建为整条未切分子
    assert [f["length"] for f in data["fragments"]] == [14]


def test_methylation_on_non_sensitive_base_leaves_cut_intact(client):
    body = _create(client, ECO_SEQ)
    # 仅偏移 0 敏感；甲基化落在偏移 1（物理坐标 5）不构成阻断
    resp = _analyze(client, body["id"], [5], [{"enzyme": "EcoRI", "sensitive_offsets": [0]}])
    assert resp.status_code == 200
    data = resp.json()
    assert data["blocked_sites"] == []
    assert data["removed_cuts"] == []
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["retained_cuts"]] == [(5, 9)]
    assert [f["length"] for f in data["fragments"]] == [5, 9]


def test_methylation_outside_any_site_does_not_block(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(client, body["id"], [0], [{"enzyme": "EcoRI", "sensitive_offsets": [0]}])
    assert resp.status_code == 200
    data = resp.json()
    assert data["blocked_sites"] == []
    assert len(data["retained_cuts"]) == 1


def test_any_one_sensitive_position_blocks_site(client):
    body = _create(client, ECO_SEQ)
    # 给偏移 0..5 全部敏感，但只有物理坐标 7（偏移 3）被甲基化，位点仍须阻断
    resp = _analyze(
        client, body["id"], [7],
        [{"enzyme": "EcoRI", "sensitive_offsets": [0, 1, 2, 3, 4, 5]}],
    )
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["methylated_positions"] == [7]
    assert len(data["removed_cuts"]) == 1


# ---------- 多切口与共享物理切口 ----------


def test_one_of_two_sites_blocked_other_retained_and_fragments_rebuilt(client):
    body = _create(client, TWO_SITE_SEQ)
    assert [(c["top_cut"], c["bottom_cut"]) for c in body["cuts"]] == [(1, 5), (11, 15)]
    # 仅阻断起于 0 的位点（偏移 0 -> 坐标 0）
    resp = _analyze(client, body["id"], [0], [{"enzyme": "EcoRI", "sensitive_offsets": [0]}])
    assert resp.status_code == 200
    data = resp.json()

    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["site"]["start"] == 0
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["retained_cuts"]] == [(11, 15)]
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["removed_cuts"]] == [(1, 5)]
    # 只剩切口 11：线性片段为 [0,11) 与 [11,16)
    assert [f["length"] for f in data["fragments"]] == [11, 5]
    assert data["fragments"][0]["left_end"]["type"] == "none"
    assert data["fragments"][0]["right_end"]["type"] == "5prime"


def test_shared_physical_cut_kept_while_one_enzyme_site_blocked(client):
    # GGCC 起于 2，EnzA/EnzB 共享物理切口 (4,4)：
    # 甲基化坐标 2 仅阻断对偏移 0 敏感的 EnzA；EnzB 不敏感，切口保留并重算
    body = _create(client, SHARED_CUT_SEQ, enzymes=[ENZ_A, ENZ_B])
    cut = body["cuts"][0]
    assert (cut["top_cut"], cut["bottom_cut"]) == (4, 4)
    assert cut["enzymes"] == ["EnzA", "EnzB"]

    resp = _analyze(
        client, body["id"], [2],
        [{"enzyme": "EnzA", "sensitive_offsets": [0]}],
    )
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["site"]["enzyme"] == "EnzA"

    # 物理切口仍由 EnzB 维持：保留并重算 enzymes/sites
    assert data["removed_cuts"] == []
    assert len(data["retained_cuts"]) == 1
    retained = data["retained_cuts"][0]
    assert retained["enzymes"] == ["EnzB"]
    assert [s["enzyme"] for s in retained["sites"]] == ["EnzB"]
    assert [f["length"] for f in data["fragments"]] == [4, 4]


def test_shared_physical_cut_removed_only_when_all_sites_blocked(client):
    body = _create(client, SHARED_CUT_SEQ, enzymes=[ENZ_A, ENZ_B])
    # 坐标 2（偏移 0）阻断 EnzA、坐标 5（偏移 3）阻断 EnzB
    resp = _analyze(
        client, body["id"], [2, 5],
        [
            {"enzyme": "EnzA", "sensitive_offsets": [0]},
            {"enzyme": "EnzB", "sensitive_offsets": [3]},
        ],
    )
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 2
    assert len(data["retained_cuts"]) == 0
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["removed_cuts"]] == [(4, 4)]
    assert [f["length"] for f in data["fragments"]] == [8]


# ---------- 反向位点与环形坐标 ----------


def test_reverse_strand_site_uses_reverse_offset_mapping(client):
    # Fiver 模式 GGATC（非回文）：序列 AAGATCCAAA 上反向命中 GATCC，起于 2，
    # 切口 (3, 6)；反向位点偏移 0 -> start + k - 1 - 0 = 6
    body = _create(client, "AAGATCCAAA", enzymes=[FIVER])
    site = body["cuts"][0]["sites"][0]
    assert site["strand"] == "-" and site["start"] == 2
    assert (body["cuts"][0]["top_cut"], body["cuts"][0]["bottom_cut"]) == (3, 6)

    resp = _analyze(client, body["id"], [6], [{"enzyme": "Fiver", "sensitive_offsets": [0]}])
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["site"]["strand"] == "-"
    assert data["blocked_sites"][0]["methylated_positions"] == [6]
    assert len(data["removed_cuts"]) == 1

    # 正向映射坐标 start + offset = 2 不是反向位点的敏感位置
    resp2 = _analyze(client, body["id"], [2], [{"enzyme": "Fiver", "sensitive_offsets": [0]}])
    assert resp2.status_code == 200
    assert resp2.json()["blocked_sites"] == []
    assert len(resp2.json()["retained_cuts"]) == 1


def test_circular_site_wraps_offset_with_modulo(client):
    # 环形 ATTCGGGGGA（长度 10）：GAATTC 跨原点起于 8；
    # 偏移 5 的物理坐标 (8 + 5) % 10 = 3
    body = _create(client, "ATTCGGGGGA", topology="circular", enzymes=[ECORI])
    assert body["cuts"][0]["sites"][0]["start"] == 8

    resp = _analyze(client, body["id"], [3], [{"enzyme": "EcoRI", "sensitive_offsets": [5]}])
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["blocked_sites"]) == 1
    assert data["blocked_sites"][0]["methylated_positions"] == [3]
    # 环形唯一切口被移除 -> 整条环状分子
    assert data["removed_cuts"] and not data["retained_cuts"]
    assert [f["length"] for f in data["fragments"]] == [10]


def test_circular_single_remaining_cut_rebuilds_full_length_fragment(client):
    # 环形 CGCAAAG（长度 7）：GC 位点起于 1（切口 (2,2)）与起于 6 跨原点
    # （切口 (0,0)）；阻断前者（偏移 0 -> 坐标 1）
    body = _create(client, "CGCAAAG", topology="circular", enzymes=[GC_ENZYME])
    assert {(c["top_cut"], c["bottom_cut"]) for c in body["cuts"]} == {(0, 0), (2, 2)}
    site_starts = sorted(s["start"] for c in body["cuts"] for s in c["sites"])
    assert site_starts == [1, 6]

    resp = _analyze(client, body["id"], [1], [{"enzyme": "GC-cutter", "sensitive_offsets": [0]}])
    assert resp.status_code == 200
    data = resp.json()
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["retained_cuts"]] == [(0, 0)]
    assert [(c["top_cut"], c["bottom_cut"]) for c in data["removed_cuts"]] == [(2, 2)]
    # 环形单切口开环：1 条全长片段
    assert len(data["fragments"]) == 1
    assert data["fragments"][0]["length"] == 7


# ---------- 稳定排序与顺序无关 ----------


def test_blocked_sites_are_stably_sorted(client):
    # 两个 EcoRI 位点全部阻断
    body = _create(client, TWO_SITE_SEQ)
    resp = _analyze(
        client, body["id"], [0, 10],
        [{"enzyme": "EcoRI", "sensitive_offsets": [0]}],
    )
    assert resp.status_code == 200
    keys = [
        (b["top_cut"], b["bottom_cut"]) for b in resp.json()["blocked_sites"]
    ]
    assert keys == sorted(keys) == [(1, 5), (11, 15)]
    removed = [
        (c["top_cut"], c["bottom_cut"]) for c in resp.json()["removed_cuts"]
    ]
    assert removed == [(1, 5), (11, 15)]


def test_result_independent_of_request_entry_order(client):
    body = _create(client, TWO_SITE_SEQ)
    # 乱序给出甲基化坐标与敏感偏移（含重复偏移），回显与结果均须确定
    rule = {"enzyme": "EcoRI", "sensitive_offsets": [3, 0, 0]}
    r1 = _analyze(client, body["id"], [10, 0], [rule]).json()
    r2 = _analyze(client, body["id"], [0, 10], [
        {"enzyme": "EcoRI", "sensitive_offsets": [0, 3, 0]}
    ]).json()
    assert r1 == r2
    assert r1["methylated_positions"] == [0, 10]


# ---------- 校验与状态码 ----------


def test_missing_record_returns_404(client):
    resp = _analyze(
        client, 999, [0], [{"enzyme": "EcoRI", "sensitive_offsets": [0]}]
    )
    assert resp.status_code == 404


def test_duplicate_methylated_positions_return_422(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(
        client, body["id"], [4, 4],
        [{"enzyme": "EcoRI", "sensitive_offsets": [0]}],
    )
    assert resp.status_code == 422
    assert "duplicate methylated positions" in resp.text


def test_out_of_bounds_positions_return_422(client):
    body = _create(client, ECO_SEQ)  # 长度 14
    for bad in (14, 100):
        resp = _analyze(
            client, body["id"], [bad],
            [{"enzyme": "EcoRI", "sensitive_offsets": [0]}],
        )
        assert resp.status_code == 422, resp.text
        assert str(bad) in resp.json()["detail"]


def test_negative_position_returns_422(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(
        client, body["id"], [-1],
        [{"enzyme": "EcoRI", "sensitive_offsets": [0]}],
    )
    assert resp.status_code == 422


def test_unknown_enzyme_name_returns_422(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(
        client, body["id"], [4],
        [{"enzyme": "BamHI", "sensitive_offsets": [0]}],
    )
    assert resp.status_code == 422
    assert "unknown enzyme names" in resp.json()["detail"]
    assert "BamHI" in resp.json()["detail"]


def test_duplicate_enzyme_rules_return_422(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(
        client, body["id"], [4],
        [
            {"enzyme": "EcoRI", "sensitive_offsets": [0]},
            {"enzyme": "EcoRI", "sensitive_offsets": [1]},
        ],
    )
    assert resp.status_code == 422
    assert "duplicate enzyme rules" in resp.text


def test_offset_out_of_range_returns_422(client):
    body = _create(client, ECO_SEQ)  # EcoRI 模式长 6，最大偏移 5
    for bad in (6, 100):
        resp = _analyze(
            client, body["id"], [4],
            [{"enzyme": "EcoRI", "sensitive_offsets": [bad]}],
        )
        assert resp.status_code == 422, resp.text
        detail = resp.json()["detail"]
        assert "EcoRI" in detail and str(bad) in detail


def test_negative_offset_returns_422(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(
        client, body["id"], [4],
        [{"enzyme": "EcoRI", "sensitive_offsets": [-1]}],
    )
    assert resp.status_code == 422


def test_empty_positions_and_rules_return_422(client):
    body = _create(client, ECO_SEQ)
    assert client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={"methylated_positions": [], "enzyme_rules": [
            {"enzyme": "EcoRI", "sensitive_offsets": [0]}
        ]},
    ).status_code == 422
    assert client.post(
        f"/digests/{body['id']}/methylation-analysis",
        json={"methylated_positions": [4], "enzyme_rules": []},
    ).status_code == 422


def test_unknown_enzyme_and_bad_position_reported_together(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(
        client, body["id"], [4, 99],
        [
            {"enzyme": "EcoRI", "sensitive_offsets": [7]},
            {"enzyme": "BamHI", "sensitive_offsets": [0]},
        ],
    )
    assert resp.status_code == 422
    detail = resp.json()["detail"]
    assert "99" in detail and "BamHI" in detail and "EcoRI" in detail


def test_duplicate_offsets_within_one_rule_are_allowed(client):
    body = _create(client, ECO_SEQ)
    resp = _analyze(
        client, body["id"], [4],
        [{"enzyme": "EcoRI", "sensitive_offsets": [0, 0, 0]}],
    )
    assert resp.status_code == 200
    assert len(resp.json()["blocked_sites"]) == 1


# ---------- 只读 ----------


def test_analysis_is_read_only(client):
    body = _create(client, TWO_SITE_SEQ)
    resp = _analyze(
        client, body["id"], [0],
        [{"enzyme": "EcoRI", "sensitive_offsets": [0]}],
    )
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{body['id']}").json() == body
    assert len(client.get("/digests").json()) == 1
