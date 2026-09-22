"""变异影响分析接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


def _eco_payload(sequence="AAAAGAATTCTTTT", topology="linear"):
    return {
        "sequence": sequence,
        "topology": topology,
        "enzymes": [
            {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}
        ],
    }


def _create(client, **overrides):
    payload = _eco_payload(**overrides)
    resp = client.post("/digests", json=payload)
    assert resp.status_code == 201
    return resp.json()


def _analyze(client, digest_id, variants):
    return client.post(
        f"/digests/{digest_id}/variant-analysis", json={"variants": variants}
    )


def test_variant_destroys_site_cut_removed_and_fragment_lengths_change(client):
    body = _create(client)  # EcoRI 位点 GAATTC 起于 4，G->C 即破坏
    resp = _analyze(client, body["id"], [{"position": 4, "alt": "C"}])
    assert resp.status_code == 200
    data = resp.json()

    assert data["variants"] == [{"position": 4, "ref": "G", "alt": "C"}]
    assert data["variant_sequence"] == "AAAACAATTCTTTT"
    assert data["sequence_length"] == 14

    cuts = data["cuts"]
    assert len(cuts["removed"]) == 1
    assert (cuts["removed"][0]["top_cut"], cuts["removed"][0]["bottom_cut"]) == (5, 9)
    assert cuts["added"] == []
    assert cuts["changed"] == []

    frags = data["fragments"]
    assert frags["before_lengths"] == [5, 9]
    assert frags["after_lengths"] == [14]
    assert len(frags["removed"]) == 1 and frags["removed"][0]["start"] == 5
    # 起点 0 的片段仍在，长度由 5 变为 14
    assert frags["changed"] == [
        {"start": 0, "before_length": 5, "after_length": 14, "length_delta": 9}
    ]


def test_variant_creates_site_cut_added_and_fragments_split(client):
    body = _create(client, sequence="AAAAGAAATCTTTT")  # 与 GAATTC 仅差一个 T->T? 构造潜在位点
    # GAAATC：将位置 8（索引：AAAA=0-3, G=4, A=5,A=6,A=7,T=8,C=9）的 A? 实际序列:
    # A A A A G A A A T C T T T T
    # 0 1 2 3 4 5 6 7 8 9 ...
    # GAATTC 需在 4..9 为 GAATTC：把 8 位 T? 现有 8=T，9=C；需要 7=T（现 7=A）
    resp = _analyze(client, body["id"], [{"position": 7, "alt": "T"}])
    assert resp.status_code == 200
    data = resp.json()
    assert data["variant_sequence"] == "AAAAGAATTCTTTT"
    assert len(data["cuts"]["added"]) == 1
    assert (data["cuts"]["added"][0]["top_cut"],
            data["cuts"]["added"][0]["bottom_cut"]) == (5, 9)
    assert data["cuts"]["removed"] == []
    frags = data["fragments"]
    assert frags["before_lengths"] == [14]
    assert frags["after_lengths"] == [5, 9]
    assert frags["added"][0]["start"] == 5
    assert frags["changed"] == [
        {"start": 0, "before_length": 14, "after_length": 5, "length_delta": -9}
    ]


def test_no_change_for_same_base_substitution(client):
    body = _create(client)
    resp = _analyze(client, body["id"], [{"position": 4, "alt": "g"}])  # 小写归一化为 G
    assert resp.status_code == 200
    data = resp.json()
    assert data["variant_sequence"] == body["sequence"]
    assert data["cuts"] == {"added": [], "removed": [], "changed": []}
    assert data["fragments"]["before_lengths"] == data["fragments"]["after_lengths"]
    assert data["fragments"]["added"] == []
    assert data["fragments"]["removed"] == []
    assert data["fragments"]["changed"] == []


def test_multiple_variants_echo_sorted_by_position(client):
    body = _create(client)
    # 故意乱序传入，回显必须按位置稳定排序
    resp = _analyze(
        client,
        body["id"],
        [{"position": 9, "alt": "A"}, {"position": 0, "alt": "T"},
         {"position": 4, "alt": "T"}],
    )
    assert resp.status_code == 200
    positions = [v["position"] for v in resp.json()["variants"]]
    assert positions == [0, 4, 9]


def test_degenerate_variant_keeps_site_but_changes_overhang(client):
    # HinfI 识别 GANTC（简并），位点 GAGTC；模式 N 与序列 G 匹配。
    # 将 overhang 区段内碱基替换为简并碱基不破坏匹配，仅改变 overhang 序列。
    payload = {
        "sequence": "AAAGAGTCAAA",
        "topology": "linear",
        "enzymes": [
            {"name": "HinfI", "pattern": "GANTC", "cut_top": 1, "cut_bottom": 4}
        ],
    }
    created = client.post("/digests", json=payload).json()
    # 切口 (4,7)：overhang = seq[4:7] = "AGT"
    assert created["cuts"][0]["overhang"] == "AGT"
    # 位置 6 的 T 替换为 W(A/T)：N∩W 仍有交集，位点存活，overhang 变 AGW
    resp = _analyze(client, created["id"], [{"position": 6, "alt": "W"}])
    assert resp.status_code == 200
    data = resp.json()
    assert data["cuts"]["added"] == [] and data["cuts"]["removed"] == []
    changed = data["cuts"]["changed"]
    assert len(changed) == 1
    assert changed[0]["before"]["overhang"] == "AGT"
    assert changed[0]["after"]["overhang"] == "AGW"
    assert changed[0]["changes"] == ["overhang", "sites"]
    assert changed[0]["after"]["sites"][0]["matched"] == "GAGWC"
    # 坐标未变、片段长度不变
    assert data["fragments"]["changed"] == []
    assert data["fragments"]["before_lengths"] == data["fragments"]["after_lengths"]


def test_circular_variant_analysis(client):
    # GAATTC 跨原点（位点起于 8），破坏其中一个碱基
    body = _create(client, sequence="ATTCGGGGGA", topology="circular")
    assert body["cuts"][0]["sites"][0]["start"] == 8
    resp = _analyze(client, body["id"], [{"position": 8, "alt": "T"}])
    assert resp.status_code == 200
    data = resp.json()
    assert len(data["cuts"]["removed"]) == 1
    assert data["cuts"]["removed"][0]["top_cut"] == 9
    # 环形无切口：一条全长片段
    assert data["fragments"]["after_lengths"] == [10]


def test_position_out_of_bounds_returns_422(client):
    body = _create(client)  # 长度 14
    for bad_pos in (14, 100):
        resp = _analyze(client, body["id"], [{"position": bad_pos, "alt": "C"}])
        assert resp.status_code == 422, resp.text


def test_illegal_base_returns_422(client):
    body = _create(client)
    for bad_alt in ("X", "Z", "U", ""):
        resp = _analyze(client, body["id"], [{"position": 4, "alt": bad_alt}])
        assert resp.status_code == 422, bad_alt


def test_duplicate_positions_returns_422(client):
    body = _create(client)
    resp = _analyze(
        client,
        body["id"],
        [{"position": 4, "alt": "C"}, {"position": 4, "alt": "T"}],
    )
    assert resp.status_code == 422
    assert "duplicate" in resp.text.lower()


def test_negative_position_returns_422(client):
    body = _create(client)
    resp = _analyze(client, body["id"], [{"position": -1, "alt": "C"}])
    assert resp.status_code == 422


def test_empty_variant_list_returns_422(client):
    body = _create(client)
    resp = client.post(
        f"/digests/{body['id']}/variant-analysis", json={"variants": []}
    )
    assert resp.status_code == 422


def test_missing_record_returns_404(client):
    resp = _analyze(client, 999, [{"position": 0, "alt": "C"}])
    assert resp.status_code == 404


def test_variant_analysis_does_not_persist(client):
    body = _create(client)
    _analyze(client, body["id"], [{"position": 4, "alt": "C"}])
    # 原记录内容不变，且没有新增任何提交
    stored = client.get(f"/digests/{body['id']}").json()
    assert stored == body
    assert len(client.get("/digests").json()) == 1


def test_result_order_is_stable(client):
    body = _create(
        client,
        sequence="GAATTC" + "AAAA" + "GAATTC" + "GGGG",
        topology="circular",
    )
    # 同时破坏两个位点中的各一个碱基
    variants = [{"position": 1, "alt": "C"}, {"position": 11, "alt": "C"}]
    r1 = _analyze(client, body["id"], variants).json()
    r2 = _analyze(client, body["id"], list(reversed(variants))).json()
    assert r1 == r2  # 与变异传入顺序无关
    removed = [(c["top_cut"], c["bottom_cut"]) for c in r1["cuts"]["removed"]]
    assert removed == sorted(removed)
