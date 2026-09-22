"""API 层测试：提交、查询、校验与原子拒绝。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


def eco_request(**overrides):
    payload = {
        "sequence": "AAAAGAATTCTTTT",
        "topology": "linear",
        "enzymes": [{"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}],
    }
    payload.update(overrides)
    return payload


def test_create_and_get_digest(client):
    resp = client.post("/digests", json=eco_request())
    assert resp.status_code == 201
    body = resp.json()
    assert body["id"] == 1
    assert body["topology"] == "linear"
    assert body["sequence_length"] == 14
    assert len(body["cuts"]) == 1
    assert body["cuts"][0]["end_type"] == "5prime"
    assert body["cuts"][0]["overhang"] == "AATT"
    assert [f["length"] for f in body["fragments"]] == [5, 9]

    got = client.get(f"/digests/{body['id']}")
    assert got.status_code == 200
    assert got.json() == body  # 查询结果与提交时完全一致


def test_list_digests(client):
    client.post("/digests", json=eco_request())
    client.post("/digests", json=eco_request(sequence="GCGCGCGC", topology="circular"))
    resp = client.get("/digests")
    assert resp.status_code == 200
    items = resp.json()
    assert len(items) == 2
    assert items[0]["id"] == 1 and items[0]["cut_count"] == 1
    assert items[1]["id"] == 2 and items[1]["cut_count"] == 0
    assert items[1]["fragment_count"] == 1  # 无切口保留整条分子


def test_get_missing_digest_returns_404(client):
    assert client.get("/digests/999").status_code == 404


def test_illegal_offset_rejects_entire_save(client):
    bad = eco_request(
        enzymes=[
            {"name": "Good", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
            {"name": "Bad", "pattern": "GATATC", "cut_top": 7, "cut_bottom": 3},
        ]
    )
    resp = client.post("/digests", json=bad)
    assert resp.status_code == 422
    assert client.get("/digests").json() == []  # 整次保存被拒绝，无任何落库


def test_negative_offset_rejected(client):
    bad = eco_request(
        enzymes=[{"name": "Neg", "pattern": "GAATTC", "cut_top": -1, "cut_bottom": 5}]
    )
    assert client.post("/digests", json=bad).status_code == 422
    assert client.get("/digests").json() == []


def test_invalid_sequence_base_rejected(client):
    resp = client.post("/digests", json=eco_request(sequence="ACGTXACGT"))
    assert resp.status_code == 422
    assert client.get("/digests").json() == []


def test_invalid_pattern_base_rejected(client):
    bad = eco_request(
        enzymes=[{"name": "Bad", "pattern": "GAUZ", "cut_top": 1, "cut_bottom": 1}]
    )
    assert client.post("/digests", json=bad).status_code == 422


def test_empty_enzyme_list_rejected(client):
    assert client.post("/digests", json=eco_request(enzymes=[])).status_code == 422


def test_circular_digest_via_api(client):
    payload = eco_request(sequence="ATTCGGGGGA", topology="circular")
    resp = client.post("/digests", json=payload)
    assert resp.status_code == 201
    body = resp.json()
    assert body["cuts"][0]["sites"][0]["start"] == 8  # 跨原点命中
    assert len(body["fragments"]) == 1
    assert body["fragments"][0]["length"] == 10


def test_same_input_same_result_order(client):
    payload = eco_request(
        sequence="GAATTCAAAAGATATCTTTTGGTACC",
        topology="circular",
        enzymes=[
            {"name": "KpnI", "pattern": "GGTACC", "cut_top": 5, "cut_bottom": 1},
            {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5},
            {"name": "EcoRV", "pattern": "GATATC", "cut_top": 3, "cut_bottom": 3},
        ],
    )
    r1 = client.post("/digests", json=payload).json()
    r2 = client.post("/digests", json=payload).json()
    assert r1["cuts"] == r2["cuts"]
    assert r1["fragments"] == r2["fragments"]
