"""凝胶条带比较接口测试。"""

import pytest
from fastapi.testclient import TestClient

from app.main import app


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    with TestClient(app) as c:
        yield c


# 回文模式 GC 的一位切口酶（cut 1/1 平端），只在 GC 位点产生切口，
# 便于按目标片段长度精确构造记录
GC_ENZYME = {"name": "GC-cutter", "pattern": "GC", "cut_top": 1, "cut_bottom": 1}
# 识别 GAATTC，在全 A 序列上无位点，用于构造"不切割"的记录
ECORI = {"name": "EcoRI", "pattern": "GAATTC", "cut_top": 1, "cut_bottom": 5}


def _seq_with_cuts(length, cut_positions):
    """构造线性序列：全 A 背景，在每个切口位置 p 处放置 GC（切口在 p-1 与 p 之间）。

    相邻切口至少相隔 2（即片段长度 >= 2），否则 GC 放置会互相覆盖。
    """
    seq = ["A"] * length
    for p in cut_positions:
        seq[p - 1] = "G"
        seq[p] = "C"
    return "".join(seq)


def _create(client, sequence, enzymes=None):
    resp = client.post(
        "/digests",
        json={
            "sequence": sequence,
            "topology": "linear",
            "enzymes": enzymes if enzymes is not None else [GC_ENZYME],
        },
    )
    assert resp.status_code == 201
    return resp.json()


def _create_with_fragments(client, fragment_lengths):
    """创建一条线性记录，使其片段长度恰为给定清单（按物理顺序，每段 >= 2）。"""
    cuts = []
    position = 0
    for frag_len in fragment_lengths[:-1]:
        position += frag_len
        cuts.append(position)
    return _create(client, _seq_with_cuts(sum(fragment_lengths), cuts))


def _create_uncut(client, length):
    """创建一条无切口的记录，唯一片段为全长分子。"""
    return _create(client, "A" * length, enzymes=[ECORI])


def _compare(client, digest_ids, ratio):
    return client.post(
        "/gel-comparisons",
        json={"digest_ids": digest_ids, "resolution_ratio": ratio},
    )


# ---------- 聚类与条带统计 ----------


def test_bands_clustered_with_lane_counts_and_diagnostic(client):
    lane_a = _create_with_fragments(client, [40, 60])      # 片段 [40, 60]
    lane_b = _create_with_fragments(client, [41, 20, 21])  # 片段 [41, 20, 21]
    resp = _compare(client, [lane_a["id"], lane_b["id"]], 0.05)
    assert resp.status_code == 200
    data = resp.json()
    assert data["digest_ids"] == [lane_a["id"], lane_b["id"]]
    assert data["resolution_ratio"] == 0.05
    assert data["lane_count"] == 2
    assert data["band_count"] == 3
    # 条带按 max_length 降序稳定排列
    assert data["bands"] == [
        {"min_length": 60, "max_length": 60, "lane_fragment_counts": [1, 0],
         "fragment_count": 1, "total_bases": 60, "diagnostic": True},
        {"min_length": 40, "max_length": 41, "lane_fragment_counts": [1, 1],
         "fragment_count": 2, "total_bases": 81, "diagnostic": False},
        {"min_length": 20, "max_length": 21, "lane_fragment_counts": [0, 2],
         "fragment_count": 2, "total_bases": 41, "diagnostic": True},
    ]


def test_clustering_anchors_to_band_max_not_neighbour(client):
    # 片段 [100, 91, 82]，ratio 0.1：91 相对 100 的跨度 0.09 可并入；
    # 82 相对条带最大值 100 的跨度 0.18 超差，必须新开条带——禁止因
    # (91-82)/91≈0.099 相邻可合并而把 82 链式归入 [100, 91] 所在条带
    lane_a = _create_with_fragments(client, [100, 91, 82])
    lane_b = _create_uncut(client, 50)
    resp = _compare(client, [lane_a["id"], lane_b["id"]], 0.1)
    assert resp.status_code == 200
    bands = resp.json()["bands"]
    assert [(b["min_length"], b["max_length"]) for b in bands] == [
        (91, 100), (82, 82), (50, 50),
    ]
    assert bands[0]["lane_fragment_counts"] == [2, 0]
    assert bands[0]["fragment_count"] == 2
    assert bands[0]["total_bases"] == 191
    assert bands[1]["lane_fragment_counts"] == [1, 0]
    assert bands[2]["lane_fragment_counts"] == [0, 1]
    assert all(b["diagnostic"] for b in bands)


def test_zero_ratio_merges_only_identical_lengths(client):
    lane_a = _create_with_fragments(client, [40, 60])
    lane_b = _create_with_fragments(client, [41, 59])
    resp = _compare(client, [lane_a["id"], lane_b["id"]], 0.0)
    assert resp.status_code == 200
    bands = resp.json()["bands"]
    # ratio=0 时 41 与 40、59 与 60 都不能合并，各自成带
    assert [(b["min_length"], b["max_length"]) for b in bands] == [
        (60, 60), (59, 59), (41, 41), (40, 40),
    ]
    assert [b["lane_fragment_counts"] for b in bands] == [
        [1, 0], [0, 1], [0, 1], [1, 0],
    ]
    assert all(b["diagnostic"] for b in bands)

    # 完全相同的长度在 ratio=0 时仍并入同一条带
    lane_c = _create_with_fragments(client, [40, 60])
    resp = _compare(client, [lane_a["id"], lane_c["id"]], 0.0)
    bands = resp.json()["bands"]
    assert [(b["min_length"], b["max_length"]) for b in bands] == [(60, 60), (40, 40)]
    assert [b["lane_fragment_counts"] for b in bands] == [[1, 1], [1, 1]]
    assert all(not b["diagnostic"] for b in bands)


def test_resolution_ratio_upper_boundary_is_inclusive(client):
    lane_a = _create_uncut(client, 100)  # 片段 [100]
    lane_b = _create_uncut(client, 80)   # 片段 [80]
    # (100-80)/100 = 0.2，恰好等于上限，仍应并入同一条带
    resp = _compare(client, [lane_a["id"], lane_b["id"]], 0.2)
    assert resp.status_code == 200
    bands = resp.json()["bands"]
    assert bands == [
        {"min_length": 80, "max_length": 100, "lane_fragment_counts": [1, 1],
         "fragment_count": 2, "total_bases": 180, "diagnostic": False}
    ]

    # 比例略小于 0.2 则拆为两条
    resp = _compare(client, [lane_a["id"], lane_b["id"]], 0.19)
    assert [(b["min_length"], b["max_length"]) for b in resp.json()["bands"]] == [
        (100, 100), (80, 80),
    ]


def test_lane_counts_follow_digest_ids_order(client):
    lane_a = _create_with_fragments(client, [40, 60])  # 片段 [40, 60]
    lane_b = _create_uncut(client, 40)                 # 片段 [40]
    forward = _compare(client, [lane_a["id"], lane_b["id"]], 0.0).json()
    reverse = _compare(client, [lane_b["id"], lane_a["id"]], 0.0).json()
    assert forward["digest_ids"] == [lane_a["id"], lane_b["id"]]
    assert reverse["digest_ids"] == [lane_b["id"], lane_a["id"]]
    # lane_fragment_counts 与请求 digest_ids 的顺序一一对应
    assert [b["lane_fragment_counts"] for b in forward["bands"]] == [[1, 0], [1, 1]]
    assert [b["lane_fragment_counts"] for b in reverse["bands"]] == [[0, 1], [1, 1]]
    # 条带范围与总碱基数本身与 id 顺序无关
    assert [
        (b["min_length"], b["max_length"], b["total_bases"]) for b in forward["bands"]
    ] == [
        (b["min_length"], b["max_length"], b["total_bases"]) for b in reverse["bands"]
    ]


def test_eight_lanes_supported(client):
    lanes = [_create_uncut(client, 50) for _ in range(8)]
    resp = _compare(client, [lane["id"] for lane in lanes], 0.1)
    assert resp.status_code == 200
    data = resp.json()
    assert data["lane_count"] == 8
    assert data["bands"] == [
        {"min_length": 50, "max_length": 50, "lane_fragment_counts": [1] * 8,
         "fragment_count": 8, "total_bases": 400, "diagnostic": False}
    ]


# ---------- 校验与状态码 ----------


def test_duplicate_digest_ids_return_422(client):
    lane = _create_uncut(client, 50)
    resp = _compare(client, [lane["id"], lane["id"]], 0.1)
    assert resp.status_code == 422
    assert "duplicate digest ids" in resp.text


def test_ratio_out_of_range_returns_422(client):
    lane_a = _create_uncut(client, 50)
    lane_b = _create_uncut(client, 60)
    for bad in (-0.01, 0.21, 1.0):
        resp = _compare(client, [lane_a["id"], lane_b["id"]], bad)
        assert resp.status_code == 422, (bad, resp.text)


def test_lane_count_limits_return_422(client):
    lane = _create_uncut(client, 50)
    # 少于 2 条
    assert _compare(client, [lane["id"]], 0.1).status_code == 422
    # 多于 8 条（校验先于记录存在性检查）
    assert _compare(client, list(range(1, 10)), 0.1).status_code == 422


def test_missing_record_returns_404(client):
    lane = _create_uncut(client, 50)
    resp = _compare(client, [lane["id"], 999], 0.1)
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_first_missing_id_in_request_order_is_reported(client):
    resp = _compare(client, [999, 998], 0.1)
    assert resp.status_code == 404
    assert "999" in resp.json()["detail"]


def test_validation_precedes_existence_check(client):
    # 重复 id 直接 422，即使 id 不存在也不报 404
    assert _compare(client, [999, 999], 0.1).status_code == 422


# ---------- 只读 ----------


def test_comparison_is_read_only(client):
    lane_a = _create_with_fragments(client, [40, 60])
    lane_b = _create_with_fragments(client, [41, 20, 21])
    resp = _compare(client, [lane_a["id"], lane_b["id"]], 0.05)
    assert resp.status_code == 200
    # 原记录内容不变，且没有新增任何提交
    assert client.get(f"/digests/{lane_a['id']}").json() == lane_a
    assert client.get(f"/digests/{lane_b['id']}").json() == lane_b
    assert len(client.get("/digests").json()) == 2
