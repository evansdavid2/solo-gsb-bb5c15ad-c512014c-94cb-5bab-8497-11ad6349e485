"""FastAPI 应用入口。

启动时（lifespan）自动完成数据库初始化；接口：
    POST /digests                          提交序列与酶列表，执行酶切并保存，返回完整结果
    GET  /digests                          列出全部提交的摘要
    GET  /digests/{id}                     按 id 查询保存的完整酶切结果
    POST /digests/{id}/variant-analysis    基于已保存记录施加单碱基替换，重跑酶切并对比差异（不落库）
    POST /digests/{id}/enzyme-panel-plan   枚举候选酶组合，寻找达到目标片段数的最少酶集合（只读，不落库）
    POST /digests/{id}/partial-digest-distribution
                                           为每个已保存切口指定独立切割概率，枚举成功/失败
                                           状态并返回片段长度多重集的概率分布（只读，不落库）
    POST /digests/{id}/methylation-analysis 只读分析甲基化对已保存酶切位点的阻断影响，
                                           重算保留切口并重建片段（不落库）
    POST /gel-comparisons                只读比较 2~8 条已保存记录的凝胶条带：汇总片段
                                           长度按降序聚类，返回条带范围、各泳道片段数、
                                           总碱基数与诊断条带（不落库）
    POST /ligation-compatibility         只读枚举 2~8 条已保存记录之间可连接的片段末端对：
                                           黏端要求一方 5' 突出端序列等于另一方 3' 突出端
                                           序列的反向互补，平端仅与平端兼容（include_blunt
                                           可关闭），返回稳定排序的配对、兼容类型、
                                           连接后突出端与总数（不落库）
    POST /ligation-assembly              只读规划 2~8 条已保存记录间的最长片段拼接链：
                                           复用末端反向互补兼容判定构建有向图，搜索不重复
                                           使用记录的最长线性链（每条记录最多贡献一个
                                           片段，可正向/反向摆入），同长度方案按请求顺序、
                                           片段索引与端方向稳定择优；include_blunt 可关闭
                                           平端连接，max_chain_length 限制链上节点数；
                                           返回节点端点、连接类型/突出端、总连接数与
                                           未能延伸的原因（不落库）
    GET  /health                           健康检查
"""

from __future__ import annotations

import os
from contextlib import asynccontextmanager

from fastapi import FastAPI, HTTPException

from .db import Database
from .digest import run_digest
from .assembly import MAX_ASSEMBLY_FRAGMENTS, plan_assembly_chain, public_node
from .distribution import MAX_PARTIAL_CUTS, partial_digest_distribution
from .gel import compare_gel_lanes
from .ligation import find_compatible_pairs
from .methylation import analyze_methylation
from .panel import plan_enzyme_panel
from .schemas import (
    DigestRequest,
    DigestResponse,
    DigestSummary,
    EnzymeIn,
    EnzymePanelPlanRequest,
    EnzymePanelPlanResponse,
    GelComparisonRequest,
    GelComparisonResponse,
    LigationAssemblyRequest,
    LigationAssemblyResponse,
    LigationCompatibilityRequest,
    LigationCompatibilityResponse,
    MethylationAnalysisRequest,
    MethylationAnalysisResponse,
    PartialDigestDistributionResponse,
    PartialDigestRequest,
    VariantAnalysisRequest,
    VariantAnalysisResponse,
)
from .variant import apply_substitutions, diff_cuts, diff_fragments, variant_summary

db: Database


@asynccontextmanager
async def lifespan(app: FastAPI):
    global db
    db = Database(os.environ.get("DATABASE_PATH", "data/digest.db"))
    yield


app = FastAPI(
    title="DNA Restriction Digest Simulator",
    version="1.0.0",
    description="DNA 限制性酶切模拟：支持线性/环形序列、IUPAC 简并碱基、"
                "双链切口偏移、黏端/平端分析与片段生成。",
    lifespan=lifespan,
)


def _build_response(record: dict) -> dict:
    return {
        "id": record["id"],
        "topology": record["topology"],
        "sequence": record["sequence"],
        "sequence_length": len(record["sequence"]),
        "enzymes": record["enzymes"],
        "cuts": record["result"]["cuts"],
        "fragments": record["result"]["fragments"],
        "created_at": record["created_at"],
    }


@app.get("/health")
def health() -> dict:
    return {"status": "ok"}


@app.post("/digests", status_code=201, response_model=DigestResponse)
def create_digest(request: DigestRequest) -> dict:
    # Pydantic 已在此前完成全部校验（含切口偏移合法性），
    # 非法请求在写库之前即以 422 拒绝，保证整次保存被拒绝。
    result = run_digest(request.sequence, request.topology, request.enzymes)
    submission_id, created_at = db.save_submission(
        request.sequence, request.topology, request.enzymes, result
    )
    return _build_response(
        {
            "id": submission_id,
            "sequence": request.sequence,
            "topology": request.topology,
            "created_at": created_at,
            "enzymes": [e.model_dump() for e in request.enzymes],
            "result": result,
        }
    )


@app.get("/digests", response_model=list[DigestSummary])
def list_digests() -> list[dict]:
    return db.list_submissions()


@app.get("/digests/{submission_id}", response_model=DigestResponse)
def get_digest(submission_id: int) -> dict:
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")
    return _build_response(record)


@app.post(
    "/digests/{submission_id}/variant-analysis",
    response_model=VariantAnalysisResponse,
)
def analyze_variants(submission_id: int, request: VariantAnalysisRequest) -> dict:
    # 先确认记录存在：不存在返回 404，且不泄露任何校验细节
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")

    sequence = record["sequence"]
    length = len(sequence)
    # 位置越界依赖记录的序列长度，只能在取到记录后校验 -> 422
    out_of_bounds = sorted(
        v.position for v in request.variants if v.position >= length
    )
    if out_of_bounds:
        raise HTTPException(
            status_code=422,
            detail=(
                f"variant positions out of bounds for sequence length {length}: "
                + ", ".join(str(p) for p in out_of_bounds)
            ),
        )

    # 仅在内存中施加替换并重跑酶切，原记录与数据库均不改动
    enzymes = [EnzymeIn(**e) for e in record["enzymes"]]
    variant_sequence = apply_substitutions(sequence, request.variants)
    before = record["result"]
    after = run_digest(variant_sequence, record["topology"], enzymes)

    return {
        "id": record["id"],
        "sequence_length": length,
        "variant_sequence": variant_sequence,
        "variants": variant_summary(sequence, request.variants),
        "cuts": diff_cuts(before["cuts"], after["cuts"]),
        "fragments": diff_fragments(before["fragments"], after["fragments"]),
    }


@app.post(
    "/digests/{submission_id}/enzyme-panel-plan",
    response_model=EnzymePanelPlanResponse,
)
def plan_panel(submission_id: int, request: EnzymePanelPlanRequest) -> dict:
    # 先确认记录存在：不存在返回 404，且不泄露任何校验细节
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")

    # 仅在内存中枚举候选组合并复用酶切计算，原记录与数据库均不改动
    plan = plan_enzyme_panel(
        record["sequence"],
        record["topology"],
        request.candidates,
        request.target_fragment_count,
    )
    if plan is None:
        raise HTTPException(
            status_code=422,
            detail=(
                f"no combination of the {len(request.candidates)} candidate "
                f"enzyme(s) yields exactly {request.target_fragment_count} "
                "fragment(s)"
            ),
        )

    return {
        "id": record["id"],
        "target_fragment_count": request.target_fragment_count,
        "selected_enzymes": [e.model_dump() for e in plan["enzymes"]],
        "fragment_count": plan["fragment_count"],
        "cuts": plan["cuts"],
    }


@app.post(
    "/digests/{submission_id}/partial-digest-distribution",
    response_model=PartialDigestDistributionResponse,
)
def partial_digest(submission_id: int, request: PartialDigestRequest) -> dict:
    # 先确认记录存在：不存在返回 404，且不泄露任何校验细节
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")

    saved_cuts = record["result"]["cuts"]

    # 枚举规模有界：超过 12 个物理切口无法在保证耗时的前提下全枚举 -> 422
    if len(saved_cuts) > MAX_PARTIAL_CUTS:
        raise HTTPException(
            status_code=422,
            detail=(
                f"submission {submission_id} has {len(saved_cuts)} saved cuts; "
                f"partial-digest enumeration supports at most {MAX_PARTIAL_CUTS}"
            ),
        )

    saved_keys = {(c["top_cut"], c["bottom_cut"]) for c in saved_cuts}
    requested = {(item.top_cut, item.bottom_cut): item.probability for item in request.cuts}

    # 每个已保存物理切口必须恰好出现一次：缺失或多余均 -> 422
    missing = sorted(saved_keys - set(requested))
    extra = sorted(set(requested) - saved_keys)
    if missing or extra:
        problems = []
        if missing:
            problems.append(
                "missing probabilities for cuts: "
                + ", ".join(f"({t}, {b})" for t, b in missing)
            )
        if extra:
            problems.append(
                "probabilities given for unknown cuts: "
                + ", ".join(f"({t}, {b})" for t, b in extra)
            )
        raise HTTPException(status_code=422, detail="; ".join(problems))

    # 概率按已保存切口的顺序（(top_cut, bottom_cut) 升序）排列，
    # 与请求条目的传入顺序无关，保证枚举与结果确定
    probabilities = [
        requested[(c["top_cut"], c["bottom_cut"])] for c in saved_cuts
    ]

    # 仅在内存中枚举成功/失败状态并重建片段，原记录与数据库均不改动
    distribution = partial_digest_distribution(
        record["sequence"], record["topology"], saved_cuts, probabilities
    )

    return {
        "id": record["id"],
        "topology": record["topology"],
        "sequence_length": len(record["sequence"]),
        "cut_count": len(saved_cuts),
        "outcome_count": len(distribution["outcomes"]),
        "total_probability": distribution["total_probability"],
        "outcomes": distribution["outcomes"],
    }


@app.post(
    "/digests/{submission_id}/methylation-analysis",
    response_model=MethylationAnalysisResponse,
)
def methylation_analysis(submission_id: int, request: MethylationAnalysisRequest) -> dict:
    # 先确认记录存在：不存在返回 404，且不泄露任何校验细节
    record = db.get_submission(submission_id)
    if record is None:
        raise HTTPException(status_code=404, detail=f"submission {submission_id} not found")

    length = len(record["sequence"])
    problems = []

    # 甲基化坐标依赖记录的序列长度，只能在取到记录后校验 -> 422
    out_of_bounds = sorted(
        {p for p in request.methylated_positions if p >= length}
    )
    if out_of_bounds:
        problems.append(
            f"methylated positions out of bounds for sequence length {length}: "
            + ", ".join(str(p) for p in out_of_bounds)
        )

    # 酶名与模式长度均来自已保存记录：未知酶名、偏移超出该酶模式范围 -> 422
    pattern_lengths: dict[str, set[int]] = {}
    for enzyme in record["enzymes"]:
        pattern_lengths.setdefault(enzyme["name"], set()).add(len(enzyme["pattern"]))

    unknown = sorted(
        {rule.enzyme for rule in request.enzyme_rules if rule.enzyme not in pattern_lengths}
    )
    if unknown:
        problems.append("unknown enzyme names: " + ", ".join(unknown))

    for rule in sorted(request.enzyme_rules, key=lambda r: r.enzyme):
        lengths = pattern_lengths.get(rule.enzyme)
        if not lengths:
            continue
        # 同名酶可能以不同模式长度保存：偏移须对其所有模式都合法
        max_offset = min(lengths) - 1
        bad = sorted({offset for offset in rule.sensitive_offsets if offset > max_offset})
        if bad:
            problems.append(
                f"sensitive offsets out of range for enzyme '{rule.enzyme}' "
                f"(max offset {max_offset}): " + ", ".join(str(o) for o in bad)
            )
    if problems:
        raise HTTPException(status_code=422, detail="; ".join(problems))

    rules = [
        (rule.enzyme, rule.sensitive_offsets) for rule in request.enzyme_rules
    ]

    # 仅在内存中评估阻断并复用片段重建，原记录与数据库均不改动
    analysis = analyze_methylation(
        record["sequence"],
        record["topology"],
        record["result"]["cuts"],
        rules,
        request.methylated_positions,
    )

    return {
        "id": record["id"],
        "topology": record["topology"],
        "sequence_length": length,
        "methylated_positions": sorted(set(request.methylated_positions)),
        **analysis,
    }


@app.post("/gel-comparisons", response_model=GelComparisonResponse)
def compare_gels(request: GelComparisonRequest) -> dict:
    # 逐条确认记录存在：任一不存在即 404，按请求顺序报告第一个缺失 id
    lanes = []
    for digest_id in request.digest_ids:
        record = db.get_submission(digest_id)
        if record is None:
            raise HTTPException(
                status_code=404, detail=f"submission {digest_id} not found"
            )
        lanes.append(
            {
                "digest_id": digest_id,
                "fragment_lengths": [
                    f["length"] for f in record["result"]["fragments"]
                ],
            }
        )

    # 仅在内存中汇总片段长度并聚类，数据库不做任何写入
    bands = compare_gel_lanes(lanes, request.resolution_ratio)
    return {
        "digest_ids": list(request.digest_ids),
        "resolution_ratio": request.resolution_ratio,
        "lane_count": len(lanes),
        "band_count": len(bands),
        "bands": bands,
    }


@app.post("/ligation-compatibility", response_model=LigationCompatibilityResponse)
def ligation_compatibility(request: LigationCompatibilityRequest) -> dict:
    # 逐条确认记录存在：任一不存在即 404，按请求顺序报告第一个缺失 id
    records = []
    for digest_id in request.digest_ids:
        record = db.get_submission(digest_id)
        if record is None:
            raise HTTPException(
                status_code=404, detail=f"submission {digest_id} not found"
            )
        records.append(
            {"digest_id": digest_id, "fragments": record["result"]["fragments"]}
        )

    # 仅在内存中枚举跨记录片段末端对，数据库不做任何写入
    pairs = find_compatible_pairs(records, request.include_blunt)
    return {
        "digest_ids": list(request.digest_ids),
        "include_blunt": request.include_blunt,
        "record_count": len(records),
        "pair_count": len(pairs),
        "pairs": pairs,
    }


@app.post("/ligation-assembly", response_model=LigationAssemblyResponse)
def ligation_assembly(request: LigationAssemblyRequest) -> dict:
    # 逐条确认记录存在：任一不存在即 404，按请求顺序报告第一个缺失 id
    records = []
    for digest_id in request.digest_ids:
        record = db.get_submission(digest_id)
        if record is None:
            raise HTTPException(
                status_code=404, detail=f"submission {digest_id} not found"
            )
        records.append(
            {"digest_id": digest_id, "fragments": record["result"]["fragments"]}
        )

    # 记忆化全枚举的状态空间随片段数与记录数膨胀：超过固定上限 -> 422
    total_fragments = sum(len(r["fragments"]) for r in records)
    if total_fragments > MAX_ASSEMBLY_FRAGMENTS:
        raise HTTPException(
            status_code=422,
            detail=(
                f"{total_fragments} candidate fragments across "
                f"{len(records)} records exceed the assembly planning limit of "
                f"{MAX_ASSEMBLY_FRAGMENTS}"
            ),
        )

    # 仅在内存中构建有向图并搜索最长拼接链，数据库不做任何写入
    plan = plan_assembly_chain(
        records, request.include_blunt, request.max_chain_length
    )
    used_positions = plan["used_positions"]
    return {
        "digest_ids": list(request.digest_ids),
        "include_blunt": request.include_blunt,
        "max_chain_length": request.max_chain_length,
        "record_count": len(records),
        "node_count": len(plan["nodes"]),
        "connection_count": len(plan["connections"]),
        "nodes": [public_node(node) for node in plan["nodes"]],
        "connections": plan["connections"],
        "used_digest_ids": [
            digest_id
            for position, digest_id in enumerate(request.digest_ids)
            if position in used_positions
        ],
        "unused_digest_ids": [
            digest_id
            for position, digest_id in enumerate(request.digest_ids)
            if position not in used_positions
        ],
        "stop_reason": plan["stop_reason"],
        "stop_reason_detail": plan["stop_reason_detail"],
    }
