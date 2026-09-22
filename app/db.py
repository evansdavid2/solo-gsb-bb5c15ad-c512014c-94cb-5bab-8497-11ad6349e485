"""SQLite 持久化层。

数据库文件路径由环境变量 DATABASE_PATH 指定（默认 data/digest.db），
应用启动时自动创建目录与表结构。所有写入在单个事务中完成，
任一失败即整体回滚，不会产生半截数据。
"""

from __future__ import annotations

import json
import os
import sqlite3
from datetime import datetime, timezone

SCHEMA = """
CREATE TABLE IF NOT EXISTS submissions (
    id          INTEGER PRIMARY KEY AUTOINCREMENT,
    sequence    TEXT NOT NULL,
    topology    TEXT NOT NULL CHECK (topology IN ('linear', 'circular')),
    created_at  TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS enzymes (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    submission_id INTEGER NOT NULL REFERENCES submissions(id) ON DELETE CASCADE,
    name          TEXT NOT NULL,
    pattern       TEXT NOT NULL,
    cut_top       INTEGER NOT NULL,
    cut_bottom    INTEGER NOT NULL
);

CREATE TABLE IF NOT EXISTS results (
    submission_id INTEGER PRIMARY KEY REFERENCES submissions(id) ON DELETE CASCADE,
    result_json   TEXT NOT NULL
);

CREATE INDEX IF NOT EXISTS idx_enzymes_submission ON enzymes(submission_id);
"""


class Database:
    def __init__(self, path: str):
        self.path = path
        parent = os.path.dirname(os.path.abspath(path))
        os.makedirs(parent, exist_ok=True)
        self._init_schema()

    def _connect(self) -> sqlite3.Connection:
        conn = sqlite3.connect(self.path)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        return conn

    def _init_schema(self) -> None:
        with self._connect() as conn:
            conn.executescript(SCHEMA)

    def save_submission(self, sequence: str, topology: str, enzymes, result: dict):
        """在一次事务中保存提交、酶列表与计算结果，返回 (id, created_at)。"""
        created_at = datetime.now(timezone.utc).isoformat()
        with self._connect() as conn:  # 异常时自动回滚
            cur = conn.execute(
                "INSERT INTO submissions (sequence, topology, created_at) VALUES (?, ?, ?)",
                (sequence, topology, created_at),
            )
            submission_id = cur.lastrowid
            conn.executemany(
                "INSERT INTO enzymes (submission_id, name, pattern, cut_top, cut_bottom)"
                " VALUES (?, ?, ?, ?, ?)",
                [
                    (submission_id, e.name, e.pattern, e.cut_top, e.cut_bottom)
                    for e in enzymes
                ],
            )
            conn.execute(
                "INSERT INTO results (submission_id, result_json) VALUES (?, ?)",
                (submission_id, json.dumps(result, ensure_ascii=False)),
            )
        return submission_id, created_at

    def get_submission(self, submission_id: int):
        """按 id 查询完整记录，不存在时返回 None。"""
        with self._connect() as conn:
            sub = conn.execute(
                "SELECT id, sequence, topology, created_at FROM submissions WHERE id = ?",
                (submission_id,),
            ).fetchone()
            if sub is None:
                return None
            enzymes = conn.execute(
                "SELECT name, pattern, cut_top, cut_bottom FROM enzymes"
                " WHERE submission_id = ? ORDER BY id",
                (submission_id,),
            ).fetchall()
            result = conn.execute(
                "SELECT result_json FROM results WHERE submission_id = ?",
                (submission_id,),
            ).fetchone()
        return {
            "id": sub["id"],
            "sequence": sub["sequence"],
            "topology": sub["topology"],
            "created_at": sub["created_at"],
            "enzymes": [dict(e) for e in enzymes],
            "result": json.loads(result["result_json"]),
        }

    def list_submissions(self) -> list[dict]:
        """列出全部提交的摘要信息，按 id 升序（结果顺序确定）。"""
        with self._connect() as conn:
            subs = conn.execute(
                "SELECT s.id, s.sequence, s.topology, s.created_at,"
                "       (SELECT COUNT(*) FROM enzymes e WHERE e.submission_id = s.id) AS enzyme_count,"
                "       r.result_json AS result_json"
                " FROM submissions s JOIN results r ON r.submission_id = s.id"
                " ORDER BY s.id"
            ).fetchall()
        summaries = []
        for row in subs:
            result = json.loads(row["result_json"])
            summaries.append(
                {
                    "id": row["id"],
                    "topology": row["topology"],
                    "sequence_length": len(row["sequence"]),
                    "enzyme_count": row["enzyme_count"],
                    "cut_count": len(result["cuts"]),
                    "fragment_count": len(result["fragments"]),
                    "created_at": row["created_at"],
                }
            )
        return summaries
