#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""T00 初始化 | 内容哈希与版本指纹采集

用途
----
本项目不自动跟随技能仓库更新，因此对每个技能目录、专家配置文件和数据文件
计算内容指纹，写入 environment/hashes_T00.json，供 protocols/SKILL_MAP.md 引用。

哈希口径
--------
* file_sha256 : 单文件 SHA-256。
* dir_sha256  : 目录指纹。按相对路径排序，逐文件计算 SHA-256，
                再对 "relpath<TAB>sha256\n" 的拼接串整体取 SHA-256。
                目录指纹随任一文件内容或文件名变化而变化。

用法
----
    python scripts/compute_hashes.py            # 写入 environment/hashes_T00.json
    python scripts/compute_hashes.py --stdout   # 只打印，不写文件

脚本只做只读计算，不修改任何被哈希的对象。
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from datetime import datetime, timezone, timedelta

# 项目根目录 = 本脚本所在目录的上一级
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "environment", "hashes_T00.json")

SKILLS = os.path.expanduser("~/.workbuddy/skills")
SAS = r"D:\Environment\Skills\scientific-agent-skills\skills"
SA = r"D:\Environment\Skills\scientific-agents\scientific-agents"
ARIS = r"D:\Environment\Skills\codex\aris-deepseek\skills"
EPI = r"D:\Environment\Skills\codex\research-profiles\epiagent"
NAT = r"D:\Environment\Skills\nature-skills\skills"
EXPERTS = r"D:\Environment\Skills\codex\research-profiles"

# (分类, 显示名, 绝对路径, 是否必需)
TARGETS = [
    # --- K-Dense scientific-agent-skills（SKILL.md 技能）---
    ("sas", "matlab", os.path.join(SKILLS, "matlab"), True),
    ("sas", "sympy", os.path.join(SKILLS, "sympy"), True),
    ("sas", "pymc", os.path.join(SKILLS, "pymc"), False),
    ("sas", "networkx", os.path.join(SKILLS, "networkx"), False),
    # --- scientific-agents（专家 AGENTS.md，不是 skill）---
    ("expert", "dynamical-systems-theorist", os.path.join(SA, "dynamical-systems-theorist"), True),
    ("expert", "mathematical-modeler", os.path.join(SA, "mathematical-modeler"), True),
    ("expert", "epidemiologist", os.path.join(SA, "epidemiologist"), True),
    ("expert", "numerical-analyst", os.path.join(SA, "numerical-analyst"), False),
    ("expert", "bayesian-statistician", os.path.join(SA, "bayesian-statistician"), False),
    # --- ARIS ---
    ("aris", "research-refine", os.path.join(SKILLS, "research-refine"), True),
    ("aris", "research-lit", os.path.join(SKILLS, "research-lit"), True),
    ("aris", "novelty-check", os.path.join(SKILLS, "novelty-check"), True),
    ("aris", "formula-derivation", os.path.join(SKILLS, "formula-derivation"), True),
    ("aris", "proof-writer", os.path.join(SKILLS, "proof-writer"), True),
    ("aris", "proof-checker", os.path.join(SKILLS, "proof-checker"), True),
    ("aris", "experiment-plan", os.path.join(SKILLS, "experiment-plan"), True),
    ("aris", "result-to-claim", os.path.join(SKILLS, "result-to-claim"), True),
    ("aris", "idea-creator", os.path.join(SKILLS, "idea-creator"), True),
    ("aris", "analyze-results", os.path.join(SKILLS, "analyze-results"), True),
    ("aris", "ablation-planner", os.path.join(SKILLS, "ablation-planner"), True),
    ("aris", "paper-plan", os.path.join(SKILLS, "paper-plan"), True),
    ("aris", "paper-compile", os.path.join(SKILLS, "paper-compile"), True),
    ("aris", "research-pipeline", os.path.join(SKILLS, "research-pipeline"), False),
    ("aris", "paper-claim-audit", os.path.join(SKILLS, "paper-claim-audit"), False),
    ("aris", "citation-audit", os.path.join(SKILLS, "citation-audit"), False),
    ("aris", "experiment-audit", os.path.join(SKILLS, "experiment-audit"), False),
    ("aris", "proof-orchestrator", os.path.join(SKILLS, "proof-orchestrator"), False),
    ("aris", "paper-write", os.path.join(SKILLS, "paper-write"), False),
    ("aris", "paper-figure", os.path.join(SKILLS, "paper-figure"), False),
    ("aris", "rebuttal", os.path.join(SKILLS, "rebuttal"), False),
    ("aris", "resubmit-pipeline", os.path.join(SKILLS, "resubmit-pipeline"), False),
    ("aris", "shared-references", os.path.join(ARIS, "shared-references"), True),
    # --- Epiagent（R/Epiverse）---
    ("epi", "epiverse-overview", os.path.join(SKILLS, "epiverse-overview"), True),
    ("epi", "data-intake", os.path.join(SKILLS, "data-intake"), True),
    ("epi", "parameters", os.path.join(SKILLS, "parameters"), True),
    ("epi", "simulation", os.path.join(SKILLS, "simulation"), True),
    ("epi", "analysis", os.path.join(SKILLS, "analysis"), True),
    ("epi", "visualisation", os.path.join(SKILLS, "visualisation"), False),
    ("epi", "reporting", os.path.join(EPI, "reporting"), False),
    # --- nature-skills ---
    ("nature", "nature-reader", os.path.join(SKILLS, "nature-reader"), True),
    ("nature", "nature-figure", os.path.join(SKILLS, "nature-figure"), True),
    ("nature", "nature-writing", os.path.join(SKILLS, "nature-writing"), True),
    ("nature", "nature-polishing", os.path.join(SKILLS, "nature-polishing"), True),
    ("nature", "nature-reviewer", os.path.join(SKILLS, "nature-reviewer"), True),
    ("nature", "nature-response", os.path.join(SKILLS, "nature-response"), True),
    ("nature", "nature-shared", os.path.join(SKILLS, "nature-shared"), True),
    ("nature", "nature-downloader", os.path.join(SKILLS, "nature-downloader"), True),
    ("nature", "nature-academic-search", os.path.join(SKILLS, "nature-academic-search"), False),
    ("nature", "nature-statistics", os.path.join(SKILLS, "nature-statistics"), False),
    ("nature", "nature-citation", os.path.join(SKILLS, "nature-citation"), False),
]

DATA_FILES = [
    "guangzhou_2025_weekly_model_data.csv",
    "guangzhou_2025_weekly_cases_official.csv",
    "guangzhou_2025_daily_cases_reconstructed.csv",
    "guangzhou_2025_daily_weather_openmeteo.csv",
    "guangzhou_2025-2026_daily_weather_extension_openmeteo.csv",
    "guangzhou_2025_monthly_cases_observed.csv",
    "guangzhou_population.csv",
    "dataset_validation.json",
]


def sha256_file(path: str, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def dir_fingerprint(path: str):
    """返回 (dir_sha256, file_count, total_bytes, entry_list)。"""
    entries = []
    for base, dirs, files in os.walk(path):
        dirs.sort()
        for name in sorted(files):
            full = os.path.join(base, name)
            rel = os.path.relpath(full, path).replace("\\", "/")
            try:
                entries.append((rel, sha256_file(full), os.path.getsize(full)))
            except OSError:
                entries.append((rel, "UNREADABLE", 0))
    entries.sort(key=lambda x: x[0])
    agg = hashlib.sha256()
    for rel, digest, _ in entries:
        agg.update(("%s\t%s\n" % (rel, digest)).encode("utf-8"))
    return agg.hexdigest(), len(entries), sum(e[2] for e in entries), entries


def frontmatter_name(skill_md: str):
    try:
        with open(skill_md, "r", encoding="utf-8", errors="replace") as f:
            head = f.read(4096)
    except OSError:
        return None
    if not head.startswith("---"):
        return None
    end = head.find("\n---", 3)
    if end < 0:
        return None
    for line in head[3:end].splitlines():
        if line.startswith("name:"):
            return line.split(":", 1)[1].strip().strip('"').strip("'")
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stdout", action="store_true", help="只打印 JSON，不写文件")
    args = ap.parse_args()

    result = {
        "generated_at": datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds"),
        "method": {
            "file_sha256": "SHA-256 of raw bytes",
            "dir_sha256": "SHA-256 over sorted '<relpath>\\t<sha256>\\n' concatenation",
        },
        "targets": [],
        "data_files": [],
        "repos": {},
    }

    for cat, name, path, required in TARGETS:
        rec = {
            "category": cat,
            "name": name,
            "path": path.replace("\\", "/"),
            "required": required,
            "status": "AVAILABLE" if os.path.exists(path) else "MISSING",
        }
        if os.path.isdir(path):
            d, n, b, entries = dir_fingerprint(path)
            rec.update({"type": "dir", "dir_sha256": d, "files": n, "bytes": b})
            md = os.path.join(path, "SKILL.md")
            if os.path.isfile(md):
                rec["frontmatter_name"] = frontmatter_name(md)
                rec["skill_md_sha256"] = sha256_file(md)
            ag = os.path.join(path, "AGENTS.md")
            if os.path.isfile(ag):
                rec["has_AGENTS_md"] = True
                rec["AGENTS_md_sha256"] = sha256_file(ag)
                rec["frontmatter_name"] = None  # 专家配置没有 frontmatter
            rec["entry_list"] = [{"p": e[0], "sha256": e[1]} for e in entries[:200]]
        elif os.path.isfile(path):
            rec.update({"type": "file", "file_sha256": sha256_file(path)})
        result["targets"].append(rec)

    for rel in DATA_FILES:
        p = os.path.join(ROOT, "data", rel)
        if os.path.isfile(p):
            result["data_files"].append(
                {
                    "path": os.path.relpath(p, ROOT).replace("\\", "/"),
                    "sha256": sha256_file(p),
                    "bytes": os.path.getsize(p),
                    "status": "AVAILABLE",
                }
            )
        else:
            result["data_files"].append(
                {"path": "data/" + rel, "sha256": None, "bytes": None, "status": "MISSING"}
            )

    text = json.dumps(result, ensure_ascii=False, indent=2)
    if args.stdout:
        print(text)
        return 0
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as f:
        f.write(text)
    print("WROTE", OUT)
    miss = [t["name"] for t in result["targets"] if t["status"] == "MISSING" and t["required"]]
    if miss:
        print("MISSING(required):", ", ".join(miss))
    return 0


if __name__ == "__main__":
    sys.exit(main())
