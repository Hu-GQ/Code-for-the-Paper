#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""T02 文献检索 | OpenAlex 检索执行器（只读、可重跑）

用途
----
按预先定义的关键词组与检索式访问 OpenAlex REST API，记录：
  * 每条检索式命中的**实际返回数量**（meta.count）
  * 时间范围过滤（无过滤 / 2022-01-01 起）
  * 命中记录的元数据：DOI、题名、作者、年份、来源、URL、被引数、类型、预印本标识
输出：
  literature/search_logs/openalex_raw_<date>.json   原始返回（可审计）
  literature/metadata/candidates_raw_<date>.csv      候选池（未筛选）
  literature/metadata/candidates_raw_<date>.md       候选池（人读，含摘要片段）

约定
----
* 只使用当前实际可访问的数据源（OpenAlex 公开 API，无需密钥）。
* 不编造 DOI / 页码 / 来源；缺失字段一律写 NA 或空。
* 预印本与正式版不做自动合并，只标记类型，版本关系由人工在筛选时判定。
* 脚本只读，不修改任何既有产物。

用法
----
    python scripts/lit_search_openalex.py
"""

from __future__ import annotations

import csv
import json
import os
import sys
import time
import urllib.parse
import urllib.request
from datetime import datetime, timezone, timedelta

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAMP = datetime.now(timezone(timedelta(hours=8))).strftime("%Y%m%d")
API = "https://api.openalex.org/works"
MAILTO = "paper01@example.org"  # OpenAlex polite pool
PER_PAGE = 8
SLEEP = 0.4

# (组ID, 组名, 检索式ID, 检索式, 是否奠基检索)
QUERIES = [
    ("A", "疾病机制·温度依赖传播", "A1", 'dengue temperature dependent transmission mechanistic model', False),
    ("A", "疾病机制·温度依赖传播", "A2", 'Aedes albopictus thermal biology dengue transmission', False),
    ("A", "疾病机制·温度依赖传播", "A3", 'Mordecai dengue R0 temperature mechanistic', True),
    ("B", "模型结构·人蚊仓室与季节丰度", "B1", 'dengue SEIR SEI vector-host compartmental model', False),
    ("B", "模型结构·人蚊仓室与季节丰度", "B2", 'seasonal mosquito abundance dengue transmission model', False),
    ("B", "模型结构·人蚊仓室与季节丰度", "B3", 'time-varying mosquito-to-human ratio dengue model', False),
    ("C", "数学性质·R0/平衡/随机灭绝", "C1", 'time-varying basic reproduction number non-autonomous epidemic model', False),
    ("C", "数学性质·R0/平衡/随机灭绝", "C2", 'next generation matrix vector-borne disease basic reproduction number', True),
    ("C", "数学性质·R0/平衡/随机灭绝", "C3", 'establishment probability imported case outbreak branching process extinction', False),
    ("D", "干预与观测·漏报/控制时点", "D1", 'dengue under-reporting reporting rate estimation', False),
    ("D", "干预与观测·漏报/控制时点", "D2", 'dengue importation outbreak stochastic extinction probability', False),
    ("D", "干预与观测·漏报/控制时点", "D3", 'dengue vector control timing seasonal intervention modelling', False),
    ("E", "区域背景·广州/中国登革热", "E1", 'Guangzhou dengue outbreak modelling China', False),
    ("E", "区域背景·广州/中国登革热", "E2", 'China dengue imported cases local transmission seasonality', False),
]

RECENT_FROM = "2022-01-01"


def get(params: dict):
    url = API + "?" + urllib.parse.urlencode(params)
    req = urllib.request.Request(url, headers={"User-Agent": "paper01-lit-search/1.0 (mailto:%s)" % MAILTO})
    with urllib.request.urlopen(req, timeout=60) as r:
        return json.load(r)


def abstract_from_index(inv):
    if not inv:
        return ""
    pos = {}
    for word, idxs in inv.items():
        for i in idxs:
            pos[i] = word
    return " ".join(pos[i] for i in sorted(pos))


def host_venue(w):
    pl = w.get("primary_location") or {}
    src = pl.get("source") or {}
    return src.get("display_name") or ""


def doi_url(w):
    doi = w.get("doi")
    if doi:
        return doi
    return (pl_landing(w) or "")


def pl_landing(w):
    pl = w.get("primary_location") or {}
    return pl.get("landing_page_url") or ""


def main() -> int:
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--queries-file", help="JSON 文件：[[组ID,组名,检索式ID,检索式,是否奠基], ...]；用于补充检索轮次")
    ap.add_argument("--suffix", default="", help="输出文件名后缀，例如 _R2")
    args = ap.parse_args()

    queries = QUERIES
    if args.queries_file:
        with open(args.queries_file, "r", encoding="utf-8") as f:
            queries = [tuple(x) for x in json.load(f)]

    out_dir_json = os.path.join(ROOT, "literature", "search_logs")
    out_dir_md = os.path.join(ROOT, "literature", "metadata")
    os.makedirs(out_dir_json, exist_ok=True)
    os.makedirs(out_dir_md, exist_ok=True)
    suffix = args.suffix

    raw = {"generated_at": datetime.now(timezone(timedelta(hours=8))).isoformat(timespec="seconds"),
           "database": "OpenAlex REST API (https://api.openalex.org/works)",
           "access": "public, no key; polite pool mailto=%s" % MAILTO,
           "per_query_fetch": PER_PAGE,
           "recent_from": RECENT_FROM,
           "runs": []}

    rows = []
    for gid, gname, qid, qtext, foundational in queries:
        entry = {"group": gid, "group_name": gname, "query_id": qid, "query": qtext,
                 "foundational": foundational, "searches": {}}
        for tag, date_filter in (("all_years", None), ("recent", RECENT_FROM)):
            params = {"search": qtext, "per-page": PER_PAGE, "mailto": MAILTO}
            if date_filter:
                params["filter"] = "from_publication_date:" + date_filter
            try:
                d = get(params)
            except Exception as e:  # 网络失败必须记录，不得静默
                entry["searches"][tag] = {"status": "FAILED", "error": "%s: %s" % (type(e).__name__, e)}
                print("FAILED", qid, tag, e)
                time.sleep(SLEEP)
                continue
            meta = d.get("meta", {})
            entry["searches"][tag] = {
                "status": "OK",
                "filter": ("from_publication_date:" + date_filter) if date_filter else "none",
                "count": meta.get("count"),
                "returned": len(d.get("results", [])),
                "request_url": API + "?" + urllib.parse.urlencode(params),
            }
            if tag == "all_years":
                for rank, w in enumerate(d.get("results", []), start=1):
                    au = [a.get("author", {}).get("display_name", "") for a in (w.get("authorships") or [])]
                    rows.append({
                        "query_id": qid,
                        "group": gid,
                        "group_name": gname,
                        "rank": rank,
                        "doi": w.get("doi") or "",
                        "title": (w.get("title") or "").replace("\n", " "),
                        "year": w.get("publication_year") or "",
                        "venue": host_venue(w),
                        "type": w.get("type") or "",
                        "is_preprint": "yes" if (w.get("type") == "preprint" or
                                                 (w.get("primary_location") or {}).get("is_preprint")) else "no",
                        "cited_by": w.get("cited_by_count") or 0,
                        "url": w.get("doi") or pl_landing(w) or "",
                        "authors": "; ".join([a for a in au if a][:6]) + ("; et al." if len(au) > 6 else ""),
                        "abstract": abstract_from_index(w.get("abstract_inverted_index"))[:600],
                    })
            print("%s %-4s %-10s count=%s" % (qid, tag, "", meta.get("count")))
            time.sleep(SLEEP)
        raw["runs"].append(entry)

    jpath = os.path.join(out_dir_json, "openalex_raw_%s%s.json" % (STAMP, suffix))
    with open(jpath, "w", encoding="utf-8") as f:
        json.dump(raw, f, ensure_ascii=False, indent=2)

    cpath = os.path.join(out_dir_md, "candidates_raw_%s%s.csv" % (STAMP, suffix))
    with open(cpath, "w", encoding="utf-8-sig", newline="") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(rows)

    mpath = os.path.join(out_dir_md, "candidates_raw_%s%s.md" % (STAMP, suffix))
    with open(mpath, "w", encoding="utf-8") as f:
        f.write("# 候选池（未筛选）— OpenAlex，%s%s\n\n" % (STAMP, suffix))
        f.write("共 %d 条（含跨检索式重复）。字段：检索式 / DOI / 题名 / 年份 / 来源 / 类型 / 被引。\n\n" % len(rows))
        for r in rows:
            f.write("### %s | %s\n" % (r["query_id"], r["title"][:160]))
            f.write("- DOI: `%s`\n- 年份/来源: %s / %s\n- 类型: %s（预印本=%s）｜被引 %s\n- 作者: %s\n- URL: %s\n\n"
                    % (r["doi"] or "NA", r["year"], r["venue"], r["type"], r["is_preprint"],
                       r["cited_by"], r["authors"], r["url"]))
    print("WROTE", jpath)
    print("WROTE", cpath)
    print("WROTE", mpath)
    print("rows:", len(rows))
    return 0


if __name__ == "__main__":
    sys.exit(main())
