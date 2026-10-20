#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""T02 文献筛选 | 对候选池执行纳排判定并生成筛选表

输入：literature/metadata/candidates_raw_20261001.csv 与 _R2.csv（由 lit_search_openalex.py 生成）
输出：literature/metadata/screening_20261001.csv    （全部候选 + 判定 + 理由 + 标签）
      literature/metadata/included_20261001.csv     （纳入条目，含角色/维度/精读优先级）
      literature/metadata/included_20261001.md      （人读版纳入清单）

判定由本轮负责人基于题名、年份、来源、类型与摘要片段（OpenAlex 提供）逐条作出。
未列入 INCLUDED / DUPLICATE 的条目统一使用通用排除理由。
"""

from __future__ import annotations

import csv
import os
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STAMP = "20261001"

# key -> (decision, role, dimension, priority, reason)
# decision: INCLUDE / EXCLUDE / DUPLICATE
# priority: 1=精读优先（T03 首批），2=次批，3=背景备查
D = {
    # ---------- A 组：疾病机制·温度依赖 ----------
    "A1#1": ("INCLUDE", "奠基·参数源", "机制", 1, "本文温度性状函数（Mordecai 2017）的直接来源，R0(T) 与 μ_v 的出处"),
    "A1#4": ("INCLUDE", "机制对照", "机制", 2, "季节温度变率影响气候适宜性，与本文'季节窗'设定直接对照"),
    "A1#5": ("INCLUDE", "参数源", "机制", 2, "白纹/埃及伊蚊温度约束，参数边界的外部依据"),
    "A1#6": ("DUPLICATE", "版本关系", "机制", 3, "Mordecai 2017 的 bioRxiv 预印本（10.1101/063735，2016）；与 A1#1 为同一工作，保留版本关系不重复计入"),
    "A1#8": ("INCLUDE", "奠基·背景", "机制", 3, "温度-登革传播的经典实证（波多黎各），奠基背景"),
    "A2#4": ("INCLUDE", "参数源", "机制", 1, "成蚊存活的温度建模（实验室与野外），直接关系 μ_v 与本文 μ_s 的生物学对照"),
    "A2#6": ("INCLUDE", "最接近工作候选", "干预", 1, "白纹伊蚊新定殖区输入基孔肯雅病例的数学评估，与本文'输入→本地传播'同构"),
    "A3#2": ("INCLUDE", "机制对照", "机制", 1, "温度经 R0 与生成间隔调制登革流行增长率，与本文 R0(T) 同题"),
    "A3#3": ("INCLUDE", "参数综述", "机制", 2, "温度-传播实验性状的系统综述，参数来源核查用"),
    "A3#4": ("INCLUDE", "机制对照", "机制", 2, "考虑温度依赖生成间隔的有效再生数估计"),
    # ---------- B 组：模型结构 ----------
    "B1#1": ("INCLUDE", "结构对照", "结构", 3, "模糊 SEIR-SEI 登革动力学，模型族对照"),
    "B1#2": ("INCLUDE", "结构对照", "结构", 2, "媒介-宿主 SEIR-SEI 登革模型分析，与本文仓室结构同构"),
    "B1#3": ("INCLUDE", "结构对照", "结构", 2, "SEIR-SEI + EnKF 估计与预测，模型+数据同化路线对照"),
    "B1#5": ("INCLUDE", "结构对照", "结构", 2, "登革传播动态建模与疾病持续性，与本文'无正地方病平衡点'结论对照"),
    "B1#6": ("INCLUDE", "数学对照", "数学", 3, "模糊 SEIR-SEI 敏感性/分岔分析"),
    "B2#1": ("INCLUDE", "机制对照", "机制", 2, "日温度波动对登革传播的影响（PNAS），与本文逐日温度强迫设定相关"),
    "B2#7": ("INCLUDE", "参数源", "机制", 2, "气象变量与蚊媒监测预测种群趋势，支撑'丰度非恒定'的观测依据"),
    # ---------- C 组：数学性质 ----------
    "C1#2": ("INCLUDE", "数学奠基", "数学", 1, "周期环境下仓室模型的阈值动力学（Wang & Zhao 型），时变 R0 的理论基础"),
    "C1#3": ("INCLUDE", "数学奠基", "数学", 1, "传染病模型再生数的统一定义，下一代矩阵方法学依据"),
    "C1#5": ("INCLUDE", "数学核心", "数学", 1, "描述登革种群动力学的**非自治**系统的阈值条件，与本文时变系统最直接相关"),
    "C1#6": ("INCLUDE", "数学对照", "数学", 3, "含直接传播与时滞的媒介传播模型"),
    "C2#3": ("INCLUDE", "数学对照", "数学", 3, "媒介传播病 R0 制图方法与案例"),
    "C2#5": ("INCLUDE", "数学对照", "数学", 3, "最终规模与 R0 的关系，与本文建立阈值（100 例）口径相关"),
    "C2#6": ("INCLUDE", "数学奠基", "数学", 1, "R0 定义与计算的奠基文献（Diekmann 等）"),
    "C3#3": ("INCLUDE", "方法指南", "数学", 2, "传染病爆发风险评估的数学方法实用指南，与本文建立概率方法学定位相关"),
    # ---------- D 组：干预与观测 ----------
    "D2#1": ("INCLUDE", "结构综述", "结构", 1, "登革传播动力学模型的结构方法系统综述（PLoS ONE），模型族竞赛的对照基准"),
    "D2#2": ("INCLUDE", "最接近工作", "干预", 2, "马德拉 2012 登革暴发：流行病学决定因素与未来流行潜能，输入引发本地暴发的实例"),
    "D2#5": ("INCLUDE", "背景", "干预", 3, "航空旅行与登革病毒洲内扩散，输入压力背景"),
    "D2#6": ("INCLUDE", "对照", "干预", 3, "气候变化对美国东南部登革传播影响的建模"),
    "D3#4": ("INCLUDE", "背景", "干预", 3, "人类移动对登革流行出现的影响"),
    "D3#5": ("INCLUDE", "背景", "机制", 3, "巴西登革时空动态：季节行波与区域同步"),
    # ---------- E 组：区域背景 ----------
    "E1#1": ("INCLUDE", "区域背景", "区域", 2, "广州极端天气事件与登革暴发（时间序列），案例区域背景"),
    "E1#3": ("INCLUDE", "最接近工作", "区域", 1, "用输入病例与气候因子预测广州 2014 异常登革暴发，与本文问题高度接近"),
    "E1#4": ("INCLUDE", "区域背景", "区域", 3, "广州登革与天气的时间序列分析（非机理模型，作为方法边界对照）"),
    "E1#5": ("INCLUDE", "背景", "机制", 3, "气候变率驱动登革动态"),
    "E1#6": ("INCLUDE", "最接近工作", "区域", 1, "气候与输入病例时点作为广州 2014 暴发决定因素，与本文'引入时点'问题最直接"),
    "E1#7": ("INCLUDE", "区域背景", "区域", 3, "中国大陆登革热概况，区域背景"),
    "E2#3": ("INCLUDE", "最接近工作", "区域", 1, "经输入病例、蚊媒与气候影响预测广州本地登革传播，同案例同问题"),
    # ---------- 第二轮补充 ----------
    "F2#6": ("INCLUDE", "最接近工作", "干预", 2, "基于种群动力学评估北意登革/基孔肯雅暴发潜在风险，输入引发本地传播风险"),
    "F2#8": ("INCLUDE", "对照", "干预", 3, "欧洲大城市登革/基孔/寨卡暴发风险估计"),
    "F3#2": ("INCLUDE", "数学核心", "数学", 1, "确定性与随机流行病模型的灭绝阈值，本文建立概率的理论对照"),
    "F3#5": ("INCLUDE", "数学对照", "数学", 2, "随机多斑块模型中疾病爆发的概率"),
    "F3#6": ("INCLUDE", "方法", "数学", 3, "随机流行病模型与推断（专著），tau-leap 与推断的方法学背景"),
    "G1#6": ("INCLUDE", "观测核心", "观测", 1, "用扩展因子估计东南亚登革负担，与本文报告率 ρ 的外部对照（漏报乘数）"),
    "G2#2": ("INCLUDE", "观测核心", "观测", 1, "隐性感染（沉默多数）主导登革传播，支撑'报告≠真实'的观测层设定"),
    "G2#8": ("INCLUDE", "观测核心", "观测", 1, "监测数据中报告延迟的建模校正，与本文延迟 d 直接对应"),
    "G3#8": ("INCLUDE", "观测对照", "观测", 2, "精化登革全球负担估计的方法，漏报修正路径"),
    "H1#2": ("INCLUDE", "机制对照", "机制", 3, "温度对埃及伊蚊种群的影响评估"),
    "H1#7": ("INCLUDE", "机制对照", "机制", 3, "基于气候的登革流行理解与预测模型"),
    "H2#2": ("INCLUDE", "机制核心", "机制", 1, "埃及伊蚊成蚊存活率在越南呈**季节变化**，与本文 μ_s（有效寿命）直接对照"),
    "H2#3": ("INCLUDE", "机制对照", "机制", 3, "异质性、混合与蚊媒传播的空间尺度"),
    "H2#5": ("INCLUDE", "机制对照", "机制", 3, "气候预测蚊媒病动态的地理与时间变异"),
    "H2#8": ("INCLUDE", "结构综述", "结构", 1, "登革传播与媒介控制数学模型的系统综述（2023），用于核对本文结构是否已被覆盖"),
    "H3#2": ("INCLUDE", "机制背景", "机制", 3, "媒介表型可塑性在传播中的作用"),
    "H3#4": ("INCLUDE", "数据资源", "数据", 3, "公开登革病例计数全球数据集（2024），可作为多年度验证的数据候选"),
    "H3#6": ("INCLUDE", "近期对照", "干预", 2, "2025：基孔肯雅早期传播动力学与流行病学特征（输入引发本地传播的近期同类工作）"),
    "H3#8": ("INCLUDE", "近期背景", "机制", 3, "2024：变化气候下温度驱动的登革传播（模式、趋势与未来）"),
}

GENERIC_EXCLUDE = "未通过纳入标准 I1–I5：题名/摘要与四大维度（机制、结构、数学、干预观测）及本文案例均不相关"


def main() -> int:
    src = []
    for suffix in ("", "_R2", "_R3"):
        p = os.path.join(ROOT, "literature", "metadata", "candidates_raw_%s%s.csv" % (STAMP, suffix))
        with open(p, encoding="utf-8-sig") as f:
            for r in csv.DictReader(f):
                r["round"] = {"":"R1","_R2":"R2","_R3":"R3"}[suffix]
                src.append(r)

    out = []
    for r in src:
        key = r.get("key") or (r["query_id"] + "#" + r.get("rank", ""))
        if r["round"] == "R3" and r.get("decision"):
            decision, role, dim, prio, reason = (
                r["decision"], r.get("role", ""), r.get("dimension", ""),
                int(r.get("read_priority") or 3), r.get("reason", ""))
        else:
            decision, role, dim, prio, reason = D.get(
                key, ("EXCLUDE", "", "", "", GENERIC_EXCLUDE))
        out.append({
            "key": key, "round": r["round"], "query_id": r["query_id"], "group_name": r["group_name"],
            "decision": decision, "role": role, "dimension": dim, "read_priority": prio,
            "doi": r["doi"], "title": r["title"], "year": r["year"], "venue": r["venue"],
            "type": r["type"], "is_preprint": r["is_preprint"], "cited_by": r["cited_by"],
            "authors": r["authors"], "url": r["url"], "reason": reason,
        })

    sdir = os.path.join(ROOT, "literature", "metadata")
    spath = os.path.join(sdir, "screening_%s.csv" % STAMP)
    with open(spath, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0].keys()))
        w.writeheader()
        w.writerows(out)

    inc = [o for o in out if o["decision"] == "INCLUDE"]
    inc.sort(key=lambda o: (o["read_priority"], -int(o["cited_by"] or 0)))
    ipath = os.path.join(sdir, "included_%s.csv" % STAMP)
    with open(ipath, "w", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(inc[0].keys()))
        w.writeheader()
        w.writerows(inc)

    mpath = os.path.join(sdir, "included_%s.md" % STAMP)
    with open(mpath, "w", encoding="utf-8") as f:
        f.write("# 纳入文献（T02 筛选后）— %s\n\n" % STAMP)
        f.write("共 %d 条（按精读优先级排序：1=首批，2=次批，3=背景备查）。\n" % len(inc))
        f.write("DOI 均来自 OpenAlex 返回；**未做任何人工补全**。\n\n")
        for p in ("1", "2", "3"):
            sub = [i for i in inc if str(i["read_priority"]) == p]
            if not sub:
                continue
            f.write("\n## 精读优先级 %s（%d 条）\n\n" % (p, len(sub)))
            for i in sub:
                f.write("- **[%s] %s** (%s, %s)\n" % (i["key"], i["title"][:150], i["year"], i["venue"][:60]))
                f.write("  - DOI: `%s`｜角色：%s｜维度：%s\n" % (i["doi"] or "NA", i["role"], i["dimension"]))
                f.write("  - 作者：%s\n  - 理由：%s\n" % (i["authors"], i["reason"]))
    print("candidates:", len(out), "| INCLUDE:", len(inc),
          "| DUPLICATE:", sum(1 for o in out if o["decision"] == "DUPLICATE"),
          "| EXCLUDE:", sum(1 for o in out if o["decision"] == "EXCLUDE"))
    print("WROTE", spath)
    print("WROTE", ipath)
    print("WROTE", mpath)
    return 0


if __name__ == "__main__":
    sys.exit(main())
