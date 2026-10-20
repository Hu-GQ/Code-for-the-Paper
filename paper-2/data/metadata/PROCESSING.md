# 数据处理说明 PROCESSING

- 版本：v1.0｜建立：2026-10-05（T06）
- 源项目：`E:\Agent Work Space\Deepseek\text1` —— **全程只读**，本轮改动数 **0**（已核验）
- 授权：研究者 2026-10-05 明示"可以复制到本工作目录，但不要修改源目录下的数据"（决策 `D2026-10-05-01` §2.4）

---

## 1. 数据流

```text
源项目（只读）
  Data Analysis/raw_sources/2025/          一级原始材料
    guangdong_cdc_dengue_weekly_reports/week_06..52.html   47 份周报原文
    open_meteo_guangzhou_2025_era5.json                     ERA5 原始响应
    guangzhou_2025_population_source.html                   人口原始网页
    source_manifest.json                                    抓取清单（URL + sha256）
  Data Analysis/data/2025/                 上游清洗产物
    *.csv / dataset_validation.json
        │
        │  scripts/preprocess/t06_import_raw.py      （只读复制 + SHA-256）
        ▼
本项目  data/raw/                          （L1 原始 50 + L2 上游清洗产物 8，均保持不变）
        │
        │  scripts/preprocess/t06_build_analysis_input.py
        ▼
本项目  data/processed/weekly_analysis_input.csv       47 行 × 11 列，正式分析输入
        data/metadata/processed_checks.json            一致性校验结果
```

**上游清洗脚本**（未修改副本）：`scripts/preprocess/upstream_build_dataset_2025.py`（741 行，
逐字复制自 `Data Analysis/build_dataset_2025.py`）。它同时充当 raw→L2 的**可复现重建链**。

---

## 2. 各步骤说明

### 2.1 导入（`t06_import_raw.py`）

- 只**读取**源目录，绝不写入。
- 逐文件计算 SHA-256；若目标已存在且哈希一致则跳过，**若不一致则报错退出（不静默覆盖）**。
- 结果写入 `data/metadata/raw_file_hashes.csv`（58 行）。
- **交叉核验**：47 份周报 HTML 的哈希与上游 `source_manifest.json` 登记值 **47/47 一致**。

### 2.2 生成分析输入（`t06_build_analysis_input.py`）

- 以官方周表（47 周）为**主干**，按 `week_start` 左连接气象（来自 53 周建模表）。
- 同时输出两套周编号（`week_index` / `iso_week`），并校验 `iso_week − week_index ≡ 5`。
- 只保留有官方病例覆盖的 47 周；**ISO 第 1–5 周与第 53 周因病例不可得而排除，不作为 0 参与拟合**。
- 校验不通过即非零退出。全部通过，见 `processed_checks.json`。

---

## 3. 必须诚实声明的六条数据性质

| # | 性质 | 含义与处置 |
|---|---|---|
| 1 | 🔴 **日尺度病例是重构的** | 日表由周计数在周内**均匀分配**（`allocate_uniform`，余数按中→外顺序补）得到，**不是真实日数据**。⇒ 可用于月度聚合与绘图，**不得**用于日尺度统计推断（如日序列相关性、日延迟分布）。本项目拟合只用**周**数据。 |
| 2 | **病例按报告日期归周** | 官方周报无发病日期（发病日期缺失率高）。⇒ 报告延迟由拟合参数 $d$ 周刻画，不由数据直接给出。 |
| 3 | **36 天无覆盖** | 2025 年第 1–5 周与第 53 周（共 36 天）无市级周报覆盖，记为空（**不可得，非零**），已从分析输入中排除。 |
| 4 | ⚠ **气象为栅格值，非点值** | 请求广州市中心 (23.1291°N, 113.2644°E)，**ERA5 实际返回格点 (23.25°N, 113.25°E，海拔 17 m)**，偏移约 0.12°（≈13 km）。⇒ 论文方法须声明"使用 ERA5 栅格再分析，非站点实测"。 |
| 5 | ⚠ **再分析产品非实测** | ERA5 温度/降水/湿度为模式再分析，存在已知系统性偏差；本项目未做站点订正。 |
| 6 | **本地与输入分列** | 本地 767 例、输入 78 例，合计 845 例。模型拟合的本地新发序列（口径 767）与输入病例（78，作为外部输入项）**用途不同**，不得混为一谈。 |

---

## 4. 上游已有的质量保证（承接，非本项目新增）

- **`EXPECTED_GUANGZHOU_CASES` 独立人工核验表**：第 6–52 周逐周的期望 (本地, 输入) 值，
  由人工独立阅读官方周报文字得出；解析值与之不符时**构建中止**（而非静默写入坏数据）。
- **`validate_and_write`**：写盘前校验总数、覆盖天数、行数一致性。
- **`validate_report_dates`**：校验周报文字中的日期与周编号一致。
- 上游 `dataset_validation.json` 的结论为 `passed_with_expected_case_coverage_gaps`（"存在预期内的覆盖缺口"）。

---

## 5. 未做的事（诚实登记）

| 项 | 状态 |
|---|---|
| 重新抓取广东 CDC 周报以核验上游解析 | **未做**（网络抓取未经授权重跑；上游已提供 47/47 哈希与人工核验表） |
| 站点气象实测数据订正 | **未做** |
| 发病日期（onset date）重建 | **未做**（源数据无此字段） |
| 2014 年数据 | **不纳入**（研究者 2026-10-05 裁定"不补多年度验证"，义务 O9） |
| 基孔肯雅热病例数据 | **不采集**（研究者 2026-10-05 裁定"不进一步取证"，义务 O8） |
| SI / 补充材料 | 未获取（CDP 代理 127.0.0.1:3456 仍 502；阻塞项 B3） |

---

## 6. 复现命令

```bash
# 1) 导入原始数据（复制 + 哈希）
python scripts/preprocess/t06_import_raw.py

# 2) 生成分析输入（含一致性校验）
python scripts/preprocess/t06_build_analysis_input.py
```

Python 解释器：`C:\Users\HuGuoQing\.workbuddy\binaries\python\envs\default\Scripts\python.exe`
（仅用标准库：`csv` / `hashlib` / `shutil` / `json` / `pathlib`，无第三方依赖）
