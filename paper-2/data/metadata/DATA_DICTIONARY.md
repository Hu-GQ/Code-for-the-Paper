# 数据字典 DATA_DICTIONARY

- 项目：传染病动力学期刊论文（广州 2025 登革热人–蚊传播模型）
- 版本：v1.0｜建立：2026-10-05（T06）
- 配套：`data/metadata/sources.csv`、`data/metadata/parameters.csv`、
  `data/metadata/raw_file_hashes.csv`、`data/metadata/PROCESSING.md`

> **单位一律写在列名中或本表"单位"列，不做隐式约定。**
> 所有 CSV 均为 UTF-8 with BOM（`utf-8-sig`）、CRLF 由 `csv.writer` 默认产生、字段全量加引号。

---

## 1. 周编号：两套口径（🔴 最易错位，务必先看）

| 口径 | 记号 | 范围 | 起点 | 用途 |
|---|---|---|---|---|
| **分析周序号** | `week_index` | **1 – 47** | 2025-02-03 | **本项目模型与全文的默认口径**；`config/project.json` 的 `peak_week_model=37`、`peak_week_observed=37` 均为此口径 |
| **ISO 周** | `iso_week` | **6 – 52** | 2025-02-03 | 广东 CDC 周报原始编号 |

**映射：`week_index = iso_week − 5`**（已由 `scripts/preprocess/t06_build_analysis_input.py` 校验为恒等式）。

| 关键事件 | `week_index` | `iso_week` | 周起始 | 本地报告 |
|---|---|---|---|---|
| 分析窗口首周 | 1 | 6 | 2025-02-03 | 0 |
| **首个非零本地病例** | **19** | 24 | 2025-06-09 | 1 |
| 观测双峰之一 | 28 | 33 | 2025-08-11 | 59 |
| **观测主峰** | **37** | 42 | 2025-10-13 | **78** |
| 分析窗口末周 | 47 | 52 | 2025-12-22 | 0 |

⚠ **写作约束**：稿件中任何"第 X 周"**必须**声明口径；建议正文统一用 `week_index`，首次出现时注明
"（`week_index = 1` 对应 2025 年 2 月 3 日起的一周，即 ISO 第 6 周）"。

---

## 2. `data/processed/weekly_analysis_input.csv`（47 行 × 11 列）—— 正式分析输入

| 列名 | 类型 | 单位 | 含义 | 缺失约定 |
|---|---|---|---|---|
| `week_index` | int | 无量纲（周） | 分析周序号 1–47 | 无 |
| `iso_week` | int | 无量纲（周） | ISO 周编号 6–52 | 无 |
| `week_start` | date (ISO 8601) | — | 周起始日（周一） | 无 |
| `week_end` | date (ISO 8601) | — | 周结束日（周日） | 无 |
| `local_reported` | int | **例（人）** | 广州**本地**报告登革热病例数（按现住址） | 0 = 真实零；无缺失 |
| `imported_reported` | int | **例（人）** | 广州**输入**报告病例数 | 0 = 真实零；无缺失 |
| `total_reported` | int | 例（人） | 本地 + 输入 | 派生列，无缺失 |
| `temperature_mean_c` | float | **摄氏度 °C** | 该周日均温（ERA5 逐日聚合） | 无缺失（气象覆盖 365 天） |
| `precipitation_sum_mm` | float | **毫米 mm** | 该周降水总量 | 无缺失 |
| `relative_humidity_mean_pct` | float | **百分比 %** | 该周平均相对湿度（逐小时→日均→周均） | 无缺失 |
| `source_url` | str | — | 该周广东 CDC 周报原文 URL | 无 |

**internal一致性（由构建脚本校验，见 `data/metadata/processed_checks.json`）**：

- 行数 47 = `official_weekly_case_rows` ✅
- `week_index` 连续 1–47 ✅；`iso_week − week_index ≡ 5` ✅
- 本地合计 **767**、输入合计 **78**，与上游 `dataset_validation.json` 一致 ✅
- 峰值 `week_index=37`（78 例），首个非零 `week_index=19` ✅

**观测过程的性质（关键）**：本表是**新增报告病例数**，按**报告日期**归周，
不是现患人数、不是住院数、不是累计数。与模型状态的关系由观测模型给出：

$$\text{报告数}(w) \;\approx\; \rho \times (\text{该周真实新发感染数}) \;\text{并延迟}\ d\ \text{周}$$

其中 $\rho$、$d$ 为待拟合参数（见 `parameters.csv`），**不作外部文献取值**。

---

## 3. `data/raw/upstream_cleaned/`（8 个上游产物，原样复制 · 只读）

| 文件 | 粒度 | 关键字段 | ⚠ 性质 |
|---|---|---|---|
| `guangzhou_2025_weekly_cases_official.csv` | 周（47） | `week, week_start, week_end, local_reported, imported_reported, source_url, data_status` | **官方一手**；本项目病例数据的唯一权威来源 |
| `guangzhou_2025_daily_cases_reconstructed.csv` | 日（365） | `date, local_reported_reconstructed, imported_reported_reconstructed, reconstruction_status` | 🔴 **日尺度是重构的，不是真实日数据**——由周计数在周内**均匀分配**得到（`allocate_uniform`）。**仅用于月度聚合与绘图，不得用于日尺度推断** |
| `guangzhou_2025_monthly_cases_observed.csv` | 月（12） | `month, local_reported, imported_reported, *data_status, source` | 由重构日表聚合；2 月仅覆盖 26/28 天 |
| `guangzhou_2025_weekly_model_data.csv` | 周（53） | `week_start, week_end, *_reconstructed, temperature_mean_c, precipitation_sum_mm, relative_humidity_mean_pct, data_status` | 53 周含气象全-year；**其中 47 周有病例**，另 6 周病例为空 |
| `guangzhou_2025_daily_weather_openmeteo.csv` | 日（365） | `date, temperature_mean_c, temperature_min_c, temperature_max_c, precipitation_mm, relative_humidity_mean_pct, source` | ERA5 再分析 |
| `guangzhou_2025-2026_daily_weather_extension_openmeteo.csv` | 日（63） | `date, temperature_mean_c, source` | 2025-12-29 起的延伸段，本项目分析窗口（至 2025-12-28）**不使用** |
| `guangzhou_population.csv` | 年（1） | `year, population, unit, source, status` | 19,101,000 人 |
| `dataset_validation.json` | — | 见 §4 | 上游自校验 |

---

## 4. `dataset_validation.json` 关键字段（上游自校验，本项目已独立复算）

| 字段 | 值 | 本项目独立复算 |
|---|---|---|
| `local_reported_total` | 767 | **767 ✅** |
| `imported_reported_total` | 78 | **78 ✅** |
| `case_coverage_start` / `end` | 2025-02-03 / 2025-12-28 | ✅ |
| `case_covered_days` | 329 | ✅（47 周 × 7 天） |
| `case_missing_days` | 36 | ✅（365 − 329） |
| `official_weekly_case_rows` | 47 | ✅ |
| `weekly_rows` | 53 | ✅ |
| `weather_days` | 365 | ✅ |
| `population` | 19,101,000 | ✅ |
| `validation` | `passed_with_expected_case_coverage_gaps` | ✅ |

---

## 5. `data/metadata/raw_file_hashes.csv`（58 行）

| 列名 | 含义 |
|---|---|
| `target_path` | 本项目内的路径（`data/raw/...`） |
| `level` | `L1_原始`（50 个：47 份周报 HTML + manifest + ERA5 JSON + 人口 HTML）／`L2_上游清洗产物`（8 个） |
| `source_path` | 源项目只读路径（绝对路径，仅作溯源记录） |
| `bytes` | 字节数 |
| `sha256` | 文件 SHA-256 |
| `status` | `copied` / `already_present_identical` |
| `note` | 说明 |

**溯源链交叉核验**：47 份周报 HTML 的 SHA-256 与上游 `source_manifest.json` 登记值
**逐条一致（47/47）** ⇒ 复制未改变任何字节。
