# config/experiments/ —— 数值实验配置（T10）

本目录是 `scripts/reproduce.m` 的**唯一**配置来源。脚本只读配置，**不改配置**
（`protocols/PROJECT_RULES.md` §4.3）。每次运行把实际使用的配置副本写入
`results/runs/<run_id>/config.json`，并为该文件内容取 SHA-256 前 8 位作为 `run_id` 的 `cfghash`。

## 0. 文件清单

| 文件 | 类型 | 对应任务 | 状态 |
|---|---|---|---|
| `model_mwspm_gen3_rhohat.json` | **模型卡**（被各实验引用，不是实验） | T07/T08 | 登记完成 |
| `exp_e0_smoke_baseline.json` | 管线冒烟检查（**不承载主张**） | T10.check | READY / IMPLEMENTED（已跑通；2026-10-07 复跑锚点逐位一致） |
| `exp_c5_integrator_grid.json` | 实验 | T11.01（CL-06/CL-07） | READY / IMPLEMENTED（**已跑**：`20261006T011600_2a8ab00f`；**2026-10-07 增补第三条模型 `mwspm_gen3_t1201`**（T12.01 全参数含季节窗，后续项 P1）⇒ 重跑 `20261007T194955_f12fe263`；`models[]` 现支持可选 `onset_doy/offset_doy/gamma/beta/rho`） |
| `exp_c5_rossmacdonald.json` | 实验 | T11.02（CL-08） | READY / IMPLEMENTED（**已跑**：`20261006T012855_1c1e8c7b`，权威版） |
| `exp_c4_rho_fit.json` | 实验 | T12.01（CL-12/CL-03） | READY / IMPLEMENTED（**已跑**：`20261006T012933_c2b62344`；G0=S3b 已解禁。⚠ **2026-10-07 清理了一个大小写重复的死键** `maxiter`（见决策 `D2026-10-07-02`）⇒ 本文件 sha256 改变，未来运行的 cfghash 为 `9b0254c5`（原 `c2b62344`），**语义未变**） |
| `exp_c4_gate_g1.json` | 实验（闸门） | T13.01（CL-01/CL-02） | READY / IMPLEMENTED（**已跑**：`20261006T020359_3295cf28`，判定 GO） |
| `exp_c4_branching.json` | 实验（随机） | T13.02（CL-01/CL-09） | READY / IMPLEMENTED（**已跑**：`20261006T103456_a29b2ba7` + `20261006T105433_a29b2ba7`） |
| `exp_c4_uncertainty.json` | 实验 | T14.01（CL-03） | READY / IMPLEMENTED（**已跑**：`20261006T104119_0570febb`） |
| `exp_t14_trait_robustness.json` | 实验（敏感性） | **T14.02**（§5-4；不承载主张） | READY / IMPLEMENTED（**已跑**：`20261007T155558_86f11a38`） |
| `exp_t14_identifiability.json` | 实验（诊断） | **T14.03**（§5-5；不承载主张） | READY / IMPLEMENTED（**已跑**：`20261007T162039_d1f013a8`；见 `results/summaries/main_results.md` §4.6） |
| `exp_t14_identifiability_fine.json` | 实验（诊断，**后续项 P4**） | **T14.03**（细网格重跑；不承载主张） | READY / IMPLEMENTED（**已跑**：`20261007T195150_93456b81`；与上一行**共用同一入口函数** `exp_t14_identifiability`，差别为 `parts=["C","D"]` 与更细网格；**先声明后运行**，见 `main_results.md` §4.6.1） |
| `exp_t14_trait_source_swap.json` | 实验（敏感性，**T14.04**） | **T14.04**（§5-4 来源替换；不承载主张） | READY / IMPLEMENTED（**已跑**：`20261007T205217_442e99f5`；同源换物种：Table A 白纹伊蚊 vs **Table B 埃及伊蚊**，两种口径并列；真正的第二篇文献来源仍 `NOT PERFORMED`。见 `main_results.md` §4.7） |
| `exp_t14_trait_source_swap_c1.json` | 实验（敏感性，**T14.05**） | **T14.05**（§5-4 **跨文献部分替换**；不承载主张） | READY / IMPLEMENTED（**已跑**：`20261007T210949_913918bd`；$a/\\sigma_v/\\mu$ 取自 Liu-Helmersson 2014、$b/c$ 保留 Mordecai ⇒ **混合来源**；用入口的 `opts.alt` 具名形式表达不同函数形式。见 `main_results.md` §4.8） |
| `exp_n1_sensitivity_sigma_h_pi.json` | 实验（敏感性，**N1**） | **N1**（R2-M6 ①π ②σ_h；不承载主张） | READY / **IMPLEMENTED**（**已预注册** `D2026-10-10-01` + **已跑** `20261010T163034_1e6494d5`（SUCCESS）；结论 = **τ\* 对 σ_h、π 均不敏感**（σ_h lag1 极差 1 周、π lag1 极差 0 周）） |

两个状态字段**正交**，必须同时读：

| 字段 | 取值 | 含义 |
|---|---|---|
| `status` | `READY` \| `TODO` \| `BLOCKED` | **授权/闸门**状态。`BLOCKED` 表示被 `PROJECT_CHECKLIST.md` 的阻塞项（B8/G0/G1 等）挡住 |
| `implementation` | `IMPLEMENTED` \| `TODO` | **代码**状态。`TODO` 表示 `entry_function` 尚未在 `src/matlab/` 落地 |

因此 `READY / TODO` 的读法是"**方法上已放行、代码还没写**"，不等同于可以直接跑。
`TODO`/`BLOCKED` 的配置**可以存在并可被校验**，但 `reproduce` 在未获授权时**拒绝执行**
（默认 `require_ready=true`）。授权运行需显式 `reproduce(id, 'require_ready', false)` 并在
运行记录里说明授权来源。

## 1. 通用字段（每个实验文件都必须具备）

| 字段 | 语义 |
|---|---|
| `experiment_id` | 与文件名（去 `.json`）一致；也是 `reproduce('<id>')` 的入口名 |
| `entry_function` | `src/matlab/` 下的 MATLAB 函数名（**实验入口**，与模型右端函数分离） |
| `status` / `implementation` | 授权状态 / 代码状态，两个字段正交，定义见 §0 |
| `task` | `PROJECT_CHECKLIST.md` 的任务 ID |
| `claims` | 该实验要支持的 `protocols/CLAIMS.csv` 中的 `claim_id` 列表；冒烟检查为空数组 |
| `model` | 引用的模型卡文件名（去 `.json`） |
| `model_version` | 必须与模型卡一致；**副本**写在实验文件里以便审计时不解析模型卡即可核对 |
| `purpose` | 一句话说明这个实验回答什么 |
| `design` | 网格/情景的显式展开（初值、参数域、终止时间、求解器） |
| `initial_conditions` | 7+1 状态初值；`null` 表示取模型卡默认 |
| `parameter_domain` | **被扫描/被拟合**的参数及其取值集合或区间 |
| `termination` | `t0_date`、`t_end_days`、`t_end_date`、`events`（事件定义与方向） |
| `solvers` | 求解器条目数组（见 §2） |
| `metrics` | 指标名列表（必须在 `src/matlab/` 中有实现） |
| `controls` | 对照/参照（参考解、现口径、解析式、外部证据） |
| `outputs` | 计划写入 `results/runs/<run_id>/data/` 的文件名 |
| `deterministic` | `true` 表示确定性 ODE，**不得用重复种子伪造置信区间** |
| `rng` | 仅随机实验需要：`{"generator": "...", "seed": ...}` |
| `honesty` | 该实验的措辞约束（禁语/不可推断），与 `CLAIMS.csv` 的"禁语"列对应 |
| `blocked_by` | 阻塞原因（对应 `PROJECT_CHECKLIST.md` 的阻塞项编号） |

## 2. 求解器条目（`solvers[]`）

两种**内部步长语义**必须显式区分——这是本项目对"不要把绘图时间点当作固定内部积分步长"的落地方式：

| `kind` | 含义 | 记录的字段 |
|---|---|---|
| `adaptive` | MATLAB 自适应求解器。积分区间只给 `[t0,tf]`（2 元素），**内部网格由求解器决定**；报告网格用 `deval` 事后取值 | `method`、`RelTol`、`AbsTol`、`NonNegative`、内部网格 `n_steps/min_step/max_step`、`termination_reason` |
| `fixed_step` | 自实现定步长（`euler`/`rk4`/`tau_leap`）。**`dt_days` 就是内部步长** | `method`、`dt_days`、`clamp_nonnegative`、`termination_reason` |

- `method` 取值：`auto` | `ode15s` | `ode23s` | `ode23t` | `ode23tb` | `ode45` | `ode113` | `ode78` | `ode89` | `euler` | `rk4` | `tau_leap`。
- ⚠ **MATLAB 没有 LSODA**：`EXPERIMENT_PLAN.md` 写的 "LSODA" 在本项目**由 `ode15s` 实现**，
  并在 `manifest.json` 的 `solver_realisation` 字段中显式登记该替代（不得缩写为"已用 LSODA"）。
- `auto` 的选择规则见 `src/matlab/mwspm_pick_solver.m`：先用有限差分 Jacobian 估计刚性比
  $\max|\lambda|/\min|\lambda|$，比值超过阈值则选刚性求解器（默认 `ode15s`），否则 `ode45`；
  判据与实测比值一并写入 manifest。
- 报告网格（`report_grid`）只影响**输出采样**：自适应求解器不因报告密而减小内部步。

## 3. 事件（`termination.events[]`）

```
{ "name": "...", "kind": "zero_crossing" | "threshold", "state_index": 8,
  "direction": 1, "threshold": 1.0, "description": "..." }
```

- `zero_crossing` → 传给 `odeset('Events', ...)`，**只有自适应求解器支持**；
- `threshold` → 由 `fixed_step` 积分器在步末检查，触发即停止并记 `termination_reason = 'event:<name>'`；
- 未触发时 `termination_reason = 't_end_reached'`。

## 4. run_id 与产物

```
run_id = <YYYYmmddTHHMMSS 本地时间 UTC+8>_<cfg_sha256 前 8 位>
results/runs/<run_id>/
├─ config.json      实际使用的配置副本（含解析后的模型参数）
├─ manifest.json    status(SUCCESS|FAILED) / env / 代码哈希 / 输入哈希 / 求解器与终止原因 / 输出清单 / rng
├─ logs/            运行日志与警告
└─ data/            数值输出（时间网格、时间序列、汇总表）
```

- 输入哈希 = 模型卡 + 数据文件 + 本实验配置的 SHA-256；代码哈希 = 实际载入的 `.m` 文件 SHA-256。
- **失败不留成功标志**：`manifest.json` 的 `status` 写 `FAILED` 并记录错误消息与栈。
- 确定性实验**不写** `ci95` 之类的伪区间；随机实验必须写 `rng.generator` 与 `rng.seed`。

## 5. 图表

`scripts/plots/` 只读 `results/runs/<run_id>/data/` 下已保存的数值，**不重算模型**
（`PROJECT_RULES.md` §4.6）。本步（T10）不制作最终论文图。
