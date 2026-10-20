# 分岔参数扫描（logit 参考解）
# 扫描1：β ∈ [3e-4, 5e-3]（log 均匀），θs=70, θv=100 → 终态 I vs β（检测 R0=1 与稳态类型切换）
# 扫描2：θs ∈ [30, 130]，β=0.003，θv=100 → 终态 I vs θs（检测"阈值钉扎"段 I=θs）
const DATADIR = normpath(joinpath(@__DIR__, "..", "data", "original"))
mkpath(DATADIR)
include(joinpath(@__DIR__, "model.jl"))
cd(DATADIR)
using .SVIRRE, DifferentialEquations, Printf, DelimitedFiles

sigmoid(z) = 1.0 / (1.0 + exp(-z))

function svir_logit!(du, u, p, t)
    S, V, I, R, zx, zy = u
    x = sigmoid(zx); y = sigmoid(zy)
    du[1] = (1 - p.pv) * p.A - p.β * (1 - x) * S * I - p.d * S
    du[2] = p.pv * p.A - (1 - p.pb) * p.β * (1 - y) * V * I - p.d * V
    du[3] = p.β * (1 - x) * S * I + (1 - p.pb) * p.β * (1 - y) * V * I - p.γ * I - p.d * I
    du[4] = p.γ * I - p.d * R
    du[5] = ((p.Ts0 + p.Ss0) * x - p.Ss0) * (I - p.θs)
    du[6] = ((p.Tv0 + p.Sv0) * y - p.Sv0) * (I - p.θv)
    nothing
end

# 单点求解：返回终态分类与 I 终值
function terminal_state(β, θs, θv; T=2000.0)
    p = Params(β=β, θs=θs, θv=θv)
    u0 = [103.0, 897.0, 10.0, 0.0, log(0.3/0.7), log(0.2/0.8)]
    prob = ODEProblem(svir_logit!, u0, (0.0, T), p)
    sol = solve(prob, Vern9(); abstol=1e-11, reltol=1e-11, saveat=T/4, maxiters=1e8)
    S, V, I, R, zx, zy = sol.u[end]
    x = sigmoid(zx); y = sigmoid(zy)
    # 分类
    typ = if I < 0.05
        :eliminated          # 疾病消除
    elseif abs(I - θs) < 0.05
        :pinned_s            # 钉扎 θs（P51/P52 类）
    elseif abs(I - θv) < 0.05
        :pinned_v            # 钉扎 θv（P41/P42 类）
    else
        :endemic             # 一般地方病（P112/P122/P212 类）
    end
    (I=I, x=x, y=y, typ=typ, R0=SVIRRE.R0(p))
end

# ---------- 扫描1：β ----------
println("== 扫描1: β 扫描 (θs=70, θv=100) ==")
betas = exp.(range(log(3e-4), log(5e-3), length=48))
res1 = []
for β in betas
    r = terminal_state(β, 70.0, 100.0)
    push!(res1, (β=β, r...))
    @printf("  β=%.5f R0=%7.4f  I=%9.4f  x=%.4f y=%.4f  %s\n", β, r.R0, r.I, r.x, r.y, String(r.typ))
end
open("bif_beta.csv", "w") do f
    println(f, "beta,R0,I_end,x_end,y_end,type")
    for r in res1; println(f, "$(r[1]),$(r[6]),$(r[2]),$(r[3]),$(r[4]),$(r[5])"); end
end

# ---------- 扫描2：θs ----------
println("\n== 扫描2: θs 扫描 (β=0.003, θv=100) ==")
θss = range(30.0, 130.0, length=48)
res2 = []
for θs in θss
    r = terminal_state(0.003, θs, 100.0)
    push!(res2, (θs=θs, r...))
    @printf("  θs=%6.1f  I=%9.4f  x=%.4f y=%.4f  %s\n", θs, r.I, r.x, r.y, String(r.typ))
end
open("bif_thetas.csv", "w") do f
    println(f, "thetas,R0,I_end,x_end,y_end,type")
    for r in res2; println(f, "$(r[1]),$(r[6]),$(r[2]),$(r[3]),$(r[4]),$(r[5])"); end
end

println("\n已保存 bif_beta.csv / bif_thetas.csv")
