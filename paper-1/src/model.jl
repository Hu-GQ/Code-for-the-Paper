# SVIR-RE 模型定义与 ODE 求解
# 对应论文 "Research on SVIR epidemic model with evolutionary game of inoculation" 公式 (4)
#
# dS/dt  = (1-p_v)A - β(1-x_A)SI - dS
# dV/dt  = p_v A - (1-p_b)β(1-y_A)VI - dV
# dI/dt  = β(1-x_A)SI + (1-p_b)β(1-y_A)VI - γI - dI
# dR/dt  = γI - dR
# dx/dt  = x(1-x)[(T_S0+S_S0)x - S_S0](I - θ_s)
# dy/dt  = y(1-y)[(T_V0+S_V0)y - S_V0](I - θ_v)

module SVIRRE

export Params, svir_ode!, run_scenario

using DifferentialEquations

struct Params
    A::Float64
    pv::Float64
    pb::Float64
    γ::Float64
    d::Float64
    β::Float64
    Ts0::Float64
    Ss0::Float64
    Tv0::Float64
    Sv0::Float64
    θs::Float64
    θv::Float64
end

# 论文 table:2 的基准参数（β、θs、θv 各场景单独给定）
# p_v=0.897（原稿 table:2 的 0.879 系笔误；修订版论文与 Table 2 现均印 0.897）。
# 注意：本模块为直接 (x,y) 变量的辅助实现，run_scenario 用 Tsit5/1e-9，仅适合快速探索；
# 论文全部图表的生产管线是 recompute_all.jl（logit 坐标 + Vern9 + 1e-11），二者 RHS 等价。
function Params(; A=10.0, pv=0.897, pb=0.7934, γ=0.143, d=0.007,
                β=0.003, Ts0=1.5, Ss0=-2.2, Tv0=1.3, Sv0=-2.0,
                θs=70.0, θv=100.0)
    Params(A, pv, pb, γ, d, β, Ts0, Ss0, Tv0, Sv0, θs, θv)
end

# 基本再生数 R_0 = β A (1 - p_b p_v) / (d(γ+d))
R0(p::Params) = p.β * p.A * (1 - p.pb * p.pv) / (p.d * (p.γ + p.d))

# ODE 右端函数（就地版本）
function svir_ode!(du, u, p, t)
    S, V, I, R, x, y = u
    Ts0, Ss0 = p.Ts0, p.Ss0
    Tv0, Sv0 = p.Tv0, p.Sv0

    du[1] = (1 - p.pv) * p.A - p.β * (1 - x) * S * I - p.d * S
    du[2] = p.pv * p.A - (1 - p.pb) * p.β * (1 - y) * V * I - p.d * V
    du[3] = p.β * (1 - x) * S * I + (1 - p.pb) * p.β * (1 - y) * V * I - p.γ * I - p.d * I
    du[4] = p.γ * I - p.d * R
    du[5] = x * (1 - x) * ((Ts0 + Ss0) * x - Ss0) * (I - p.θs)
    du[6] = y * (1 - y) * ((Tv0 + Sv0) * y - Sv0) * (I - p.θv)
    nothing
end

# 运行一个场景，返回解对象
function run_scenario(p::Params; u0=[103.0, 897.0, 10.0, 0.0, 0.3, 0.2],
                       tspan=(0.0, 5000.0), saveat=1.0, solver=Tsit5())
    prob = ODEProblem(svir_ode!, u0, tspan, p)
    solve(prob, solver; saveat=saveat, abstol=1e-9, reltol=1e-9, maxiters=1e7)
end

end # module
