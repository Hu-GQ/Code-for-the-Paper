# Generate baseline SVIR and VIR trajectories in exact logit coordinates.
using DifferentialEquations, DelimitedFiles, Printf

# Resolve the data directory relative to this file, so the script can be run from any cwd.
const DATADIR = normpath(joinpath(@__DIR__, "..", "data", "original"))
mkpath(DATADIR)
cd(DATADIR)

sigmoid(z) = 1.0/(1.0+exp(-z))

# ---------- SVIR-RE logit 系统（变量 S,V,I,R,zx,zy） ----------
function svir_logit!(du, u, p, t)
    A,pv,pb,γ,d,β,Ts0,Ss0,Tv0,Sv0,θs,θv = p
    S,V,I,R,zx,zy = u
    x = sigmoid(zx); y = sigmoid(zy)
    du[1] = (1-pv)*A - β*(1-x)*S*I - d*S
    du[2] = pv*A - (1-pb)*β*(1-y)*V*I - d*V
    du[3] = β*(1-x)*S*I + (1-pb)*β*(1-y)*V*I - (γ+d)*I
    du[4] = γ*I - d*R
    du[5] = ((Ts0+Ss0)*x - Ss0)*(I-θs)
    du[6] = ((Tv0+Sv0)*y - Sv0)*(I-θv)
    nothing
end
# Baseline recruitment vaccination fraction: p_v = 0.897.
P0 = [10.0, 0.897, 0.7934, 0.143, 0.007, 0.003, 1.5, -2.2, 1.3, -2.0, 70.0, 100.0]

function run_svir(p, u0; tspan=(0.0,400.0), saveat=0.5)
    prob = ODEProblem(svir_logit!, u0, tspan, p)
    solve(prob, Vern9(); abstol=1e-11, reltol=1e-11, saveat=saveat, maxiters=1e8)
end
u0 = [103.0, 897.0, 10.0, 0.0, log(0.3/0.7), log(0.2/0.8)]

# ---------- VIR-RE logit 系统（变量 V,I,R,zy） ----------
function vir_logit!(du, u, p, t)
    A,pb,γ,d,β,Tv0,Sv0,θv = p
    V,I,R,zy = u
    y = sigmoid(zy)
    du[1] = A - (1-pb)*β*(1-y)*V*I - d*V
    du[2] = (1-pb)*β*(1-y)*V*I - (γ+d)*I
    du[3] = γ*I - d*R
    du[4] = ((Tv0+Sv0)*y - Sv0)*(I-θv)
    nothing
end
function run_vir(p, u0; tspan=(0.0,600.0), saveat=0.5)
    prob = ODEProblem(vir_logit!, u0, tspan, p)
    solve(prob, Vern9(); abstol=1e-11, reltol=1e-11, saveat=saveat, maxiters=1e8)
end

# ================= 场景 1：P110（疾病消除） =================
p = copy(P0); p[6]=0.00032; p[11]=30.0; p[12]=40.0   # β, θs, θv
s = run_svir(p, u0; tspan=(0.0,3000.0))
writedlm("rev_P110.csv", hcat(s.t, [u[1] for u in s.u],[u[2] for u in s.u],
    [u[3] for u in s.u],[u[4] for u in s.u],[sigmoid(u[5]) for u in s.u],[sigmoid(u[6]) for u in s.u]), ',')
println("P110 done: I_end=", s.u[end][3])

# ================= 场景 2：P112（地方病，完全无合作） =================
p = copy(P0); p[6]=0.003; p[11]=70.0; p[12]=100.0
s = run_svir(p, u0; tspan=(0.0,400.0))
writedlm("rev_P112.csv", hcat(s.t, [u[1] for u in s.u],[u[2] for u in s.u],
    [u[3] for u in s.u],[u[4] for u in s.u],[sigmoid(u[5]) for u in s.u],[sigmoid(u[6]) for u in s.u]), ',')
println("P112 done: I_end=", s.u[end][3], " x_end=", sigmoid(s.u[end][5]))

# ================= 场景 3：P51（钉扎 θs，部分合作） =================
p = copy(P0); p[6]=0.003; p[11]=56.0; p[12]=76.0
s = run_svir(p, u0; tspan=(0.0,5000.0))
writedlm("rev_P51.csv", hcat(s.t, [u[1] for u in s.u],[u[2] for u in s.u],
    [u[3] for u in s.u],[u[4] for u in s.u],[sigmoid(u[5]) for u in s.u],[sigmoid(u[6]) for u in s.u]), ',')
println("P51 done: I_end=", s.u[end][3], " x_end=", sigmoid(s.u[end][5]))

# ================= 图4：θs/θv 影响（SVIR-RE） =================
# (a) disease-free: A=3, pv=0.98, pb=0.9, γ=0.2, β=0.003
pbase_a = [3.0, 0.98, 0.9, 0.2, 0.007, 0.003, 1.5, -2.2, 1.3, -2.0, 30.0, 40.0]
for (i,θs) in enumerate([10.0, 30.0, 50.0])
    p = copy(pbase_a); p[11]=θs; p[12]=θs+10
    s = run_svir(p, u0; tspan=(0.0,2000.0))
    writedlm("rev_fig4a_$i.csv", hcat(s.t, [u[3] for u in s.u]), ',')
end
# (b) endemic: A=10, pv=0.897, pb=0.7934, γ=0.143, β=0.003
pbase_b = copy(P0)
for (i,θs) in enumerate([40.0, 56.0, 70.0])
    p = copy(pbase_b); p[11]=θs; p[12]=θs+20
    s = run_svir(p, u0; tspan=(0.0,5000.0))
    writedlm("rev_fig4b_$i.csv", hcat(s.t, [u[3] for u in s.u]), ',')
end
println("fig4 done")

# ================= 图5：VIR-RE θv 影响 =================
# (a) disease-free: A=7, pb=0.9, γ=0.3, β=0.003
pb_a = [7.0, 0.9, 0.3, 0.007, 0.003, 1.3, -2.0, 40.0]
u0v = [1000.0, 10.0, 0.0, log(0.2/0.8)]
for (i,θv) in enumerate([20.0, 50.0, 80.0])
    p = copy(pb_a); p[8]=θv
    s = run_vir(p, u0v; tspan=(0.0,800.0))
    writedlm("rev_fig5a_$i.csv", hcat(s.t, [u[2] for u in s.u]), ',')
end
# (b) endemic: A=10, pb=0.7934, γ=0.143, β=0.003
pb_b = [10.0, 0.7934, 0.143, 0.007, 0.003, 1.3, -2.0, 76.0]
for (i,θv) in enumerate([19.0, 38.0, 57.0, 76.0, 95.0])
    p = copy(pb_b); p[8]=θv
    s = run_vir(p, u0v; tspan=(0.0,20000.0), saveat=1.0)
    writedlm("rev_fig5b_$i.csv", hcat(s.t, [u[2] for u in s.u]), ',')
end
println("fig5 done")
