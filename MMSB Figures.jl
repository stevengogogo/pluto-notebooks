### A Pluto.jl notebook ###
# v0.14.5

using Markdown
using InteractiveUtils

# ╔═╡ 287b9e43-2313-4724-b873-d9bd87a8948d
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="Plots", version="1"),
        Pkg.PackageSpec(name="PlutoUI", version="0.7"),
        Pkg.PackageSpec(name="PyPlot", version="2"),
		Pkg.PackageSpec(name="DifferentialEquations", version="6"),
		Pkg.PackageSpec(name="Parameters", version="0.12"),
        Pkg.PackageSpec(name="LabelledArrays", version="1"),
        Pkg.PackageSpec(name="Setfield", version="0.7"),
    ])
    using Plots, PlutoUI, LinearAlgebra, DifferentialEquations, Parameters, LabelledArrays, Setfield
	
	import PyPlot as plt
	Plots.gr(fmt=:png, lw=2)
	
	TableOfContents()
end

# ╔═╡ 5ba678b7-7251-4003-9f6e-df6339335635
md"""
# Fig 1.07 Collins toggle switch

For Figures 1.7, 7.13, 7.14, 7.15
"""

# ╔═╡ 8fbe196a-0cb1-408f-bb89-d5a2bcdb3a40
md"""
# Fig 1.09 Hodgkin-Huxley model
"""

# ╔═╡ eb14ab5a-e7d7-415a-8cbc-7d5152c4e01b
md"""
# Fig 2.04 Exponential decay
"""

# ╔═╡ 333892b7-e248-422f-b16e-0cde638dd9ee
plot([t-> 3 * exp(-t) t->3 * exp(-2t) t-> 3 * exp(-3t)], 0.0, 5.0, xlim = (0, 5), ylim=(0, 3.2), xlabel="Time", ylabel="Concentration", label = ["exp(-t)" "exp(-2t)" "exp(-3t)"], title= "Figure 2.4")

# ╔═╡ a4e81834-0827-4d0a-a2c8-6fff6002ddc9
md"""
# Fig 2.07 Euler method
"""

# ╔═╡ 478ee763-2b39-4716-86fb-ad49194b5125
function fig0207()
	
	# Forward Euler stepper
	step_euler(u, p, t, h, f) = u + f(u, p, t) * h
	
	model(u, p, t) = p * u
	
	
	tspan = (0.0, 2.0)
	p = -2.0
	u0 = 1.0
	
	plot(x->exp(p*x), tspan[1], tspan[end], lab="True solution", legend=:topright)
	
	for h in (1//16, 1//8, 1//4, 1//2)
		ts = range(tspan[1], step=h, stop=tspan[2])
		us = fill(u0, length(ts))
		u = us[1]
		for i in 2:length(ts)
			us[i] = step_euler(us[i-1], p, ts[i-1], h, model)
		end
		plot!(ts, us, lab=string("h = ", float(h)))
	end
	
 	plot!(title="Fig 2.07 Euler method", xlabel="Time (arbitrary units)", ylabel="Concentration (arbitrary units)")
end

# ╔═╡ 68786419-44ab-48b9-b46c-135af4ffaf95
fig0207()

# ╔═╡ 805ba136-caea-4391-9e9e-b3275341613e
md"""
# Fig 2.09 Numerical Simulation of a metabolic network
"""

# ╔═╡ 76aab80e-40eb-473c-8f82-49c67c304086
function fig0209()
	function model!(du, u, p, t)
		@unpack a, b, c, d = u
		v1 = 2.0a
		v2 = 2.5*a*b
		v3 = 3.0c
		v4 = 3.0d

		du.a = 3.0 - v1 - v2
		du.b = v1 - v2
		du.c = v2 - v3
		du.d = v2 - v4
		return du
	end
	
	u0 = LVector(a=0.0, b=0.0, c=0.0, d=0.0)
	tend = 10.0
	sol = solve(ODEProblem(model!, u0, tend))
	plot(sol, xlims=(0.0, 4.0), ylims=(0.0, 1.0), 
     xlabel="Time (sec)", ylabel="Concentration (mM)", title="Figure 2.09",
     label=["A" "B" "C" "D"], legend=:bottomright)
end

# ╔═╡ 7623e4ca-d9e6-4e97-9dee-95efe616f959
fig0209()

# ╔═╡ dcff0bfd-6f4b-4ff8-8532-59888754bb52
md"""
# Figure 2.11-14 Model reduction
"""

# ╔═╡ 037facee-27ca-4a42-a79c-aa8cad6a1187
function fig0211()
	
	# Fig 2.11: Full model
	function fullmodel!(du, u, p, t)
		@unpack K0, K1, KM1, K2 = p
		@unpack a, b = u
		vAB = K1 * a - KM1 * b
		du.a = K0 - vAB
		du.b = vAB - K2 * b
		return du
	end
	
	tspan = 3.0
	p1 = (K0=0, K1=9, KM1=12, K2=2)
	u0 = LVector(a=0.0, b=10.0)
	sol1 = solve(ODEProblem(fullmodel!, u0, tspan, p1))
	
	pl1 = plot(sol1, xlabel="Time (arbitrary units)", ylabel="Concentration (arbitrary units)", title="Fig. 2.11 (Full model)", label=["A" "B"])
	
	# Fig 2.12: rapid equilibrium
	function remodel!(du, u, p, t)
		@unpack K0, K1, KM1, K2 = p
		kb = K2 * K1 / (KM1 + K1)
		du.b = K0 - kb * u.b
	end
	
	
	u0re = LVector(b=sum(u0))
	sol2 = solve(ODEProblem(remodel!, u0re, tspan, p1))
	ts = 0.0:0.1:tspan
	btilde = sol2(ts, idxs=1)

	pl2 = plot(sol1, line=(:dash, 1),label=["A (full solution)" "B (full solution)"])
	plot!(pl2, ts, (p1.KM1 / (p1.KM1 + p1.K1)) .* btilde, lab="A (rapid equilibrium)")
	plot!(pl2, ts, (p1.K1 / (p1.KM1 + p1.K1))  .* btilde, lab="B (rapid equilibrium)")
	plot!(pl2, title="Fig. 2.12 (Rapid equilibrium model)", xlabel="Time (arbitrary units)", ylabel="Concentration (arbitrary units)")
	
	# Fig 2.13: another rapid equilibrium
	p2 = (K0=9, K1=20, KM1=12, K2=2)
	u0 = LVector(a=8.0, b=4.0)
	u0re = LVector(b=sum(u0))
	sol3full = solve(ODEProblem(fullmodel!, u0, tspan, p2))
	sol3re = solve(ODEProblem(remodel!, u0re, tspan, p2))
	
	btilde = sol3re(ts, idxs=1)
	pl3 = plot(title="Fig. 2.13 Ref vs Fast equlibrium", xlabel="Time (arbitrary units)", ylabel="Concentration (arbitrary units)")
	plot!(pl3, sol3full, line=(:dash, 1),label=["A (full solution)" "B (full solution)"])
	plot!(pl3, ts, (p2.KM1 / (p2.KM1 + p2.K1)) .* btilde, lab="A (rapid equilibrium)")
	plot!(pl3, ts, (p2.K1 / (p2.KM1 + p2.K1))  .* btilde, lab="B (rapid equilibrium)")
	
	# Figure 2.14: Quasi-steady state assumption(QSSA)
	
	function qssmodel!(du, u, p, t)
		@unpack K0, K2 = p
		du.b = K0 - K2 * u.b
	end
	
	u0_qss = LVector(b = (p2.K1 * sum(u0) - p2.K0) / (p2.K1 + p2.KM1))
	sol4 = solve(ODEProblem(qssmodel!, u0_qss, tspan, p2))
	btilde = sol4(ts, idxs=1)
	
	pl4 = plot(sol3full, line=(:dash, 2), xlims=tspan,
     xlabel="Time (arbitrary units)",
     ylabel="Concentration (arbitrary units)",
     title="Figure 2.14: Ref vs QSSA")

	plot!(pl4, sol4, label="B (QSSA)", line=(2, :red))
	plot!(pl4, ts, (p2.K0 .+ p2.KM1 .* btilde) ./ p2.K1, label="A (QSSA)", line=(2, :blue))
	
	return (pl1, pl2, pl3, pl4)
end

# ╔═╡ 1a4984d2-1e8c-4f3b-94c4-fb937e9de6f1
let
	figs = fig0211()
	plot(figs..., size=(800, 800))
end

# ╔═╡ 5021a671-13d3-428e-bc92-4e3cf930fa0b
md"""
# Problem 2.4.6 Continuation diagram
"""

# ╔═╡ 0ebae16b-3b6c-4298-9812-fac0227e8ed7
function figp246()
	f(u, p, t) = p * (1.0-u)
	p = 1.0
	u0 = 0.0
	tspan = (0.0, 10.0)
	prob = ODEProblem(f, u0, tspan, p)
	sol = solve(prob)
	plot(sol, xlabel="Time (arbitrary units)", ylabel="Concentration (arbitrary units)", legend=:bottomright)
end

# ╔═╡ ca209a4f-ca44-4d0c-8c2a-83954ccd964a
figp246()

# ╔═╡ da0b9be3-36d8-4f9b-b03d-30270111e862
md"""
# Figure 3.03 Michaelis-Menten kinetics
"""

# ╔═╡ 91b7d696-109a-44b1-9631-69e514a290ac
md"""
# Problem 3.7.5 Reduced Michaelis-Menten
"""

# ╔═╡ a4a64652-a4b9-4345-9652-252daff344d8
md"""
# Figure 4.1, 4.2, and 4.3

Steady states and phase plots
"""

# ╔═╡ f3c1efd0-b76d-496a-83bd-f2c8c57d2b94
md"""
# Figure 4.4, 4.5 Vector fields
"""

# ╔═╡ 04fc83f3-1d1b-4ba7-a08c-7fc646a73de3
md"""
# Figure 4.7, 4.8, 4.9, and 4.19A

Symmetric (bistable) biological networks.
"""

# ╔═╡ 9519b1ba-726f-4a25-8993-fbe897e998cd
md"""
# Figure 4.11

Surface plots.

Reference: [PlotsGallery.jl](https://goropikari.github.io/PlotsGallery.jl/src/surface.html)
"""

# ╔═╡ 96e7a435-3411-4412-89d9-fb3b42e62a7c
let
	z1(x, y) = x^2 + 0.5y^2
	z2(x, y) = (.2x^2-1)^2 + y^2
	
	x1 = range(-1.0, 1.0, length=51)
	y1 = range(-1.0, 1.0, length=51)
	
	x2 = range(-2.75, 2.75, length=80)
	y2 = range(-0.75, 0.75, length=80)
	
	p1 = surface(x1, y1, z1, title="Single-well potential")
	p2= contourf(x1, y1, z1)
	p3 = surface(x2, y2, z2, title="Double-well potential")
	p4 = contourf(x2, y2, z2)

	plot(p1, p2, p3, p4,size=(1000, 800))
end

# ╔═╡ 7b747fe2-05bc-40b5-9fcf-f6e11125cc91
md"""
# Figure 4.15, 4.16, and 4.17

Oscillatory network.
"""

# ╔═╡ ec439af8-b24b-4167-8b79-756955e1375d
"""
Model of oscillatory network from Figure 4.14. This code generates Figures
4.15, 4.16, and 4.17
"""
function model0415!(du, u, p, t)
    @unpack s1, s2 = u
    @unpack K0, K1, K2, N = p
    v0 = K0
    v1 = K1 * s1 * (1 + s2^N)
    v2 = K2 * s2
	du.s1 = v0 - v1
	du.s2 = v1 - v2
    return du
end

# ╔═╡ 998ba19b-63ca-4079-ab4b-f6de5d0d61b0
md"""
# Figure 4.18 Continuation diagram

[Bifurcations.jl](https://github.com/tkf/Bifurcations.jl) does not work on Julia v1.6. And [BifurcationKit.jl](https://github.com/rveltz/BifurcationKit.jl) might be too complex for this example.

"""

# ╔═╡ 32825174-50e3-42b9-9c5d-ea744b36fefb
md"""
# Figure 4.22 Tangent line
"""

# ╔═╡ c94276d8-e844-437d-8ede-742bba44a2fd
let 
	curve(t) = 3 / (t-2)
	tange(t) = 1.5 - (t - 4) * 0.75

	plot(curve, 2.2, 8.0, lab="Curve")
	plot!(tange, 2.7, 5.3, lab="Tangent line")
	plot!(title="Fig 4.22", xlabel="Reaction rate", ylabel="Inhibitor concentration", 
		  xlims=(2.0, 8.0), ylims=(0.0, 4.0))
end

# ╔═╡ 6ede64c4-4851-4506-88d9-7e58d49c35e9
begin
	# Convenience functions
    hill(x, k) = x / (x + k)
    hill(x, k, n) = hill(x^n, k^n)
    exprel(x) = ifelse(x≈zero(x), one(x), x / expm1(x))
	md"""
	# Appendix
	Running environment and some auxillary functions
	"""
end

# ╔═╡ 0a52a098-9f37-4d89-96d1-813a0abc2f3b
"""
model of Collins toggle switch
from Gardiner et al. (2000) Nature 403, pp. 339-342
Figures 1.7, 7.13, 7.14, 7.15
"""
function fig0107()
	# The model
	function collins!(du, u, p, t)
		@unpack a1, a2, β, γ, i1, i2 = p
		@unpack s1, s2 = u
		du.s1 = a1 * hill(1 + i2, s2, β) - s1
		du.s2 = a2 * hill(1 + i1, s1, γ) - s2
    	return du
	end
	
	# Events and callbacks
    on_i1!(integrator) = integrator.p = @set p.i1 = 10.0
    off_i1!(integrator) = integrator.p = @set p.i1 = 0.0
    on_i2!(integrator) = integrator.p = @set p.i2 = 10.0
    off_i2!(integrator) = integrator.p = @set p.i2 = 0.0

    cbs = CallbackSet(PresetTimeCallback([10.0], on_i2!),
                  PresetTimeCallback([20.0], off_i2!),
                  PresetTimeCallback([30.0], on_i1!),
                  PresetTimeCallback([40.0], off_i1!))
	u0 = LVector(s1=0.075, s2=2.5)
    tend = 50.0
    p = (a1 = 3.0, a2 = 2.5, β = 4.0, γ = 4.0, i1 = 0.0, i2 = 0.0)
	
	sol = solve(ODEProblem(collins!, u0, tend, p), callback=cbs)
	plot(sol, legend=:right, xlabel="Time", ylabel="Concentration", title="Figure 1.7")
end

# ╔═╡ 33edf7e8-d3d6-4ced-844b-9fd2c0665293
fig0107()

# ╔═╡ 3a7e2fc4-05ca-4e80-b2b2-48a0513dfc4d
function fig0109()
	"ODE system of HH model"
	function hh!(du, u, p, t)
		@unpack E_N, E_K, E_LEAK, G_N_BAR, G_K_BAR, G_LEAK, C_M, I_STIM = p
		@unpack v, m, h, n = u

		mαV = -0.10 * (v + 35)
		mα = exprel(mαV)
		mβ = 4.0 * exp(-(v + 60) / 18.0)

		hα = 0.07 * exp(-(v+60)/20)
		hβ = 1 / (exp(-(v+30)/10) + 1)

		nαV = -0.1 * (v+50)
		nα = 0.1 * exprel(nαV)
		nβ = 0.125 * exp( -(v+60) / 80)

		iNa = G_N_BAR * (v - E_N) * (m^3) * h
		iK  = G_K_BAR * (v - E_K) * (n^4)
		iLeak = G_LEAK * (v - E_LEAK)
		du.v = -(iNa + iK + iLeak + I_STIM) / C_M
		du.m = -(mα + mβ) * m + mα
		du.h = -(hα + hβ) * h + hα
		du.n = -(nα + nβ) * n + nα
		return du
	end
	
	# Events and callbacks
	on_i1!(integrator) = integrator.p = @set p.I_STIM = -6.80
    off_i1!(integrator) = integrator.p = @set p.I_STIM = 0.0
    on_i2!(integrator) = integrator.p = @set p.I_STIM = -6.90

    cbs = CallbackSet(PresetTimeCallback([20.0], on_i1!),
                  PresetTimeCallback([21.0, 61.0], off_i1!),
                  PresetTimeCallback([60.0], on_i2!))
	
	# Parameters for Hodgkin-Huxley model
    p = (E_N = 55.0,       # Reversal potential of Na (mV)
		 E_K = -72.0,      # Reversal potential of K (mV)
		 E_LEAK = -49.0,   # Reversal potential of leaky channels (mV)
		 G_N_BAR = 120.0,  # Max. Na channel conductance (mS/cm^2)
		 G_K_BAR = 36.0,   # Max. K channel conductance (mS/cm^2)
		 G_LEAK = 0.30,    # Max. leak channel conductance (mS/cm^2)
		 C_M = 1.0,        # membrane capacitance (uF/cm^2)
		 I_STIM = 0.0     # Stimulus current (uA/cm^2)
    )
	
	
	# Initial conditions, time span, and parameters
	u0 = LVector(v=-59.8977, m=0.0536, h=0.5925, n=0.3192)
	tend = 100.0
	sol = solve(ODEProblem(hh!, u0, tend, p), callback = cbs)
	plot(sol, vars=(0, 1), ylabel="Membrane potential (mV)", xlabel="Time (ms)", legend=false, title="Figure 1.9")
end

# ╔═╡ f0fc6b02-c070-49d8-bce3-777758c6d482
fig0109()

# ╔═╡ 58d5e228-be35-4b1f-b9f2-368e9c99ee20
function fig0303(; 	tend = 1.0, 
				   	p = (ET = 1.0, K1 = 30.0, KM1 = 1.0, K2 = 10.0),
					u0 = LVector(S= 5.0, ES = 0.0, P = 0.0))
	
	
	# Full enzyme catalysis model
	function full_model!(du, u, p, t)
		@unpack ET, K1, KM1, K2 = p
		@unpack S, ES, P = u
		freeEnz = ET - ES
		v1 = K1 * S * freeEnz - KM1 * ES
		v2 = K2 * ES
		du.S = - v1
		du.ES = v1 - v2
		du.P = v2
		return du
	end
	
	# Apply QSSA on ES complex to reduce the model complexity
	function reduced_model!(du, u, p, t)
		@unpack ET, K1, KM1, K2 = p
		@unpack S = u
		du.S = -K2 * ET * hill(S, (KM1 + K2) / K1)
		return du
	end
	
	sol = solve(ODEProblem(full_model!, u0, tend, p))
	
	p1 = plot(sol)
	plot!(p1, sol, vars=((t, es) -> (t, p.ET - es), 0, 2), label="E")
	plot!(p1, xlabel="Time (arbitrary units)", ylabel="Concentration (arbitrary units)", legend=:right)
	
	u0R = LVector(S = sum(u0))
	solr = solve(ODEProblem(reduced_model!, u0R, tend, p))
	
	p2 = plot(sol, vars=(0, [1, 3]), line=(:dash), label=["S (full)", "P (full)"])
	plot!(p2, solr, line=(:blue), lab="S (reduced)")
	plot!(p2, solr, vars=((t, s)->(t, u0R[1] - s), 0, 1), line=(:red), lab="P (reduced)")
	plot!(p2, xlabel="Time (arbitrary units)",  ylabel="Concentration (arbitrary units)", xlims=(0.0,1.0), ylims=(0.0,5.0), legend = :right)
	
	return (p1, p2)
end

# ╔═╡ 6a30b6c0-867d-4e0b-bd0d-30f58944a11f
let
	figs = fig0303()
	plot(figs..., size=(800, 400))
end

# ╔═╡ 6a1d0241-e493-4319-ab7d-c4c17db18a82
function figp375()
	
	reducedmm(x, k) = x / k

	function model!(du, u, p, t)
		@unpack V0, VM1, VM2, VM3, KM1, KM2, KM3, _mm = p
		@unpack s1, s2, s3 = u
		v1 = VM1 * _mm(s1, KM1)
		v2 = VM2 * _mm(s2, KM2)
		v3 = VM3 * _mm(s3, KM3)
		du.s1 = V0 - v1
		du.s2 = v1 - v2
		du.s3 = v2 - v3
		return du
	end
	
	# first set of ICs
	u0 = LVector(s1=0.3, s2=0.2, s3=0.1)
	tend = 2.0
	
	p = (V0 = 2.0,
		VM1 = 9.0,
		VM2 = 12.0,
		VM3 = 15.0,
		KM1 = 1.0,
		KM2 = 0.4,
		KM3 = 3.0,
		_mm = hill)
	
	sol1 = solve(ODEProblem(model!, u0, tend, p))
	p1 = plot(sol1, ylims=(0.0, 0.8),
     title="Problem 3.7.5 (1)",
     xlabel="Time (arbitrary units)",
     ylabel="Concentration (arbitrary units)")
	
	# second set of ICs
	u02 = LVector(s1=6.0, s2=4.0, s3=4.0)
	
	sol2 = solve(ODEProblem(model!, u02, 4.0, p))
	p2 = plot(sol2, ylims=(0.0, 6.0),
     title="Problem 3.7.5 (2)",
     xlabel="Time (arbitrary units)",
     ylabel="Concentration (arbitrary units)")
	
	# Reduced model
	param2 = @set p._mm = reducedmm
	sol3 = solve(ODEProblem(model!, u0, tend, param2))
	p3 = plot(sol1, ylims=(0.0, 0.8),
     title="Problem 3.7.5 (1) (full vs reduced)",
     xlabel="Time (arbitrary units)",
     ylabel="Concentration (arbitrary units)",
     labels=["S1 " "S2 " "S3 "], ls=:dash)

	plot!(p3, sol3, labels=["S1 (reduced)" "S2 (reduced)" "S3 (reduced)"] )
	
	sol4 = solve(ODEProblem(model!, u02, 4.0, param2))
	
	p4 = plot(sol2, ylims=(0.0, 8.0),
     title="Problem 3.7.5 (2) (full vs reduced)",
     xlabel="Time (arbitrary units)",
     ylabel="Concentration (arbitrary units)",
     labels=["S1 " "S2 " "S3 "], ls=:dash)

	plot!(p4, sol4, labels=["S1 (reduced)" "S2 (reduced)" "S3 (reduced)"] )
	
	return (p1, p2, p3, p4)
end

# ╔═╡ 091c6dcf-343d-41ee-92b1-2585aecf12a6
let
	figs = figp375()
	plot(figs..., size=(800, 800))
end

# ╔═╡ 3b20843f-539d-455e-abcf-fe15e5fcf86a
function model401!(du, u, p, t)
    @unpack K1, K2, K3, K4, K5, N = p
    @unpack A, B = u
    
    v1 = K1 * hill(1, B, N)
    v5 = K5 * A
    
    du.A = v1 - v5 - K3 * A
    du.B = K2 + v5 - K4 * B
end

# ╔═╡ e97c3472-bb06-4968-92eb-336520e1c455
function figure0401()
	tend = 1.5
	params = (K1=20.0, K2=5.0, K3=5.0, K4=5.0, K5=2.0, N=4)
	
	u0s = (LVector(A=0.0, B=0.0), 
       LVector(A=0.5, B=0.6),
       LVector(A=0.17, B=1.1),
       LVector(A=0.25, B=1.9),
       LVector(A=1.85, B=1.70))
	
	sols = map(u0 -> solve(ODEProblem(model401!, u0, tend, params)), u0s)
	
	p1 = plot(sols[1], xlabel="Time", ylabel="Concentration", title="Fig. 4.2A (Time series)")
	
	p2 = plot(sols[1], vars=(1, 2), xlabel="[A]", ylabel="[B]", aspect_ratio=:equal,
     title="Fig. 4.2B (Phase plot)", ylims=(0.0, 2.0), xlims=(0.0, 2.0), 
     legend=nothing)
	
	p3 = plot()

	for sol in sols
		plot!(p3, sol, linealpha=0.5, legend = nothing)
	end

	plot!(p3, xlabel="Time", ylabel="Concentration", title="Fig. 4.3A (Time series)")
	
	p4 = plot()

	for sol in sols
		plot!(p4, sol, vars=(1, 2), linealpha=0.7, legend = nothing)
	end

	plot!(p4, aspect_ratio=:equal, title="Fig. 4.3B (Phase plot)", xlabel="[A]", ylabel="[B]", ylims=(0.0, 2.0), xlims=(0.0, 2.0), size=(600, 600))
	
	return (p1, p2, p3, p4)
	
end

# ╔═╡ e342c29c-24c0-4f7e-95b9-5e3cec57210b
let
	figs = figure0401()
	plot(figs..., size=(800, 800))
end

# ╔═╡ 1ead4f46-cad1-4916-bef6-23bed88059bd
function figure0404(; params = (K1=20.0, K2=5.0, K3=5.0, K4=5.0, K5=2.0, N=4),
	 				  tend = 1.5, 
	                  u0s = (LVector(A=0.0, B=0.0), 
							 LVector(A=0.5, B=0.6),
							 LVector(A=0.17, B=1.1),
							 LVector(A=0.25, B=1.9),
							 LVector(A=1.85, B=1.70)))
	
	# Functions for nullclines
	nullcline_a(b, p) = p.K1 / (p.K5 + p.K4)  * hill(1, b, p.N)
	nullcline_a(b) = nullcline_a(b, params)
	nullcline_b(b, p) = (p.K4*b - p.K2) / p.K5
	nullcline_b(b) = nullcline_b(b, params)
	
	# Tweaking arrow length
	scale_du(du, scale=20) = du ./ (norm(du)^0.5 * scale)
	
	# vector field
	function df(x, y)
		u = LVector(A=x, B=y)
		du = similar(u)
		model401!(du, u, params, 0.0)
		return scale_du(du)
	end
	
	sols = map(u0 -> solve(ODEProblem(model401!, u0, tend, params)), u0s)
	
	xx = [x for y in 0.0:0.1:2.0, x in 0.0:0.1:2.0]
	yy = [y for y in 0.0:0.1:2.0, x in 0.0:0.1:2.0]
	p1 = quiver(xx, yy, quiver=df, line=(:lightblue))

	for sol in sols
		plot!(p1, sol, vars=(1, 2), legend = nothing, line=(:blue))
	end
	
	plot!(p1, aspect_ratio=:equal, title="Fig. 4.4A (Phase Plot with vector field)", 
      xlabel="[A]", ylabel="[B]", xlim=(0.0, 2.0), ylim=(0.0, 2.0), size=(600, 600))
	
	# Figure 4.5A
	p2 = plot(aspect_ratio=:equal, title="Fig. 4.5A, Phase plot with nullclines")

	# Phase plots
	for sol in sols
		plot!(p2, sol, vars=(1, 2), linealpha=0.7, lab=nothing)
	end

	# Parametric plotting for nullcline
	plot!(p2, nullcline_a, identity, 0.0, 2.0, label="A nullcline", line=(:black, :dot))
	plot!(p2, nullcline_b, identity, 0.0, 2.0, label="B nullcline", line=(:black, :dash))
	plot!(p2, xlim=(0.0, 2.0), ylim=(0.0, 2.0), legend=:bottomright, size=(600, 600), xlabel="[A]", ylabel="[B]")

	
	p3 = quiver(xx, yy, quiver=df, line=(:lightblue), title="Fig. 4.5B, Vector field with nullclines", xlabel="[A]", ylabel="[B]")
	plot!(p3, nullcline_a, identity, 0.0, 2.0, label="A nullcline", line=(:black, :dot))
	plot!(p3, nullcline_b, identity, 0.0, 2.0, label="B nullcline", line=(:black, :dash))
	plot!(p3, aspect_ratio=1.0, xlim=(0.0, 2.0), ylim=(0.0, 2.0), legend=:bottomright, size=(600, 600))

	return (p1, p2, p3)
end

# ╔═╡ 41b0b13d-045e-4d14-aaa3-e676fa1a23e5
fig0404, fig0405a, fig0405b = figure0404()

# ╔═╡ f1453c64-f07f-4a6c-bac6-8b56373ac434
fig0404

# ╔═╡ c6a20d1f-e3c4-4dc4-802c-6b06d0968f68
fig0405a

# ╔═╡ 563c736a-9b2f-42ec-9d7e-e2868122f7eb
fig0405b

# ╔═╡ c662c3a3-e83c-4f78-897b-693e76de8de8
"Model of symmetric network from Figure 4.6. This code generates Figures 4.7, 4.8, 4.9, and 4.19A"
function model0407!(du, u, p, t)
    @unpack K1, K2, K3, K4, N1, N2 = p
    @unpack s1, s2 = u
    du.s1 = K1 * hill(1, s2, N2) - K3 * s1
    du.s2 = K2 * hill(1, s1, N1) - K4 * s2
    return du
end

# ╔═╡ af0f7d8a-22e7-4224-a56e-17fc48a385de
let
	tend = 4.0
	p1 = (K1=20.0, K2=20.0, K3=5.0, K4=5.0, N1=1.0, N2=4.0)

	u0s = (LVector(s1=3.0, s2=1.0), 
		   LVector(s1=1.0, s2=3.0))

	sols = map(u0 -> solve(ODEProblem(model0407!, u0, tend, p1)), u0s)
	
	p1 = plot(sols[1], xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.7A (1)")
    p2 = plot(sols[2], xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.7A (2)")
    fig47a = plot(p1, p2, layout=(2, 1), size=(600, 600))
end

# ╔═╡ 6ebffaf6-cf39-43a2-a29b-bdd9498bb375
let
	tend = 4.0
	p1 = (K1=20.0, K2=20.0, K3=5.0, K4=5.0, N1=4.0, N2=4.0)

	u0s = (LVector(s1=3.0, s2=1.0), 
		   LVector(s1=1.0, s2=3.0))

	sols = map(u0 -> solve(ODEProblem(model0407!, u0, tend, p1)), u0s)
	
	p1 = plot(sols[1], xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.8A (1)")
    p2 = plot(sols[2], xlabel="Time", ylabel="Concentration", legend=:right, title= "Fig 4.8A (2)")
    fig47a = plot(p1, p2, layout=(2, 1), size=(600, 600))
end

# ╔═╡ f53cbbe1-2366-4c38-a84f-dc75c9f48e3d
let
	p1 = (K1=20.0, K2=20.0, K3=5.0, K4=5.0, N1=1.0, N2=4.0)
	
	function df(x, y)
		u = LVector(s1=x, s2=y)
		du = similar(u)
		model0407!(du, u, p1, 0.0)

		# Tweaking arrow length
		du ./ (norm(du)^0.5 * 20)
	end
	nullclineS1(s2, p) = p.K1 / p.K3 * hill(1, s2, p.N2)
	nullclineS1(s2) = nullclineS1(s2, p1)
	nullclineS2(s1, p) = p.K2 / p.K4 * hill(1, s1, p.N1)
	nullclineS2(s1) = nullclineS2(s1, p1)
	
	r = LinRange(0.0, 5.0, 20)
	xx = [x for y in r, x in r]
	yy = [y for y in r, x in r]
	pl = quiver(xx, yy, quiver=df, line=(:lightblue))

	plot!(pl, nullclineS1, identity, 0.0, 5.0, lab="Nullcline S1", line=(:dash, :red))
	plot!(pl, identity, nullclineS2, 0.0, 5.0, lab="Nullcline S2", line=(:dash, :blue))
	plot!(pl, title="Fig 4.7 B", xlim=(0.0, 5.0), ylim=(0.0, 5.0), aspect_ratio = 1.0, size = (600, 600))
	
	pl
end

# ╔═╡ 51332af1-7e4a-4f4a-b68d-7380344ad010
let
	p1 = (K1=20.0, K2=20.0, K3=5.0, K4=5.0, N1=4.0, N2=4.0)
	
	function df(x, y)
		u = LVector(s1=x, s2=y)
		du = similar(u)
		model0407!(du, u, p1, 0.0)

		# Tweaking arrow length
		du ./ (norm(du)^0.5 * 20)
	end
	nullclineS1(s2, p) = p.K1 / p.K3 * hill(1, s2, p.N2)
	nullclineS1(s2) = nullclineS1(s2, p1)
	nullclineS2(s1, p) = p.K2 / p.K4 * hill(1, s1, p.N1)
	nullclineS2(s1) = nullclineS2(s1, p1)
	
	r = LinRange(0.0, 5.0, 20)
	xx = [x for y in r, x in r]
	yy = [y for y in r, x in r]
	
	pl = quiver(xx, yy, quiver=df, line=(:lightblue))

	plot!(pl, nullclineS1, identity, 0.0, 5.0, lab="Nullcline S1", line=(:dash, :red))
	plot!(pl, identity, nullclineS2, 0.0, 5.0, lab="Nullcline S2", line=(:dash, :blue))
	plot!(pl, title="Fig 4.8 B", xlim=(0.0, 5.0), ylim=(0.0, 5.0), aspect_ratio = 1.0, size = (600, 600))
	
	r2 = LinRange(1.0, 1.5, 20)
	xx2 = [x for y in r2, x in r2]
	yy2 = [y for y in r2, x in r2]

	pl2 = quiver(xx2, yy2, quiver=(x, y)->df(x,y) ./ 5, line=(:lightblue))

	plot!(pl2, nullclineS1, identity, r2[1], r2[end], lab="Nullcline S1", line=(:dash, :red))
	plot!(pl2, identity, nullclineS2, r2[1], r2[end], lab="Nullcline S2", line=(:dash, :blue))
	plot!(pl2, title="Fig 4.8 B (close up)", xlim=(r2[1], r2[end]), ylim=(r2[1], r2[end]), 
		  aspect_ratio = 1.0, size = (600, 600), xlabel="[S1]", ylabel="[S2]")
	
	plot(pl, pl2, size=(1000, 500))
end

# ╔═╡ b9444e52-5b7a-41da-9885-62af0fc350e7
let
	nullclineS1(s2, p) = p.K1 / p.K3 * hill(1, s2, p.N2)
	nullclineS1(s2) = nullclineS1(s2, p1)
	nullclineS2(s1, p) = p.K2 / p.K4 * hill(1, s1, p.N1)
	nullclineS2(s1) = nullclineS2(s1, p1)
	
	# Fig 4.19 A
	pls = map((8.0, 16.0, 20.0, 35.0)) do k1
		p = (K1=k1, K2=20.0, K3=5.0, K4=5.0, N1=4.0, N2=4.0)
		plot(s2 -> nullclineS1(s2, p), s2 -> s2, 0.0, 7.0, lab="Nullcline S1")
		plot!(s1 -> s1, s1 -> nullclineS2(s1, p), 0.0, 7.0, lab="Nullcline S2")
		plot!(title = "K1 = $k1", xlim=(0.0, 7.0), ylim=(0.0, 7.0), 
		  aspect_ratio = 1.0, size = (800, 800), xlabel="[S1]", ylabel="[S2]")
	end

	plot(pls...)
end

# ╔═╡ e3bb6d7f-9f17-4bd2-a40e-783114404b96
function figure0415(; param = (K0 = 8.0, K1 = 1.0, K2 = 5.0, N = 2),
	                  r = LinRange(0.0, 4.0, 20),
	                  tend = 8.0,
	                  figtitle="Fig 4.15")

	u0s = ( LVector(s1=1.5, s2=1.0), LVector(s1=0.0, s2=1.0),
        	LVector(s1=0.0, s2=3.0), LVector(s1=2.0, s2=0.0))
	sols = map(u0 -> solve(ODEProblem(model0415!, u0, tend, param)), u0s)

	# Fig 4.15 A
	p1 = plot(sols[1], xlabel="Time", ylabel="Concentration", title ="$figtitle (A)", xlims=(0.0, 8.0))
	
	# Fig 4.15 B
	"Vetor field"
	function ∂F(x, y)
		u = LVector(s1=x, s2=y)
		dxdy = model0415!(similar(u), u, param, 0.0)
		return dxdy ./ (norm(dxdy)^0.5 * 20)
	end
	
	nullcline_s1(s2, p=param) = (p.K0 / p.K1) * hill(1, s2, p.N)
	nullcline_s2(s2, p=param) = (p.K2 * s2) / (p.K1 * (1 + s2^p.N))
	
	
	xx = [x for y in r, x in r]
	yy = [y for y in r, x in r]
	p2 = quiver(xx, yy, quiver=∂F, line=(:lightblue))

	for sol in sols
		plot!(p2, sol, vars=(1, 2), label=nothing)
	end
	
	rMin, rMax = r[1], r[end]
	
	plot!(p2, nullcline_s1, identity, rMin, rMax, label="Nullcline S1", line=(:dash, :red))
	plot!(p2, nullcline_s2, identity, rMin, rMax, label="Nullcline S2", line=(:dash, :blue))
	plot!(p2, title = "$figtitle (B)", xlabel="[S1]", ylabel="[S2]", 
      xlims=(rMin, rMax), ylims=(rMin, rMax), aspect_ratio=:equal, size=(700, 700))
	
	return (p1, p2)
end

# ╔═╡ a0772af3-3aa3-42d3-b017-7ff99780ecfe
fig415a, fig415b = figure0415()

# ╔═╡ 97ed8540-3804-41f5-9178-d7c65b85e75c
fig415a

# ╔═╡ 167c703e-4948-487e-94c5-701f3cc15ebe
fig415b

# ╔═╡ 64cb5664-bf93-4828-b2eb-d689c548d147
fig416a, fig416b = figure0415(param = (K0 = 8.0, K1 = 1.0, K2 = 5.0, N = 2.5), tend = 1000.0, figtitle="Fig 4.16")

# ╔═╡ 03c7c066-246a-406a-9c4c-1aa83f3bfd90
fig416a

# ╔═╡ 6b60cf0c-097e-450a-8135-c77e055835c1
fig416b

# ╔═╡ 9743e59f-14c7-4336-99cc-5a2e0381e408
let
	param = (K0 = 8.0, K1 = 1.0, K2 = 5.0, N = 2.5)
	sol417 = solve(ODEProblem(model0415!, LVector(s1=2.0, s2=1.5), 10.0, param))
	r = LinRange(0.0, 4.0, 20)
	xx = [x for y in r, x in r]
	yy = [y for y in r, x in r]
	
	"Vetor field"
	function ∂F(x, y)
		u = LVector(s1=x, s2=y)
		dxdy = model0415!(similar(u), u, param, 0.0)
		return dxdy ./ (norm(dxdy)^0.5 * 20)
	end
	
	nullcline_s1(s2, p=param) = (p.K0 / p.K1) * hill(1, s2, p.N)
	nullcline_s2(s2, p=param) = (p.K2 * s2) / (p.K1 * (1 + s2^p.N))
	
	
	xx = [x for y in r, x in r]
	yy = [y for y in r, x in r]
	p2 = quiver(xx, yy, quiver=∂F, line=(:lightblue))
	
	
	quiver(xx, yy, quiver=∂F, line=(:lightblue))
	plot!(sol417, vars=(1, 2), label=nothing, line=(:black), arrow=0.4)


	plot!(nullcline_s1, identity, 0.0, 4.0, label="Nullcline S1", line=(:dash, :red))
	plot!(nullcline_s2, identity, 0.0, 4.0, label="Nullcline S2", line=(:dash, :blue))
	plot!(title = "Fig. 4.17", xlabel="[S1]", ylabel="[S2]", 
		  xlims=(1.0, 3.0), ylims=(1.0, 3.0), aspect_ratio=:equal, size=(700, 700))
end

# ╔═╡ 23451b22-a3c0-453f-94bb-bc2b6a97d19b
function model0418!(du, u, p, t)
    @unpack K1, K2, K3, K4, K5, N = p
    @unpack A, B = u
    
    v1 = K1 * hill(1, B, N)
    v5 = K5 * A
    
    du.A = v1 - v5 - K3 * A
    du.B = K2 + v5 - K4 * B
end

# ╔═╡ 0018fe66-ded6-4033-8d84-6e5196051316
let
	p = (K1=20.0, K2=5.0, K3=5.0, K4=5.0, K5=2.0, N=4)
	k1Range = LinRange(0.0, 1000.0, 100)
	u0 = LVector(A=0.0, B=0.0)
	
	# See also Ensemble analysis: https://diffeq.sciml.ai/stable/features/ensemble/
	a = map(k1Range) do k1
		p1 = @set p.K1 = k1
		prob = SteadyStateProblem(model0418!, u0, p1)
		sol = solve(prob)
		sol[1]
	end
	
	plot(k1Range, a, title = "Fig 4.18 Continuation diagram", xlabel = "K1" , ylabel= "Steady state [A]", leg=nothing, ylim=(0.0, 4.0))
end

# ╔═╡ Cell order:
# ╠═5ba678b7-7251-4003-9f6e-df6339335635
# ╠═0a52a098-9f37-4d89-96d1-813a0abc2f3b
# ╠═33edf7e8-d3d6-4ced-844b-9fd2c0665293
# ╠═8fbe196a-0cb1-408f-bb89-d5a2bcdb3a40
# ╠═3a7e2fc4-05ca-4e80-b2b2-48a0513dfc4d
# ╠═f0fc6b02-c070-49d8-bce3-777758c6d482
# ╠═eb14ab5a-e7d7-415a-8cbc-7d5152c4e01b
# ╠═333892b7-e248-422f-b16e-0cde638dd9ee
# ╠═a4e81834-0827-4d0a-a2c8-6fff6002ddc9
# ╠═478ee763-2b39-4716-86fb-ad49194b5125
# ╠═68786419-44ab-48b9-b46c-135af4ffaf95
# ╠═805ba136-caea-4391-9e9e-b3275341613e
# ╠═76aab80e-40eb-473c-8f82-49c67c304086
# ╠═7623e4ca-d9e6-4e97-9dee-95efe616f959
# ╠═dcff0bfd-6f4b-4ff8-8532-59888754bb52
# ╠═037facee-27ca-4a42-a79c-aa8cad6a1187
# ╠═1a4984d2-1e8c-4f3b-94c4-fb937e9de6f1
# ╠═5021a671-13d3-428e-bc92-4e3cf930fa0b
# ╠═0ebae16b-3b6c-4298-9812-fac0227e8ed7
# ╠═ca209a4f-ca44-4d0c-8c2a-83954ccd964a
# ╠═da0b9be3-36d8-4f9b-b03d-30270111e862
# ╠═58d5e228-be35-4b1f-b9f2-368e9c99ee20
# ╠═6a30b6c0-867d-4e0b-bd0d-30f58944a11f
# ╠═91b7d696-109a-44b1-9631-69e514a290ac
# ╠═6a1d0241-e493-4319-ab7d-c4c17db18a82
# ╠═091c6dcf-343d-41ee-92b1-2585aecf12a6
# ╠═a4a64652-a4b9-4345-9652-252daff344d8
# ╠═3b20843f-539d-455e-abcf-fe15e5fcf86a
# ╠═e97c3472-bb06-4968-92eb-336520e1c455
# ╠═e342c29c-24c0-4f7e-95b9-5e3cec57210b
# ╠═f3c1efd0-b76d-496a-83bd-f2c8c57d2b94
# ╠═1ead4f46-cad1-4916-bef6-23bed88059bd
# ╠═41b0b13d-045e-4d14-aaa3-e676fa1a23e5
# ╠═f1453c64-f07f-4a6c-bac6-8b56373ac434
# ╠═c6a20d1f-e3c4-4dc4-802c-6b06d0968f68
# ╠═563c736a-9b2f-42ec-9d7e-e2868122f7eb
# ╠═04fc83f3-1d1b-4ba7-a08c-7fc646a73de3
# ╠═c662c3a3-e83c-4f78-897b-693e76de8de8
# ╠═af0f7d8a-22e7-4224-a56e-17fc48a385de
# ╠═f53cbbe1-2366-4c38-a84f-dc75c9f48e3d
# ╠═6ebffaf6-cf39-43a2-a29b-bdd9498bb375
# ╠═51332af1-7e4a-4f4a-b68d-7380344ad010
# ╠═b9444e52-5b7a-41da-9885-62af0fc350e7
# ╠═9519b1ba-726f-4a25-8993-fbe897e998cd
# ╠═96e7a435-3411-4412-89d9-fb3b42e62a7c
# ╠═7b747fe2-05bc-40b5-9fcf-f6e11125cc91
# ╠═ec439af8-b24b-4167-8b79-756955e1375d
# ╠═e3bb6d7f-9f17-4bd2-a40e-783114404b96
# ╠═a0772af3-3aa3-42d3-b017-7ff99780ecfe
# ╠═97ed8540-3804-41f5-9178-d7c65b85e75c
# ╠═167c703e-4948-487e-94c5-701f3cc15ebe
# ╠═64cb5664-bf93-4828-b2eb-d689c548d147
# ╠═03c7c066-246a-406a-9c4c-1aa83f3bfd90
# ╠═6b60cf0c-097e-450a-8135-c77e055835c1
# ╠═9743e59f-14c7-4336-99cc-5a2e0381e408
# ╠═998ba19b-63ca-4079-ab4b-f6de5d0d61b0
# ╠═23451b22-a3c0-453f-94bb-bc2b6a97d19b
# ╠═0018fe66-ded6-4033-8d84-6e5196051316
# ╠═32825174-50e3-42b9-9c5d-ea744b36fefb
# ╠═c94276d8-e844-437d-8ede-742bba44a2fd
# ╠═6ede64c4-4851-4506-88d9-7e58d49c35e9
# ╠═287b9e43-2313-4724-b873-d9bd87a8948d
