### A Pluto.jl notebook ###
# v0.14.6

using Markdown
using InteractiveUtils

# ╔═╡ 36d178f5-95fd-455f-954e-c85a5c92eebc
begin
    import Pkg
    Pkg.activate(mktempdir())
    Pkg.add([
        Pkg.PackageSpec(name="PlutoUI", version="0.7"),
        Pkg.PackageSpec(name="Plots", version="1"),
    ])
	using PlutoUI
	using Plots
	using Random
	PlutoUI.TableOfContents()
end

# ╔═╡ 8a1e5a8b-7fbd-4818-af91-a52df41e0e51
md"""
# Line plots

## Line plots (I)
"""

# ╔═╡ f27bb3ad-ba52-4334-9aac-a9fa48b07f47
x = 0:0.1:2pi

# ╔═╡ 7a692a81-a0fe-4d67-b6b1-697d976f89af
y1 = cos.(x)

# ╔═╡ 03e58a6b-02b0-4133-91f4-775a6aee6a39
y2 = sin.(x)

# ╔═╡ 9b2ea567-bdec-4211-a5d2-670a50b82c7f
begin
	p1 = plot(x, y1, c="blue", linewidth=3)
	plot!(p1, x, y2, c="red", line=:dash)
	title!(p1, "Trigonometric functions")
	xlabel!(p1, "angle")
	ylabel!(p1, "sin(x) and cos(x)")
	plot!(p1, xlims=(0,2pi), ylims=(-2, 2), size=(600, 600))
	p1
end

# ╔═╡ a8aa4f51-3e82-47d2-aa09-89ece6c7bbd4
begin
	p2 = plot(x, y1,
    		  c="blue",
    		  linewidth=3,
    		  title="Trigonometric functions",
    		  xlabel="angle",
    		  ylabel="sin(x) and cos(x)")
	plot!(p2, x, y2, c="red", line=:dash)

	plot!(p2, xlims=(0,2pi), ylims=(-2, 2), size=(600, 600))
end

# ╔═╡ e649e45c-7e3c-4b59-b0a4-7ca650077fac
md"""
## Line plots (II)
"""

# ╔═╡ 15d78933-4ac8-4c60-b080-e8d32a5022d2
begin
	Random.seed!(2018)
	time = 30
	walker1 = cumsum(randn(time))
	walker2 = cumsum(randn(time))
	walker3 = cumsum(randn(time))
	walker4 = cumsum(randn(time))
	walker5 = cumsum(randn(time))
	plot(1:time, [walker1 walker2 walker3 walker4 walker5],
    xlabel="time",
    ylabel="position",
    label=["walker1" "walker2" "walker3" "walker4" "walker5"],
    leg=:bottomleft)
end

# ╔═╡ be3813ec-be35-11eb-0f75-956c2874118b
md"""
# Appendix
Running environment and some auxillary functions.
	
**You can use [this helper tool](https://fonsp.com/article-test-3/pkghelper.html) to generate these commands!**
"""

# ╔═╡ Cell order:
# ╠═8a1e5a8b-7fbd-4818-af91-a52df41e0e51
# ╠═f27bb3ad-ba52-4334-9aac-a9fa48b07f47
# ╠═7a692a81-a0fe-4d67-b6b1-697d976f89af
# ╠═03e58a6b-02b0-4133-91f4-775a6aee6a39
# ╠═9b2ea567-bdec-4211-a5d2-670a50b82c7f
# ╠═a8aa4f51-3e82-47d2-aa09-89ece6c7bbd4
# ╠═e649e45c-7e3c-4b59-b0a4-7ca650077fac
# ╠═15d78933-4ac8-4c60-b080-e8d32a5022d2
# ╠═be3813ec-be35-11eb-0f75-956c2874118b
# ╠═36d178f5-95fd-455f-954e-c85a5c92eebc
