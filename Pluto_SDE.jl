### A Pluto.jl notebook ###
# v0.12.20

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 5b4cd610-68df-11eb-148b-8f94bd42f945
begin 
	rootdir = @__DIR__; 
	# @__DIR__ is the directory of current .jl file
	srcdir = joinpath(rootdir, "src");
end

# ╔═╡ a3a0bfee-e845-11ea-004f-df9b39aade97
begin
	using Pkg
	Pkg.activate(rootdir)
	Pkg.instantiate();
	using FFTW
	using Plots
	using PlutoUI
	using DifferentialEquations
	using Formatting
	using LaTeXStrings
	using Printf
end

# ╔═╡ 366389d0-68e3-11eb-19dd-d78bd237ab71
PlutoUI.TableOfContents()

# ╔═╡ 7acb7db0-06ac-11eb-3218-b1d3fd2d72c7
md"### Particular solution using Euler method (Cyganowski, 2001)"

# ╔═╡ 054caa20-795a-11eb-2eba-ddc9b516b4dc
@bind frictiontype Radio(["Coulomb","viscous"], default = "Coulomb")

# ╔═╡ f8fc0520-699d-11eb-0d45-5932beaa149f
md"
duration $(@bind totalTimePower Slider(0.1:0.5:3.5, default=2));
diffusion $(@bind D Slider(0.05:0.05:0.2, show_value = true,default=2));

friction $(@bind F_C Slider(0.05:0.05:0.2, show_value = true,default=2));
driving force $(@bind F_ext_to_F_C_ratio Slider(0:0.05:0.95, show_value = true,default=0.5))

datapoints $(@bind NwPower Slider(1:4, default=3,show_value = true))
"

# ╔═╡ 6bfa7a80-8018-11eb-2c0f-3b411d6291d4
md"#### variables"

# ╔═╡ f6cc0882-716e-11eb-3380-95579d37d720
md"
### Numerically solving the Stochastic Differential Equation (SDE) using `DifferentialEquations.jl`
"

# ╔═╡ 464ce522-717c-11eb-2810-7b244052747d
md"
For 

$\mathrm{d}v = f(v) \mathrm{d}t + g(v) \mathrm{d}W$,
we have the following solution if $f(v) = \alpha v$ and $g(v) = \beta v$

$v(t, W_t) = v_0 \exp{((\alpha-\frac{\beta^2}{2})t+\beta W_t)}$

> see the [official example](https://diffeq.sciml.ai/stable/tutorials/sde_example/)

"

# ╔═╡ 06cd3970-801a-11eb-1065-49fad5c8ba21
md"Run this session: $(@bind RunIt2 CheckBox(default=false))"

# ╔═╡ 3b7d82e0-8012-11eb-2bea-fbcfebe6ffe4
md"
### Analytical solution
"

# ╔═╡ 5e9c7290-8012-11eb-3191-79e66024157d
md"
#### from Fokker-Planck equation (FPE)
"

# ╔═╡ 6f1de950-8012-11eb-399a-432019e28e00
md"
#### according to Ito calculus:

For the Langevin equation

$\frac{\mathrm{d}v(t)}{\mathrm{d}t} = -{\gamma}v(t) + \Gamma(t)$

its solution is 

$v(t) = \mathrm{e}^{-t/\tau_B}v(0) + \int_0^t \mathrm{e}^{-(t-s)/\tau_B}\sqrt{2D}\mathrm{d} W(s)$ 

where $\tau_B \sim 1/\gamma$ is the relaxation of the particle velocity.
> see Eq. 6.3, 6.5, 6.9 of [this](http://physics.gu.se/~frtbm/joomla/media/mydocs/LennartSjogren/kap6.pdf).
"

# ╔═╡ e0acaee0-84c1-11eb-245b-4b0042d2642c
md"### Spectral Analysis"

# ╔═╡ f300c2c0-84c1-11eb-335f-816c1c15fe2a
md"#### boxcar-shaped velocity history"

# ╔═╡ d3de7090-84cb-11eb-1e8c-635d7db04e99
@bind numcars Slider(1:5, show_value = true, default = 3)

# ╔═╡ 1a829760-84c2-11eb-10d9-adc0e83b8a45
let
	include(joinpath(srcdir,"randboxcars.jl"));
	x = collect(range(0,10,length = 5000));
	ys = fill(0, (length(x),numcars))
	ys = randboxcars(x, ys);
	p_array = []; # or p_arr1 = Plots.Plot{Plots.GRBackend}[] to create an empty array of plot
	for i = 1:numcars
		p = plot(x,ys[:,i]);
		push!(p_array, p);
	end
	p = plot(x, sum(ys,dims = 2),title="superposition of boxcars above");
	push!(p_array, p);
	subplotheights = fill(1,numcars+1);
	subplotheights[end] = numcars;
	subplotheights = subplotheights./sum(subplotheights); # normalize subplot heights
	lout = grid(numcars+1,1, heights = subplotheights);
	plot(p_array..., layout = lout, legend = false)
end

# ╔═╡ 80d83300-897a-11eb-369c-e7a5bfa489fa
@bind numcars2 Slider(1:500, show_value = true, default = 50)

# ╔═╡ 4e1ee580-84ca-11eb-2693-49048d8f1dd3
let
	include(joinpath(srcdir,"randboxcars.jl"));
	x = collect(range(0,10,length = 5000));
	ys = randboxcars(x, numcars2; boxwidthstd=0.1);
	p1 = plot(x, ys, title="$numcars2 boxcars");
	p2 = plot(abs.(fft(ys)), xaxis=:log, yaxis=:log, title="amplitude spectrum");
	plot(p1,p2,layout = (1,2),size=(600,250), legends= false);
end

# ╔═╡ 4ec09530-8013-11eb-2870-dd7421db7daa
md"### Common variables"

# ╔═╡ 703d4ec0-f728-11ea-150b-e9e101b7acbd
begin
	totalTime = Float64(10^totalTimePower);
	Nw = Int64(floor(10^NwPower));
	Y0 = 0.0;
	F_ext = Float64(F_C*F_ext_to_F_C_ratio);
	titlefSz = 10;
end

# ╔═╡ 63564260-8018-11eb-058f-292b305207bc
begin
	dt = Float64(totalTime/(Nw-1));
	traceT = collect(range(1,step = dt,length = Nw));
end

# ╔═╡ 72ddc3b0-79c0-11eb-0b19-19eb81a48e41
let
	include(joinpath(srcdir,"SDE.jl"))
	include(joinpath(srcdir,"circles.jl"));
	drift_x(v) = -F_C*v+F_ext;
	drift_y(v) = -F_C*v;
	traceY = SDE(traceT, [drift_x, drift_y], D, [Y0, Y0]);
	absY = traceY[:,1].^2 .+ traceY[:,2].^2;
	p1 = plot(traceT, traceY[:,1], xlabel = "", ylabel="v_x(t)", 
		titlefontsize = titlefSz,
		title="sample path (1D)\n in velocity space",legend = false);
	p2 = plot(traceT, traceY[:,2], xlabel = "t",ylabel="v_y(t)",legend = false);
	p3 = plot(traceY[:,1], traceY[:,2], seriestype=:path, # default: seriestype=:line
		xlabel="v_x", ylabel="v_y", titlefontsize = titlefSz,
		title="sample path (2D)\n in velocity space",legend = false); 
	d = [0.0 0.0];
	dts = fill(NaN,size(traceY))
	for i = 1:length(traceT)
		yi = traceY[i:i,:];
		d = d + yi*dt; 
		dts[i,:] = d;
	end
	p4 =plot(circles.(dts[:,1],dts[:,2],absY.*3), 
		seriestype=:shape, c=:red, fillalpha = 0.05, lw = 0,
		aspectratio = 1, titlefontsize = titlefSz, 
		title = "sample path (2D)\n in space");
## since there are too many circles, if specifying labels in the following way will be extremely time consuming.
	# p4legend = fill("", length(absY)+1);
	# p4legend[end-1:end] =  ["local slip magitude" "sample path in space"];
	plot!(p4, dts[:,1],dts[:,2],seriestype=:path,xlabel="x",ylabel = "y",
			legend = false);
		  # label = p4legend)
	lout1 = @layout [
		            [p1
			         p2] p3 p4
		           ];
	# lout2 = @layout [
	# 	            grid(2,1) grid(1,1) grid(1,1)
	# 	           ]; # both lout1 and lout2 gives almost the same result!
	plot(p1,p2,p3,p4,layout = lout1, size=(690,250));
end

# ╔═╡ 2cc390b2-69a0-11eb-2e8d-273d9dd6bee6
# noted that the space after "$" is very important!
begin
strT = @sprintf("%.2f",totalTime);
strdt = @sprintf("%.2E",dt);
md"
$ 
T = $strT;
D = $D;
F\_C = $(F_C); 
F\_{ext} = $(F_ext);
\text{datapoints} = $Nw; 
\mathrm{d}t = $strdt
$
"
end

# ╔═╡ a50f5ac0-8009-11eb-2235-5d7e75f62a1f
md"
With the drift term $f(v)=-F_C|v|+F_\text{ext}$ and the diffusion $g(v) = \sqrt{2D}$, by  using `DifferentialEquations.jl` we can get one sample path (*left*), ensemble averaged sample paths (*middle*), and distribution of ensemble solutions at $t=$ $(strT) (*right*) as:
"

# ╔═╡ cf170040-e845-11ea-145d-f716ad472b84
begin
if frictiontype=="viscous"
	drift(v) = -F_C*v+F_ext;
elseif frictiontype=="Coulomb"
	drift(v) = -F_C*sign(v)+F_ext;
	P0(fric,Fdrive,D) = (Fdrive^2- fric^2)/(-2*fric*D); 
	# Pst(fric,Fdrive,D,v) = P0(fric,Fdrive,D)*exp.(Fdrive*v/D).*exp.(-fric*v.*sign.(v)/D);
	P00 = P0(F_C,F_ext,D);
	Pst(fric,Fdrive,D,v) = P00*exp.(Fdrive*v/D).*exp.(-fric*v.*sign.(v)/D);
	# positive_expectation() = 
end
end

# ╔═╡ fe09c7b0-e846-11ea-3458-3d72f807da89
let
	include(joinpath(srcdir,"SDE.jl"));
	traceY = SDE(traceT,drift,Float64(D),Y0);
	predictat = totalTime;
	ensembleY = SDE(dt,predictat,10000,drift,D,Y0); 
	v = range(minimum(traceY),maximum(traceY),length=5000);
	p1 = plot(traceT,traceY,xlabel = "t",ylabel= "v(t)",legend= false);
	p2 = histogram(traceY,bins = 50,normalize=:pdf,
		xlabel = "v",ylabel= "P(t)", labels = "");
	plot!(p2,v,Pst(F_C,F_ext,D,v),labels = "FPE sol.",legend=true);
	# :pdf, :density, :probability or :none
	p3 = histogram(ensembleY, bins=50, normalize =:pdf, xlabel = "v", 
		ylabel = "P(t)", labels = "");
	v2 = range(minimum(ensembleY),maximum(ensembleY),length=5000);
	plot!(p3,v2,Pst(F_C,F_ext,D,v2),labels = "FPE sol.");
	plot(p1,p2,p3,layout = (1,3),
	title = ["one sample path" "time averaged\n distribution" "ensemble averaged\n distribution"],titlefontsize=titlefSz, size = (690, 190))
end

# ╔═╡ 0d59e7b0-7e63-11eb-13d2-178ea4c4d09e
let 
	# defining the problem
	if RunIt2
	u₀=0
	dt2 = totalTime/40;
	alg = SRIW1();# algorithm: EM(), SRIW1(), Tsit5(), etc...
	f(u,p,t) = -F_C*sign(u)+F_ext;
	g(u,p,t) = sqrt(2*D);
	tspan = (0.0,totalTime)
	prob = SDEProblem(f,g,u₀,tspan)
	sol = solve(prob,alg,dt=dt2)
	
	# single sample path
	v = range(minimum(sol.u),maximum(sol.u),length=5000);
	p1 = plot(sol.t,sol.u,xlabel = "t",ylabel= "v(t)", 
		title = "one sample path",legend = false);
	
	# ensemble prediction
	using DifferentialEquations.EnsembleAnalysis
	ensembleprob = EnsembleProblem(prob);
	solE = solve(ensembleprob,alg,EnsembleThreads();
		trajectories=1000,dt=dt2) # EM(),EnsembleThreads() is optional
	
	ensum1 = EnsembleSummary(solE);
	ensum2 = EnsembleSummary(solE, quantiles = [0.25, 0.75]);
	p2 = plot(ensum1,labels = "middle 95%");
	plot!(p2, ensum2, title = "Ensemble summary", labels = "middle 50%",legend=true)
	
	ensembleu = map(x-> x.u[end],solE);
	p3 = histogram(ensembleu, bins = 50, normalize = :pdf, 
		xlabel = "v", ylabel = "P(v)",labels = "");
	v = range(minimum(ensembleu),maximum(ensembleu),length=5000);
	plot!(p3, v,Pst(F_C,F_ext,D,v), labels = "FPE sol.")
	plot(p1,p2,p3,layout = (1,3), size= (690,190), titlefontsize=titlefSz)
	end

end

# ╔═╡ 6b13e4b0-8010-11eb-0826-f74b42b9c37a
md"
# CHECKPOINT: 

- check how to use multithreading [ref](https://diffeq.sciml.ai/stable/features/ensemble/#ensemble)

- read julia debugger (learn how to use @run and @enter) [here it is](https://julialang.org/blog/2019/03/debuggers/)

- see if there is a better way to pass the function `drift` or an array of drift functions (e.g. `[drift1, drift2]`) with assertion into the function `SDE.jl`.

- 
	- read it: https://stackoverflow.com/questions/63618640/type-of-function-in-julia
	- https://stackoverflow.com/questions/46028069/how-do-i-specify-a-functions-type-signature-in-julia

"

# ╔═╡ Cell order:
# ╠═366389d0-68e3-11eb-19dd-d78bd237ab71
# ╠═5b4cd610-68df-11eb-148b-8f94bd42f945
# ╠═a3a0bfee-e845-11ea-004f-df9b39aade97
# ╟─7acb7db0-06ac-11eb-3218-b1d3fd2d72c7
# ╟─054caa20-795a-11eb-2eba-ddc9b516b4dc
# ╟─72ddc3b0-79c0-11eb-0b19-19eb81a48e41
# ╟─fe09c7b0-e846-11ea-3458-3d72f807da89
# ╟─f8fc0520-699d-11eb-0d45-5932beaa149f
# ╟─2cc390b2-69a0-11eb-2e8d-273d9dd6bee6
# ╟─6bfa7a80-8018-11eb-2c0f-3b411d6291d4
# ╠═63564260-8018-11eb-058f-292b305207bc
# ╟─f6cc0882-716e-11eb-3380-95579d37d720
# ╟─464ce522-717c-11eb-2810-7b244052747d
# ╟─a50f5ac0-8009-11eb-2235-5d7e75f62a1f
# ╟─06cd3970-801a-11eb-1065-49fad5c8ba21
# ╠═0d59e7b0-7e63-11eb-13d2-178ea4c4d09e
# ╟─3b7d82e0-8012-11eb-2bea-fbcfebe6ffe4
# ╟─5e9c7290-8012-11eb-3191-79e66024157d
# ╠═cf170040-e845-11ea-145d-f716ad472b84
# ╟─6f1de950-8012-11eb-399a-432019e28e00
# ╟─e0acaee0-84c1-11eb-245b-4b0042d2642c
# ╟─f300c2c0-84c1-11eb-335f-816c1c15fe2a
# ╠═d3de7090-84cb-11eb-1e8c-635d7db04e99
# ╠═1a829760-84c2-11eb-10d9-adc0e83b8a45
# ╟─80d83300-897a-11eb-369c-e7a5bfa489fa
# ╟─4e1ee580-84ca-11eb-2693-49048d8f1dd3
# ╟─4ec09530-8013-11eb-2870-dd7421db7daa
# ╠═703d4ec0-f728-11ea-150b-e9e101b7acbd
# ╠═6b13e4b0-8010-11eb-0826-f74b42b9c37a
