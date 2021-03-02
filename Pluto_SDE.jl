### A Pluto.jl notebook ###
# v0.12.21

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
	using Plots
	using PlutoUI
	using DifferentialEquations
	using Formatting
	using LaTeXStrings
	using Printf
end

# ╔═╡ 366389d0-68e3-11eb-19dd-d78bd237ab71
PlutoUI.TableOfContents()

# ╔═╡ b4359b10-7a26-11eb-151d-f95a6b998a79
md"# CHECKPOINT

read array slicing
- https://discourse.julialang.org/t/trouble-understanding-slicing/36821/5
read julia debugger (learn how to use @run and @enter)
"

# ╔═╡ fdcd11d0-7a8b-11eb-191b-316fdb08c3d9
function circles(h,k,r)
	# see this for more information: https://discourse.julialang.org/t/plot-a-circle-with-a-given-radius-with-plots-jl/23295/5
theta = collect(range(0,2*π,length = 500));
return h .+ r*sin.(theta), k .+ r*cos.(theta);
end

# ╔═╡ 7acb7db0-06ac-11eb-3218-b1d3fd2d72c7
md"Particular solution using Euler method (Cyganowski, 2001)"

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

# ╔═╡ 921c7722-06ad-11eb-380e-074417e055ed
md"Particular solution using DifferentialEquation.jl"

# ╔═╡ ceee2ff0-06ac-11eb-1d87-4b7dd46d71ad
# begin
# 	f(u,p,t) = -F_C*sign(u) + F_ext;
# 	g(u,p,t) = sqrt(D);
# 	tspan = (0.0,totalTime);
# 	prob = SDEProblem(f,g,Y0,tspan);
# 	sol = solve(prob,alg,dt=dt,adaptive=false);
# 	p1 = plot(sol,xlabel="t",ylabel="v(t)");
# 	v = range(minimum(sol.u), length=1000, stop=maximum(sol.u));
# 	p2 = histogram(sol.u,bins = 50,normalize=:pdf,xlabel = "v",ylabel= "v(t)");
# 	plot!(p2,v,Pst(F_C,F_ext,D,v));
# 	plot(p1,p2,layout = (1,2),legend= false);
# 	# or plotSDE(sol.t,sol.u,Pst,F_ext,F_C,D)
# end

# ╔═╡ dced9a40-06ad-11eb-0200-6790a3dee68e
md"select the algorithm for solving the SDE problem:"

# ╔═╡ b5500040-06ad-11eb-0799-cb644eb1589a
@bind alg_str Select(["SRIW1","EM"])

# ╔═╡ 24df4420-06ae-11eb-1416-bb77f494eaa2
begin
	algs_str = ["SRIW1","EM"];
	algs = [SRIW1(),EM()];
	LocB = indexin(algs_str,[alg_str]);
	LocB[isnothing.(LocB)].=0;
	Lia = LocB.!=0;
	alg = algs[Lia][1]
end

# ╔═╡ 703d4ec0-f728-11ea-150b-e9e101b7acbd
Nw = Int64(floor(10^NwPower));

# ╔═╡ 24030a60-06ad-11eb-1b35-6d49087b3ab1
Y0 = 0.0;

# ╔═╡ a7306060-e455-11ea-0160-c781aeb2956e
F_ext = Float64(F_C*F_ext_to_F_C_ratio);

# ╔═╡ 14da7410-e457-11ea-2fa8-ef26a42c3ef5
md"External force = $F_ext"

# ╔═╡ f0b666c0-e843-11ea-1bde-450568599264
totalTime = Float64(10^totalTimePower);

# ╔═╡ 193acb80-e458-11ea-02df-2360eed46602
md"duration of the sample path = $totalTime"

# ╔═╡ 0619dd80-e844-11ea-3506-81aac6c830e2
dt = Float64(totalTime/(Nw-1))

# ╔═╡ 2cc390b2-69a0-11eb-2e8d-273d9dd6bee6
# noted that the space after "$" is very important!
let
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

# ╔═╡ 92f12440-79c4-11eb-3878-eda56482eec0
traceT = collect(range(1,step = dt,length = Nw));

# ╔═╡ 72ddc3b0-79c0-11eb-0b19-19eb81a48e41
let
	include(joinpath(srcdir,"SDE.jl"))
	titlefSz = 11;
	drift_x(v) = -F_C*v+F_ext;
	drift_y(v) = -F_C*v;
	traceY = SDE(traceT, [drift_x, drift_y], D, Y0; dim=2);
	absY = traceY[:,1].^2 .+ traceY[:,2].^2;
	p1 = plot(traceT, traceY[:,1], xlabel = "", ylabel="v_x(t)", 
		titlefontsize = titlefSz,
		title="sample path (1D)\n in velocity space",legend = false);
	p2 = plot(traceT, traceY[:,2], xlabel = "t",ylabel="v_y(t)",legend = false);
	p3 = plot(traceY[:,1], traceY[:,2], seriestype=:path, # default: seriestype=:line
		xlabel="v_x", ylabel="v_y", titlefontsize = titlefSz,
		title="sample path (2D)\n in velocity space"); 
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
	plot(p1,p2,p3,p4,layout = lout1, size=(700,300));
end

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
	
	titlefSz = 11;
	v = range(minimum(traceY),maximum(traceY),length=5000);
	p1 = Plots.plot(traceT,traceY,xlabel = "t",ylabel= "v(t)");
	p2 = histogram(traceY,bins = 50,normalize=:pdf,xlabel = "v",ylabel= "v(t)");
	plot!(p2,v,Pst(F_C,F_ext,D,v))
	# :pdf, :density, :probability or :none
	p3 = histogram(ensembleY, bins=50, normalize =:pdf, xlabel = "v", 
		ylabel = "v(t)");
	v2 = range(minimum(ensembleY),maximum(ensembleY),length=5000);
	Plots.plot!(p3,v2,Pst(F_C,F_ext,D,v2));
	Plots.plot(p1,p2,p3,layout = (1,3),legend= false,
	title = ["one sample path" "time averaged\n distribution" "ensemble averaged\n distribution"],titlefontsize=titlefSz, size = (700, 200))
end

# ╔═╡ f6cc0882-716e-11eb-3380-95579d37d720
md"
### Numerically solving the Stochastic Differential Equation (SDE)
"

# ╔═╡ 454dfea0-716f-11eb-1736-a70c4949c21b
L"
\mathrm{d}v(t) = \text{drift}(t) + \sqrt{2D} \mathrm{d}W
"

# ╔═╡ 464ce522-717c-11eb-2810-7b244052747d
md"
See if 

$v(t) = \mathrm{e}^{-t/\tau_B}v(0) + \frac{1}{m}\int_0^t \mathrm{e}^{-(t-s)/\tau_B}\mathrm{d} W(s)$ 

increase the efficiency of calculating SDE.
This equation comes from Eq. 6.19 from [here](http://physics.gu.se/~frtbm/joomla/media/mydocs/LennartSjogren/kap6.pdf).
"

# ╔═╡ 4f235140-e84a-11ea-2bf7-438cfd72dae5
# function SDE2(F_C,F_ext,D,dt,totalTime)
# 	u₀=0.5;
# 	f(u,p,t) = -F_C*sign(u) + F_ext
# 	g(u,p,t) = sqrt(D);
# #	dt = 1//2^(4)
# 	tspan = (0,totalTime)
# 	prob = SDEProblem(f,g,u₀,tspan)
# 	sol = solve(prob,EM(),dt=dt);
# end

# ╔═╡ 4bc12890-e84c-11ea-3cf4-85c3d9d481f0
# sol = SDE2(F_C,F_ext,D,dt,totalTime);

# ╔═╡ Cell order:
# ╠═366389d0-68e3-11eb-19dd-d78bd237ab71
# ╠═5b4cd610-68df-11eb-148b-8f94bd42f945
# ╠═b4359b10-7a26-11eb-151d-f95a6b998a79
# ╠═a3a0bfee-e845-11ea-004f-df9b39aade97
# ╠═fdcd11d0-7a8b-11eb-191b-316fdb08c3d9
# ╟─7acb7db0-06ac-11eb-3218-b1d3fd2d72c7
# ╟─72ddc3b0-79c0-11eb-0b19-19eb81a48e41
# ╟─fe09c7b0-e846-11ea-3458-3d72f807da89
# ╟─054caa20-795a-11eb-2eba-ddc9b516b4dc
# ╟─f8fc0520-699d-11eb-0d45-5932beaa149f
# ╟─2cc390b2-69a0-11eb-2e8d-273d9dd6bee6
# ╟─921c7722-06ad-11eb-380e-074417e055ed
# ╠═ceee2ff0-06ac-11eb-1d87-4b7dd46d71ad
# ╟─dced9a40-06ad-11eb-0200-6790a3dee68e
# ╟─b5500040-06ad-11eb-0799-cb644eb1589a
# ╟─24df4420-06ae-11eb-1416-bb77f494eaa2
# ╟─193acb80-e458-11ea-02df-2360eed46602
# ╟─14da7410-e457-11ea-2fa8-ef26a42c3ef5
# ╠═0619dd80-e844-11ea-3506-81aac6c830e2
# ╠═703d4ec0-f728-11ea-150b-e9e101b7acbd
# ╠═24030a60-06ad-11eb-1b35-6d49087b3ab1
# ╠═a7306060-e455-11ea-0160-c781aeb2956e
# ╠═f0b666c0-e843-11ea-1bde-450568599264
# ╠═92f12440-79c4-11eb-3878-eda56482eec0
# ╠═cf170040-e845-11ea-145d-f716ad472b84
# ╠═f6cc0882-716e-11eb-3380-95579d37d720
# ╠═454dfea0-716f-11eb-1736-a70c4949c21b
# ╠═464ce522-717c-11eb-2810-7b244052747d
# ╠═4f235140-e84a-11ea-2bf7-438cfd72dae5
# ╠═4bc12890-e84c-11ea-3cf4-85c3d9d481f0
