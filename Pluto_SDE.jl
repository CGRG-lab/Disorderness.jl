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

# ╔═╡ a3a0bfee-e845-11ea-004f-df9b39aade97
begin
	using PlutoUI
	using Pkg
	Pkg.activate(pwd())
	Pkg.add(["Plots","PlutoUI","DifferentialEquations","Formatting"])
	using Plots
	using DifferentialEquations
	using Formatting
end

# ╔═╡ 7acb7db0-06ac-11eb-3218-b1d3fd2d72c7
md"Particular solution using Euler method (Cyganowski, 2001)"

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

# ╔═╡ ee6fbf50-e457-11ea-1b58-279baf4f9051
@bind totalTimePower html"<input type=range min=1 max=5 step=1>"

# ╔═╡ 7eec6dc0-e846-11ea-1485-976fd1ec0891
@bind D html"<input type=range min=0.1 max=100 step=0.1>"

# ╔═╡ 8ab53e70-e846-11ea-2a1e-5f96a440c04e
md"Diffusion = $D"

# ╔═╡ 7219ada0-e455-11ea-1c15-ff8d9028e833
@bind F_C html"<input type=range min=0.5 max=10 step=0.5>"

# ╔═╡ e0ee7b20-e455-11ea-1055-f7137767ff32
md"Coulomb friction = $F_C "

# ╔═╡ 87fd02c0-e455-11ea-030d-51d9aded687d
@bind F_ext_to_F_C_ratio html"<input type=range min=0 max=0.95 step=0.05>"

# ╔═╡ f4cc20a0-e844-11ea-1015-8981a146f885
@bind NwPower html"<input type=range min = 1 max = 4 step = 1>"

# ╔═╡ 703d4ec0-f728-11ea-150b-e9e101b7acbd
Nw = Int64(floor(10^NwPower))

# ╔═╡ 8291b0c0-f728-11ea-3d89-abe11a4442e2
md"Nw = $Nw"

# ╔═╡ 24030a60-06ad-11eb-1b35-6d49087b3ab1
Y0 = 0

# ╔═╡ a7306060-e455-11ea-0160-c781aeb2956e
F_ext = F_C*F_ext_to_F_C_ratio;

# ╔═╡ 14da7410-e457-11ea-2fa8-ef26a42c3ef5
md"External force = $F_ext"

# ╔═╡ f0b666c0-e843-11ea-1bde-450568599264
totalTime = Float64(10^totalTimePower);

# ╔═╡ 193acb80-e458-11ea-02df-2360eed46602
md"duration of the sample path = $totalTime"

# ╔═╡ 0619dd80-e844-11ea-3506-81aac6c830e2
dt = Float64(totalTime/(Nw-1))

# ╔═╡ 96faa410-e844-11ea-1522-e12b589c2b7b
traceT = collect(range(1,step = dt,length = Nw));

# ╔═╡ 873ef8e0-e845-11ea-1a24-218c57f5ebe5
@bind frictiontype Select(["Coulomb","viscous"])

# ╔═╡ cf170040-e845-11ea-145d-f716ad472b84
if frictiontype=="viscous"
	drift(v) = -F_C*sign(v);
elseif frictiontype=="Coulomb"
	drift(v) = -F_C*sign(v);
	P0(r,Fc,D) = (Fc^2- r^2)/(-2*r*D); 
	Pst(r,Fc,D,v) = P0(r,Fc,D)*exp.(Fc*v/D).*exp.(-r*v.*sign.(v)/D);
end

# ╔═╡ 84c43400-e457-11ea-35fb-fffce102e3f9
function SDE(traceT,drift,D,Y0)
	b=sqrt(2*D);
    dts = diff(traceT);
    steps = length(dts);
    sigma=sqrt.(dts);  # variance=sigma^2
    dW = sigma.*randn((steps,)) .+ 0
    
    traceY = fill(NaN,size(traceT))
    traceY[1] = Y0;
    y = Y0;
    for i = 2:length(traceT)
        y = y + drift(y)*dts[i-1] + F_ext*dts[i-1] + b*dW[i-1]
        traceY[i] = y;
    end
    return traceY
end

# ╔═╡ fffdf5b0-e458-11ea-36af-5555b6bb783d
function plotSDE(traceT,traceY,Pst,F_ext,F_C,D)
	v = range(minimum(traceY),maximum(traceY),length=5000);
	p1 = plot(traceT,traceY,xlabel = "t",ylabel= "v(t)");
	p2 = histogram(traceY,bins = 50,normalize=:pdf,xlabel = "v",ylabel= "v(t)");
	plot!(p2,v,Pst(F_C,F_ext,D,v))
	# :pdf, :density, :probability or :none
	plot(p1,p2,layout = (1,2),legend= false)
end

# ╔═╡ cff1da70-e846-11ea-106f-9f72f2b84762
traceY = SDE(traceT,drift,D,Y0);

# ╔═╡ fe09c7b0-e846-11ea-3458-3d72f807da89
plotSDE(traceT,traceY,Pst,F_ext,F_C,D)

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
# ╠═a3a0bfee-e845-11ea-004f-df9b39aade97
# ╟─7acb7db0-06ac-11eb-3218-b1d3fd2d72c7
# ╠═fe09c7b0-e846-11ea-3458-3d72f807da89
# ╟─921c7722-06ad-11eb-380e-074417e055ed
# ╠═ceee2ff0-06ac-11eb-1d87-4b7dd46d71ad
# ╟─dced9a40-06ad-11eb-0200-6790a3dee68e
# ╟─b5500040-06ad-11eb-0799-cb644eb1589a
# ╟─24df4420-06ae-11eb-1416-bb77f494eaa2
# ╟─193acb80-e458-11ea-02df-2360eed46602
# ╟─ee6fbf50-e457-11ea-1b58-279baf4f9051
# ╟─8ab53e70-e846-11ea-2a1e-5f96a440c04e
# ╟─7eec6dc0-e846-11ea-1485-976fd1ec0891
# ╟─e0ee7b20-e455-11ea-1055-f7137767ff32
# ╟─7219ada0-e455-11ea-1c15-ff8d9028e833
# ╟─14da7410-e457-11ea-2fa8-ef26a42c3ef5
# ╟─87fd02c0-e455-11ea-030d-51d9aded687d
# ╟─8291b0c0-f728-11ea-3d89-abe11a4442e2
# ╟─f4cc20a0-e844-11ea-1015-8981a146f885
# ╟─0619dd80-e844-11ea-3506-81aac6c830e2
# ╟─703d4ec0-f728-11ea-150b-e9e101b7acbd
# ╟─24030a60-06ad-11eb-1b35-6d49087b3ab1
# ╠═a7306060-e455-11ea-0160-c781aeb2956e
# ╠═f0b666c0-e843-11ea-1bde-450568599264
# ╠═96faa410-e844-11ea-1522-e12b589c2b7b
# ╠═873ef8e0-e845-11ea-1a24-218c57f5ebe5
# ╠═cf170040-e845-11ea-145d-f716ad472b84
# ╠═84c43400-e457-11ea-35fb-fffce102e3f9
# ╠═fffdf5b0-e458-11ea-36af-5555b6bb783d
# ╠═cff1da70-e846-11ea-106f-9f72f2b84762
# ╠═4f235140-e84a-11ea-2bf7-438cfd72dae5
# ╠═4bc12890-e84c-11ea-3cf4-85c3d9d481f0
