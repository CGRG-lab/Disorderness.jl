
	u₀=0
	f(u,p,t) = -F_C*sign(u)+F_ext;
	g(u,p,t) = sqrt(2*D);
	# dt = 1//2^(4)
	tspan = (0.0,totalTime)
	prob = SDEProblem(f,g,u₀,tspan)
	sol = solve(prob,SRIW1(),dt=dt)
	
	# single sample path
	v = range(minimum(sol.u),maximum(sol.u),length=5000);
	p2 = histogram(sol.u, bins = 50, normalize = :pdf, xlabel = "v", ylabel = "P(v)");
	p1 = plot(sol.t,sol.u,xlabel = "t",ylabel= "v(t)");
	plot!(p2, v,Pst(F_C,F_ext,D,v))
	
	# ensemble prediction
	using DifferentialEquations.EnsembleAnalysis
	ensembleprob = EnsembleProblem(prob);
	solE = solve(ensembleprob,SRIW1(),EnsembleThreads();
		trajectories=1000);#,dt=totalTime/50) # EM(),EnsembleThreads() is optional
	summ = EnsembleSummary(solE, 0:0.5:totalTime)
	plot(summ)
	
	ensembleu = map(x-> x.u[end],solE);
	p3 = histogram(ensembleu, bins = 50, normalize = :pdf, xlabel = "v", ylabel = "P(v)");
	v = range(minimum(ensembleu),maximum(ensembleu),length=5000);
	plot!(p3, v,Pst(F_C,F_ext,D,v))
	plot(p1,p2,p3,layout = (1,3), size= (700,200))
	
	@time solve(ensembleprob,SRIW1(),EnsembleThreads(); trajectories=1000, dt=totalTime/50);
	@time solve(ensembleprob,SRIW1(); trajectories=1000, dt=totalTime/50);
		