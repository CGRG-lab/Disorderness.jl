function circles(h,k,r)
	# see this for more information: https://discourse.julialang.org/t/plot-a-circle-with-a-given-radius-with-plots-jl/23295/5
theta = collect(range(0,2*π,length = 500));
return h .+ r*sin.(theta), k .+ r*cos.(theta);
end