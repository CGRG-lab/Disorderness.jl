# example
# numcars = 5;
# x = collect(range(0,10,length = 100));	
# ys = randboxcars(x, numcars)
# boxwidthstd = 0.4 unit of x's length
function randboxcars(x, numcars::Int; boxwidthstd = 0.4)
	halfwidthstd = boxwidthstd*0.5;
	boxcentr = rand(x, numcars);
	halfboxwidth = halfwidthstd .* (x[end]-x[1]) .* abs.(randn(numcars));
	ys = fill(0, length(x));
	for i = 1:numcars
		lb, centr, rb = get_rect(x,boxcentr,halfboxwidth,i);
		ys[lb .< x .< rb] .= ys[lb .< x .< rb] .+ 1;
	end
	return ys
end
# number of rows of ys has to be the same as length(x).
function randboxcars(x, ys::Array; boxwidthstd = 0.4)
	halfwidthstd = boxwidthstd*0.5;
	numcars = size(ys, 2);
	boxcentr = rand(x, numcars);
	halfboxwidth = halfwidthstd .* (x[end]-x[1]) .* abs.(randn(numcars));
	ys = fill(0, (length(x), numcars));
	for i = 1:numcars
		lb, centr, rb = get_rect(x,boxcentr,halfboxwidth,i);
		ys[lb .< x .< rb,i] .= 1;
	end
	return ys
end

function get_rect(x,boxcentr,halfboxwidth,i)
	centr = findlast(x .<= boxcentr[i]);
	rb = x[centr] + halfboxwidth[i]; # right bound
	lb = x[centr] - halfboxwidth[i]; # left bound
	return lb, centr, rb
end