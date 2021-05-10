# to get one realization

function SDE(traceT::Array{Float64,1},drift,D::Union{Float64},Y0::Union{Float64,Array{Float64,1}})
	b=sqrt(2*D);
    dts = diff(traceT);
    steps = length(dts);
    sigma=sqrt.(dts);  # variance=sigma^2
    dim = length(Y0);
    dW = sigma.*randn((steps,dim)) .+ 0;
    if ~isa(Y0, Array)
        Y0 = [Y0];
    end
    if ~isa(drift, Array)
        drift = [drift];
    end
    
    if length(drift)!= dim
        error("Number of `drift` function should be the same as the number of the initial condition `Y0`");
    end

    traceY = fill(NaN,(length(traceT),dim));
    traceY[1,:] = Y0;
    
    for j = 1:dim 
        y = Y0[j];
        for i = 2:length(traceT)
            y = y + drift[j](y)*dts[i-1] + b*dW[i-1, j]
            traceY[i, j] = y;
        end
    end
    return traceY
end

function SDE(dt,drift::Function,energyFunc::Function,D,Y0)
# CHECKPOINT:
# 1. Create Union for ::Function and Array{Function,1}
# 2. energyFunc is for the way calculating energy of an event.
# 3. Try to make it compatible with 2D process

	b=sqrt(2*D);
    dts = diff(traceT);
    steps = length(dts);
    sigma=sqrt.(dts);  # variance=sigma^2
    traceY2 = fill(0,size(traceT))
    for k = 1:numrealization
        dW = sigma.*randn((steps,)) .+ 0
        traceY = fill(NaN,size(traceT))
        traceY[1] = Y0;
        y = Y0;
        for i = 2:length(traceT)
            y = y + drift(y)*dts[i-1] + b*dW[i-1]
            traceY[i] = y;
        end
        traceY2 = traceY .+ traceY2;
    end
    return traceY2
end

# ensemble prediction
function SDE(dt,predictat,numrealizations::Int,drift::Function,D,Y0)
	b=sqrt(2*D);
    steps = floor(predictat/dt) |> Int64;
    dt = predictat/steps; # adjust dt a little bit to make dt*steps exactly the time we want to predict.
    sigma = sqrt(dt);  # variance=sigma^2
    lenY = steps +1;
    ensembleY = [];
    for k = 1:numrealizations
        dW = sigma.*randn((steps,)) .+ 0       
        y = Y0;
        for i = 2:lenY
            y = y + drift(y)*dt + b*dW[i-1]
        end
        push!(ensembleY, y);
    end
    return ensembleY
end

