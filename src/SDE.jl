# to get one realization
function SDE(traceT,drift,D,Y0; dim=1)
	b=sqrt(2*D);
    dts = diff(traceT);
    steps = length(dts);
    sigma=sqrt.(dts);  # variance=sigma^2
    dW = sigma.*randn((steps,dim)) .+ 0;
    
    if dim == 1
        drift = [drift];
    end
    traceY = fill(NaN,(length(traceT),dim));
    traceY[1,:] .= Y0;
    
    for j = 1:dim 
        y = Y0;
        for i = 2:length(traceT)
            y = y + drift[j](y)*dts[i-1] + b*dW[i-1, j]
            traceY[i, j] = y;
        end
    end
    return traceY
end

function SDE(traceT,drift::Function,D,Y0,numrealization)
# NOT FINISHED YET
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

