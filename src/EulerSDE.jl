function EulerSDE(traceT,drift,D,Y0)
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