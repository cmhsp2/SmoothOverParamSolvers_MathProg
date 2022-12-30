function [x,fvals] = func_PrimalDual(D,proxG,proxFs,sig,tau,niter,mfunc,record_interval,init)

[p2,p] = size(D);
disp('running Primal dual ...')

xi = init.xi;%zeros(p2,1)*0.01;
x = init.x;%zeros(p,1)*0.01;
barx = x;

theta = 1;
fvals = [];
fvals(1) = mfunc(x);

for its = 1:niter

    xi = proxFs(xi + sig*D*barx);
    xold = x;
    x = proxG(x - tau*D'*xi);
    barx = x + theta*(x  - xold);
    
    if mod(its,record_interval)==0
        fvals(end+1) = mfunc(x);        
    end
    
    
    n_tmp = length(fvals);
    if n_tmp>3 && abs(fvals(n_tmp)-fvals(n_tmp-1))<1e-7
        break
    end
    
    
end

fvals(end+1) = mfunc(x);


end

