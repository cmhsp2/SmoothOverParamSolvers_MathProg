
function [x,fvals] = func_ADMM_overlap_group(proxG,D,K,b,gamma,niter,mfunc)

Ktb = K'*b;
y = zeros(size(D,1),size(b,2));
psi = y;

fvals=zeros(niter,1);

[m,n] = size(K);   
W = spdiags(1./diag(D'*D),0,n,n);
P = (gamma*speye(m)+K*(W*K'));

for i=1:niter
    v = ( Ktb + D'*psi + gamma* D'*y);
    
         
    Wv = W*v;
    x =1/gamma* Wv;
    x2 =K'*(P\(K*Wv));
    x = x - 1/gamma*W*x2;    

    y = proxG(D*x - psi/gamma);
    psi = psi + gamma*(y-D*x);
    fvals(i) = mfunc(x);
    
    
    if i>10 && abs(fvals(i)-fvals(i-1))<1e-7
        break
    end
end
fvals = fvals(1:i-1);
% x = x(i2,:);
end


