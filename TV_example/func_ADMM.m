%%ADMM
%	Solve min_{x,z} max_psi 1/2 |Kx - b|_2^2 + G(z) s.t z = D*x
%   Inputs:     parameter, gamma >0
%               matrices, K, D
%               maximum number of iterations, niter
%               objective function, mfunc
%               record_interval, how often to record the function value
%               proxG, proximal operator of G
%   Outputs: the minimiser x, function values fvals

function [x,fvals] = func_ADMM(proxG,D,K,b,gamma,niter,mfunc,record_interval)
DtD = (D'*D);
KtK = K'*K;
Ktb = K'*b;
M = KtK +gamma* DtD;

if nargin<8
    record_interval = 1;
end



idx = symrcm(M);
R = chol( M(idx,idx));

[~,i2 ]= sort(idx);
%     R = chol(M);

y = zeros(size(D,1),size(b,2));
psi = y;

fvals=zeros(niter,1);


D = D(:,idx);
Ktb = Ktb(idx,:);
r = 1;
for i=1:niter
    v = ( Ktb + D'*psi + gamma* D'*y);
    x = R\(R'\v);
    y = proxG(D*x - psi/gamma);
    psi = psi + gamma*(y-D*x);
    
    if mod(i,record_interval)==0
    fvals(r) = mfunc(x(i2,:));
    r = r+1;
    end
    
    
    if i>10 && abs(fvals(i)-fvals(i-1))<1e-12
        break
    end
end
fvals(r) = mfunc(x(i2,:));

fvals = fvals(1:r);
x = x(i2,:);

end