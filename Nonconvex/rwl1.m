% solve min_beta lambda*|beta|_{2,q}+norm(X*beta - y, 'fro')^2/2


function [x,fval] = rwl1(X,y,lambda,q)


fval = [];

nm = @(x) sqrt(sum(abs(x).^2,2));
mfunc = @(x) lambda*sum(nm(x).^q) + norm(X*x - y, 'fro')^2/2;
[m,n] = size(X);
eps= 1;
x = func_l1_varpro(X,y,lambda);
kmax = 5;
fval(end+1) = mfunc(x);

for k=1:kmax
    it = 1;
    res = 1;
    while it<20 && res >0.001
        xm = x;
        
        w = (sqrt(sum(abs(x).^2,2))+eps).^(1-q);
        b = func_l1_varpro(X*spdiags(w, 0, n,n),y,lambda);
        x = w.*b;
        res =  norm(xm-x)/norm(x);
        it = it+1;
    end
    fval(end+1) = mfunc(x);
    
    eps = eps/10;
end
