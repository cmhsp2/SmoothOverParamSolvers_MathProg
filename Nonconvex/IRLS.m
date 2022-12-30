% solve min_beta lambda*|beta|_{2,q}+norm(X*beta - y, 'fro')^2/2
function [x,fval] = IRLS(X,y,lambda,q)

fval = [];
nm = @(x) sqrt(sum(abs(x).^2,2));
mfunc = @(x) lambda*sum(nm(x).^q) + norm(X*x - y, 'fro')^2/2;
[m,n] = size(X);
warning off
eps= 1;
x = pinv(full(X))*y;
kmax = 16;
fval(end+1) = mfunc(x);
tic
for k=1:kmax
    it = 1;
    res = 1;
    while it<10000 && res >sqrt(eps)/100
        xm = x;
        w = (sum(abs(x).^2,2)+eps).^(q/2-1);
        x = (1./w).*X'*((lambda*speye(m)+(X*spdiags(1./w, 0,n,n)*X'))\y);
        res =  norm(xm-x)/norm(x);
        it = it+1;
    end
    fval(end+1) = mfunc(x);
    
    eps = eps/10;
end

