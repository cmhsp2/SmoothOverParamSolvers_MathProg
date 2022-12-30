% solve min_b |X*b- y|_2^2/(2*lambda)+ |x|_q/q
% q should be in (1/2,2/3]
% lambda >= 0
% Apply varpro, splitting into 3 factors with 2 factors on the outside
%
function [b,f] = func_lq_varpro(X,y,lambda,q)
[m,n] = size(X);
p = q/(2-2*q);
Efun = @(v)  getFG(v,X,y,lambda,p);

v0 = rand()*ones(2*n,1);

warning off
lb = -inf(2*n,1);
ub = inf(2*n,1);
opts    = struct('x0',v0,'printEvery', -1, 'm', 20, 'maxIts', 100 );
opts.factr = 1e-4;
[v, ~, R] = lbfgsb(Efun, lb, ub, opts );

f = R.err(:,1);
V2 = (v(1:n).*v(n+1:2*n)).^2;
alpha = (lambda*speye(m)+ X*spdiags(V2,0,n,n)*X')\(-y);
b = - V2 .*(X'*alpha);

warning on

end

function [fval,Grad] = getFG(v,X,y,lambda,p)
[m,n] = size(X);
v1 = v(1:n);
v2 = v(n+1:2*n);

V2 = (v1.*v2).^2;

alpha = (lambda*speye(m)+ X*spdiags(V2,0,n,n)*X')\(-y);
Xta2 = (X'*alpha).^2;


gradv = [sign(v1).*abs(v1).^(2*p-1); v2];
Swapv = [v2; v1];
Grad = gradv - Swapv.^2.*v.*repmat(sum(Xta2,2),2,1);


normv =   sum(abs(v1).^(2*p))/(2*p)+ sum(v2.^2)/2;
sumv  = @(v) sum(v(:));
fval = normv  ...
    - 1/2*sumv(V2.*Xta2) - sumv(alpha.*y)-lambda*norm(alpha, 'fro')^2/2 ;



end
