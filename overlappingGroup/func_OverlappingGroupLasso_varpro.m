% Solve min_{beta} |D*beta|_{1,2} + 1/(2*lambda) |X*beta - y|_2^2
% with lambda>0
% The operator D and index set idx are defined s.t.  Dx = (sqrt(n_i)* y_i), 
% where (y_i)_{i=idx(j), ... idx(j+1)-1} are the elements in the jth group of
% x, and n_i = idx(j+1) - idx(j)
%

function [beta_pro,R] = func_OverlappingGroupLasso_varpro(D,X,y,lambda, idx,proxcalc)
XtX = X'*X;
Xty = X'*y;
DXty = D*Xty;
DDt = D*D';


S2 = @(x) Sum2(x,idx);
vec = @(x) x(:);


getab = @(v) get_AB(v,D, DDt, DXty, Xty, XtX,X,lambda,proxcalc,idx);
f = @(v2,a,b) sum(v2)/2 - sum(v2.*S2(a))/2 ...
    + 1/2/lambda * norm(X*b - y, 'fro')^2 ...
    + sum(vec((D*b).*a) );
Gradf = @(v)  getGrad(v,getab,f,S2);


t0 = randn(length(idx)-1,1);
tic
warning off;
lb = -inf(size(t0));
ub = inf(size(t0));
opts    = struct('x0',t0,'printEvery', 1, 'm', 5, 'maxIts', 100 );
% opts.factr = 1e-4;
[v, ~, R] = lbfgsb(Gradf, lb, ub, opts );
v2 = v.^2;
[a,beta_pro] = getab(v2);

warning on

end


function y = Sum2(x,idx)
x = x.^2;
y = zeros(length(idx)-1,1);
for i=1:length(idx)-1
    y(i) = sum(x(idx(i):idx(i+1)-1));
end
end


function [a,b] = get_AB(v20,D, DDt, DXty, Xty, XtX,X,lambda,proxcalc, Idx)

v2 = zeros(length(Idx)-1,1);
for i=1:length(Idx)-1
    v2(Idx(i):Idx(i+1)-1) = v20(i);
end
p = length(v2);

if proxcalc
    M = (spdiags(v2,0,p,p) + lambda*DDt);
    a = M\ DXty;
    b = Xty - lambda*D'*a;
else
    
    [m,n] = size(X);
    
    W = spdiags(1./diag((D')*(spdiags(1./v2,0,p,p)*D)),0,n,n);
    WXty = W*Xty;
    b =1/lambda* WXty;
    b2 =X'*((lambda*speye(m)+X*(W*X'))\(X*WXty));
    b = b - 1/lambda*W*b2;    
    
    Db = D*b;
    a = Db./v2;
    
    
end

end

function [f,gradF] = getGrad(v,get_ab,f,S2)
v2 = v.^2;
[a,b] = get_ab(v2);

gradF = v - v.*S2(a);

f = f(v2,a,b);
end