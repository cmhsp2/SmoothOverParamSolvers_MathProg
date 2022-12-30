% Use Varpro to solve
%        min_beta |D beta|_{1,2} + 1/(2*lambda)*|X*beta - y|_2^2
%
% When y in R^{m times T}, this is group total variation.
% 
% beta in R^{n times T} and |D beta|_{1,2} is a group norm, 
% with groups of size 2*T for isotropic TV and size T for anisotropic mode
% t0 is the starting point is is either size n (for isotropic) or 2*n
% (anisotropic).
function [beta_pro,R,a,v] = func_TV2D_pro(X,y,lambda,proxcalc,D,t0,isotropic)

% VarPro
XtX = X'*X;
Xty = X'*y;
DXty = D*Xty;
DDt = D*D';

if isotropic
    S2 = @(x) sum(x(1:end/2,:).^2 + x(1+end/2:end,:).^2,2);
else
    S2 = @(x) sum(x.^2,2);
end

vec = @(x) x(:);

getab = @(v) get_AB(v,D, DDt, DXty, Xty, XtX,lambda,proxcalc,isotropic);
f = @(v2,a,b) sum(v2)/2 - sum(v2.*S2(a))/2 ...
    + 1/2/lambda * norm(X*b - y, 'fro')^2 ...
    + sum(vec((D*b).*a) );
Gradf = @(v)  getGrad(v,getab,f,S2);


warning off;
lb = -inf(size(t0));
ub = inf(size(t0));
opts    = struct('x0',t0,'printEvery', 1, 'm', 5, 'maxIts', 100 );
% opts.factr = 1e-4;
[v, ~, R] = lbfgsb(Gradf, lb, ub, opts );
v2 = v.^2;
[a,beta_pro] = getab(v2);

warning on;

end






%linear system solve
function [a,b] = get_AB(v2,D, DDt, DXty, Xty, XtX,lambda,proxcalc,isotropic)
if isotropic
    v2 = repmat(v2,2,1);
end
p = length(v2);

if proxcalc
    M = (spdiags(v2,0,p,p) + lambda*DDt);
    a = M\ DXty;
    
    b = Xty - lambda*D'*a;
    
else
    b = (XtX+lambda*(D')*(spdiags(1./v2,0,p,p)*D))\Xty;
    Db = D*b;
    a = Db./v2;
end

end

%get function value and gradient
function [f,g] = getGrad(v,get_ab,f,S2)
v2 = v.^2;

[a,b] = get_ab(v2);
gradF = v - v.*S2(a);


f= f(v2,a,b);
g= gradF;
end
