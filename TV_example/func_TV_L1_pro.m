% Use Varpro to solve
%        min_beta |D beta|_{1,2} + 1/(2*lambda)*|X*beta - y|_1
%
% When y in R^{m times T}, this is group total variation.
% 
% beta in R^{n times T} and |D beta|_{1,2} is a group norm, 
% with groups of size 2*T for isotropic TV and size T for anisotropic mode
% t0 is the starting point is is either size n (for isotropic) or 2*n
% (anisotropic).
function [x,R] = func_TV_L1_pro(A,y,lambda,D,t0,isotropic, Aeye)



[m,n] = size(A);
Gradf = @(wv)  getFGrad(wv, D, A,y,isotropic,lambda,Aeye);

warning off;
lb = -inf(size(t0));
ub = inf(size(t0));
opts    = struct('x0',t0,'printEvery', 1, 'm', 5, 'maxIts', 100 );
% opts.factr = 1e-4;
[wv, ~, R] = lbfgsb(Gradf, lb, ub, opts );
w = wv(1:m);
v = wv(m+1:end);
v2 = v.^2;
w2 = w.^2;
[xi,alpha,x] = linearSolve(v2,w2, D, A,y,isotropic,lambda,Aeye);


warning on;

end






%linear system solve
function [xi,alpha,x] = linearSolve(v2,w2, D, A,y,isotropic,lambda, Aeye)
if isotropic
    v2 = repmat(v2,2,1);
end
p = length(v2);
[m,n] = size(A);

if ~Aeye

M = [spdiags(v2,0,p,p), sparse(p,m), D; sparse(m,p), lambda*spdiags(w2,0,m,m), A; ...
    D', A', sparse(n,n)];
[m,d] = size(y);
b = sparse(p+n+m,d);
b(p+1:p+m,:) = y;

res = M\b;
alpha = res(1:p,:);
xi = res(p+1:p+m,:);
x = res(p+m+1:end,:);
% 
% M = D'*spdiags(1./v2,0,p,p)*D +1/lambda* A'*spdiags(1./w2,0,m,m)*A;
% x = 1/lambda*(M\(A'*((1./w2).*y)));
% alpha = -(1./v2).*(D*x);
% xi = -(1./w2).*(A*x-y);

else
    
    M = lambda*D*spdiags(w2,0,m,m)*D'+spdiags(v2,0,p,p);
    alpha = -M\(D*y);
    xi = -D'*alpha;
    x = y -lambda*w2.*xi;
    
end



end


function [fval,grad] = getFGrad(wv, D, A,y,isotropic,lambda,Aeye)
if isotropic
    S2 = @(x) sum(x(1:end/2,:).^2 + x(1+end/2:end,:).^2,2);
else
    S2 = @(x) sum(x.^2,2);
end

[m,d] = size(y);
w = wv(1:m);
v = wv(m+1:end);
v2 = v.^2;
w2 = w.^2;
[xi,alpha,x] = linearSolve(v2,w2, D, A,y,isotropic,lambda,Aeye);


% S = @(x) sum(sqrt(sum(x(1:end/2,:).^2+x(end/2+1:end,:).^2,2)));
% fval =  S(D*x) + sum(sqrt(sum((A*x - y).^2, 2)));
fval = sum(v2)/2+sum(w2)/2/lambda - sum(v2.*S2(alpha))/2 - lambda*sum(w2.*sum(xi.^2,2))/2 + sum(xi(:).*y(:));
grad = [w/lambda-lambda*w.*sum(xi.^2,2); v-v.*S2(alpha)];

end
