% solve min_b |X*b- y|_2^2/(2*lambda)+ |x|_q/q
% q should be in (1/2,2/3]
% lambda >= 0
% Apply varpro, splitting into 3 factors with 2 factors on the inside
% So the inner solver is varpro on l1 
function [x_u,f] = func_lq_nonlin_varpro(X,y,lambda, q)

[m,n] = size(X);
p = q/(2-2*q);
Efun = @(u) getFG(u, X,y,lambda,p);
lb = -inf(n,1);
ub = inf(n,1);
u0 = rand()*ones(n,1);
% %
warning off
opts    = struct('x0',u0,'printEvery', -1, 'm', 20, 'maxIts', 100 );
% opts.factr = 1e-4;
[u, ~, R] = lbfgsb(Efun, lb, ub, opts );

b = func_l1_varpro(X*spdiags(u,0,n,n),y,lambda);
x_u = u.*b;
warning on
f =  R.err(:,1);


end

function [f,Grad] = getFG(u, X,y,lambda,p)
[m,n] = size(X);

[b,~,a,fval] = func_l1_varpro(X*spdiags(u,0,n,n),y,lambda);
f = sum(abs(u).^(2*p))/(2*p) + fval;

Grad = sign(u).*abs(u).^(2*p-1)  + sum(b.*(X'*a),2);

end
