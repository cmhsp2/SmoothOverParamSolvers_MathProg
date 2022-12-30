%solve multitask lasso, min_b  |b|_{1,2} + 1/(2*lambda) |X*b - y|_F^2
function [b,R,alpha,fval] = func_l1_varpro(X,y,lambda,v0, XtX, Xty)
[m,n] = size(X);
if nargin<4
    v0 =  randn(n,1)*0.1;
end


sumv = @(x) sum(x(:));

% getalpha = @(v)(lambda*speye(m)+ X*spdiags(v.^2,0,n,n)*X')\(-y);
nmy = norm(y(:))^2;
if m<n
    Efun = @(v) getGF(v, lambda, y, X, [], [],nmy);
else    
    if nargin<6
     XtX = X'*X;
     Xty = X'*y;
    end
     Efun = @(v) getGF(v, lambda, y, X, XtX, Xty,nmy);

end

options.niter = 5000;
options.bfgs_memory = 20;
options.verb = 0;

%compute l1 solution

[v,R] =  perform_bfgs(Efun , v0, options);

v2 = v.^2;

if m<n
    alpha =  (X*spdiags(v2,0,n,n)*X' + lambda * speye(m))\(-y);
    b = - v2 .*(X'*alpha);
    
	fval =  sum(v2)/2  ...
    - 1/2*sumv(v2.*(X'*alpha).^2) - sumv(alpha.*y)-lambda*norm(alpha(:))^2/2;
    fba = [fval;b(:);alpha(:)];
else    

    
    u = (lambda*speye(n)+spdiags(v,0,n,n)*XtX*spdiags(v,0,n,n))\(v.*Xty);
    b = v.*u;
%     alpha = X*b-y;
    Xta = (XtX*b - Xty)/lambda;
    S = sum(Xty(:).*b(:));
	fval =  norm(v)^2/2  ...
    - 1/2*sumv(v2.*Xta.^2) - (S - nmy)/lambda -(b'*(XtX*b) - 2*S + nmy)/lambda/2;
    fba = [fval;b(:)];
    alpha = [];
end



end

function [fval,grad] = getGF(v, lambda, y, X, XtX, Xty,nmy)
[m,n] = size(X);
v2 = v.^2;
sumv = @(x) sum(x(:));
if m<n
    alpha =  (X*spdiags(v2,0,n,n)*X' + lambda * speye(m))\(-y);
    grad =  v - v.*sum((X'*alpha).^2,2) ;
	fval =  sum(v2)/2  ...
    - 1/2*sumv(v2.*(X'*alpha).^2) - sumv(alpha.*y)-lambda*norm(alpha(:))^2/2;

else    

    
    u = (lambda*speye(n)+spdiags(v,0,n,n)*XtX*spdiags(v,0,n,n))\(v.*Xty);
    b = v.*u;
    Xta = (XtX*b - Xty)/lambda;
    S = sum(Xty(:).*b(:));
    nma2 = (b'*(XtX*b) - 2*S + nmy)/lambda^2;
    ya = (S - nmy)/lambda;
    
    grad =  v - v.*sum((Xta).^2,2) ;
	fval =  sum(v2)/2  ...
    - 1/2*sumv(v2.*Xta.^2) - ya-lambda*nma2/2;

end

end