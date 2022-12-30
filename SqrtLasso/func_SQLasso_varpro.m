function [beta_varpro, fval_varpro] = func_SQLasso_varpro(X,y,lam)
[m,n] = size(X);
nmy = norm(y)^2;
if m<n
    Efun = @(w) getFG(w, lam, y, X, [], [],nmy);
else    
    XtX = X'*X;
    Xty = X'*y;
    Efun = @(w) getFG(w, lam, y, X, XtX, Xty,nmy);
end



lb = -inf(n+1,1);
ub = inf(n+1,1);
xinit = ones(n+1,1);
opts    = struct('x0',xinit,'printEvery', 1, 'm', 20, 'maxIts', 1000 ); % opts.factr = 1e-4;
[w, ~, R] = lbfgsb(Efun, lb, ub, opts );

v = w(1:end-1);
z = w(end);
if m<n
    alpha = (X*spdiags(v.^2,0,n,n)*X' + lam*z^2 * speye(m))\(-y);
    beta_varpro = -v.^2.*(X'*alpha);

else        
    u = (lam*z^2*speye(n)+spdiags(v,0,n,n)*XtX*spdiags(v,0,n,n))\(v.*Xty);
    beta_varpro = v.*u; 
end


fval_varpro = R.err(:,1);



end

function [fval,grad] = getFG(w, lam, y, X, XtX, Xty,nmy)
v = w(1:end-1);
v2 = v.^2;
z = w(end);
[m,n] = size(X);

if m<n
    alpha =  (X*spdiags(v2,0,n,n)*X' + lam*z^2 * speye(m))\(-y);
    Xta = X'*alpha;
    Xta2 = Xta.^2;
    nma2 = sum(alpha.^2);
    grad = [v - v.*Xta2; z/lam - z*nma2*lam];
    fval =  sum(v2)/2 +z^2/2/lam - sum(v2.*Xta2)/2 - lam*z^2*nma2/2 - sum(y.*alpha);

    
    
else    
    u = (lam*z^2*speye(n)+spdiags(v,0,n,n)*XtX*spdiags(v,0,n,n))\(v.*Xty);
    b = v.*u;
%     alpha = (X*b - y)/(lam*z^2);
    Xta = (XtX*b - Xty)/(lam*z^2);
    Xta2 = Xta.^2;
    S = sum(Xty.*b);
    nma2 = (b'*(XtX*b) - 2*S + nmy)/(lam*z^2)^2;
    ya = (S - nmy)/(lam*z^2);
    grad = [v - v.*Xta2; z/lam - z*nma2*lam];
    fval =  sum(v2)/2 +z^2/2/lam - sum(v2.*Xta2)/2 - lam*z^2*nma2/2 - ya;

end
end