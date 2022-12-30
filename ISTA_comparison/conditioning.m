clearvars
% clf
%fourier
n=  100;
fc = 8;
t = linspace(0,1,n);
Xm =  exp(-2*pi*1i*(-fc:fc)'*t)/sqrt(2*fc+1);
X0 = [real(Xm); imag(Xm)];
m = size(X0,1);
X = X0;
L = norm(X)^2;

% X=randn(m,n);

beta0 = zeros(n,1);
supp0 = [n/2,ceil(n/4)];
supp0(:);
beta0(supp0) = [1;-1];
y = X*beta0;

lambda_vals = norm(X'*y,'inf').*linspace(1/2000,1/2,20);%[2000,1500, 1000,500,100,50,10,2];
ratio = [];
for i=1:length(lambda_vals)

lambda = lambda_vals(i);
mfun = @(beta) 0.5*norm(X*beta - y)^2/lambda+norm(beta,1);

%% VarPro
Matv = @(v) lambda*speye(m)+X*spdiags(v.^2,0,length(v),length(v))*X';

a_fun = @(v) -(Matv(v)\y);
mfun_2 =  @(v,a) norm(v)^2/2 - lambda*norm(a)^2/2 - sum(a.*y) - norm(v.*(X'*a))^2/2;
Efun0 = @(v, a) deal(mfun_2(v,a), v - v.*(X'*a).^2);
Efun = @(v) Efun0(v, a_fun(v));

lb = -inf(n,1);
ub = inf(n,1);
v0 = ones(n,1)/sqrt(n);
% %
warning off
opts    = struct('x0',v0,'printEvery', -1, 'm', 20, 'maxIts', 1000 );
opts.factr = 1e-4;
tic
[v, ~, R] = lbfgsb(Efun, lb, ub, opts );
time_varpro = toc;
beta_varpro = -(v.^2).*(X'*a_fun(v));
fvals_varpro = R.err(:,1);



%% Hadamard parameterization


gradv = @(u,v) lambda*v +  u.* X'*( X*(u.*v) - y );
gradu = @(u,v) lambda*u +  v.* X'*( X*(u.*v) - y );
grad = @(u,v) [gradu(u,v); gradv(u,v)  ];
Efun = @(uv) deal(mfun(uv(1:end/2).*uv(end/2+1:end) ), grad(uv(1:end/2),uv(end/2+1:end) ) );


lb = -inf(n*2,1);
ub = inf(n*2,1);
uv0 = randn(n*2,1);

warning off
opts    = struct('x0',uv0,'printEvery', -1, 'm', 20, 'maxIts', 5000 );
opts.factr = 1e-3;
tic
[uv, ~, R] = lbfgsb(Efun, lb, ub, opts );
time_hada = toc;
beta_hada = uv(1:end/2).*uv(end/2+1:end);

fvals = R.err(:,1);

%% alternating minimisation
% lambda= .1;
if 0
u = ones(n,1);
v = ones(n,1);
maxiter = 1000;
XtX = X'*X;
Xty = X'*y;
Mtx = @(v) lambda*speye(n)+XtX.*(v.*v');


Inv = @(v) v.*(X' *((lambda*speye(m) + X*spdiags(v.^2,0, n,n)*X'   )\y));

obj = [];
for i = 1:maxiter
%     u = Mtx(v)\(v.*Xty);
%     v = Mtx(u)\(u.*Xty);
    
    u = Inv(v);
    v = Inv(u);
    obj(i) = mfun(u.*v);
end
semilogy(obj/lambda-min(fvals_varpro))
end

%% plot objective
objbest = min(min(fvals_varpro), min(fvals));
semilogy(linspace(0,time_varpro, length(fvals_varpro)),fvals_varpro - objbest, 'b')
hold on
semilogy(linspace(0,time_hada, length(fvals)),fvals - objbest, 'r')

legend('varpro', 'hadamard')

%% check optimality
eta = X'*(X*beta_hada-y)/lambda;
max(eta)
min(eta)

eta_v = X'*(X*beta_varpro-y)/lambda;
max(eta_v)
min(eta_v)
% 
% supp_c = abs(beta_varpro)<1e-5
% sum(supp_c)
% 1/(1-max(eta_v(supp_c)))

%% check conditioning of Hessian
XtX = X'*X;
u = -v.*(X'*a_fun(v));

eta = X'*(X*(u.*v)-y)/lambda;
A = 1/lambda*(diag(u)*XtX*diag(u))+eye(n);
D =  1/lambda*(diag(v)*XtX*diag(v))+eye(n);
B = 1/lambda*(diag(u)*XtX*diag(v) + lambda*diag(eta));

Hess = [A, B; B', D];

Hess_varPro = A- B*inv(D)*B';
cond_varp = cond(Hess_varPro);
cond_hada = cond(Hess);

ratio(i) = cond_varp/cond_hada;
end
% suppc = abs(v)<1e-5;
% w = 1-abs(eta).^2;
% min(w(suppc))
% min(eig(Hess_varPro))
clf
plot(lambda_vals,ratio,'x')
xlabel('\lambda')
ylabel('COND_{VarPro}/COND_{Hadamard}')