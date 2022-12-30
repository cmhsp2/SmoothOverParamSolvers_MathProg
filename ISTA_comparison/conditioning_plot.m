%plot the condition number of nabla^2 f(v) at the optimum v, for various
%ill-conditioned Fourier data matrices

clearvars
% clf


fc = 15;
J=8;
jmin = 4;
C_X = [];
C_V = [];
for j=jmin:J
    n = 2^j;
    
    %data matrix
    t = 0:1/n:fc/n;
    X =  exp(-2*pi*1i*(0:fc)'*t)/sqrt(2*fc+1);
    m = size(X,1);
    L = norm(X)^2;
    
    C_X(end+1) =cond(X);
    % X=randn(m,n);
    
    beta0 = zeros(fc+1,1);
    supp0 = [ceil(fc/2),1];
    supp0(:);
    beta0(supp0) = [1;1];
    y = X*beta0;
    
    
    lambda =100/n^3; 
    mfun = @(beta) 0.5*norm(X*beta - y)^2/lambda+norm(beta,1);
    
    %% VarPro
    Matv = @(v) lambda*speye(m)+X*spdiags(v.^2,0,length(v),length(v))*X';
    
    a_fun = @(v) -(Matv(v)\y);
    mfun_2 =  @(v,a) norm(v)^2/2 - lambda*norm(a)^2/2 - real(sum(a'*y)) - norm(v.*abs(X'*a))^2/2;
    Efun0 = @(v, a) deal(mfun_2(v,a), v - v.*abs(X'*a).^2);
    Efun = @(v) Efun0(v, a_fun(v));
    
    lb = -inf(fc+1,1);
    ub = inf(fc+1,1);
    v0 = ones(fc+1,1)/sqrt(fc+1);
    % %
    warning off
    opts    = struct('x0',v0,'printEvery', -1, 'm', 20, 'maxIts', 1000 );
    opts.factr = 1e-4;
    tic
    [v, ~, R] = lbfgsb(Efun, lb, ub, opts );
    time_varpro = toc;
    beta_varpro = -(v.^2).*(X'*a_fun(v));
    fvals_varpro = R.err(:,1);
    
    %% conditioning of Hessian
    XtX = X'*X;
    a = a_fun(v);
    u = -v.*(X'*a);
    eta = X'*(X*(u.*v)-y)/lambda;
    A = 1/lambda*(diag(u)*XtX*diag(u))+eye(fc+1);
    D =  1/lambda*(diag(v)*XtX*diag(v))+eye(fc+1);
    B = 1/lambda*(diag(u)*XtX*diag(v) + lambda*diag(eta));
    
    Hess = [A, B; B', D];
    
    Hess_varPro = A- B*inv(D)*B';
    cond_varp = cond(Hess_varPro);
    
    C_V(end+1) = cond_varp;
end
%%
clf
semilogy((jmin:J),C_X,'-x','linewidth',2)
hold on
semilogy((jmin:J),C_V,'-o','linewidth',2)
xlabel('log(N)', 'fontsize', 18)
ylabel('Condition Number', 'fontsize', 18)
ax = gca;
ax.FontSize = 18;
xlim([jmin,J])
ylim([0,1e21])

set(gca,'fontname','times')  % Set it to times

