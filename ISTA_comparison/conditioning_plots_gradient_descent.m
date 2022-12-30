%plot the condition number of the inner linear system and nabla^2 f(v) as
%we iterate through gradient descent (for an ill-conditioned Fourier data
%matrix

clearvars
clf
%fourier
fc = 15;
J=8;
jmin = 4;
for j =jmin:J
    %define data matrix
    n = 2^j;
    t = 0:1/n:fc/n;
    X =  exp(-2*pi*1i*(0:fc)'*t)/sqrt(2*fc+1);
    L = norm(X)^2;
    m = size(X,1);
    
    
    %true signal (2 sparse)
    beta0 = zeros(fc+1,1);
    supp0 = [ceil(fc/2),1];
    supp0(:);
    beta0(supp0) = [1;1];
    y = X*beta0;
    
    
    %regularisaton parameter (scaling to ensure we recover 2 spikes)
    lambda = 100/n^3; 
    mfun = @(beta) 0.5*norm(X*beta - y)^2/lambda+norm(beta,1);
    
    % VarPro
    Matv = @(v) lambda*speye(m)+X*spdiags(v.^2,0,length(v),length(v))*X';
    a_fun = @(v) -(Matv(v)\y);
    obj_fn =  @(v,a) norm(v)^2/2 - lambda*norm(a)^2/2 - sum(a'*y) - norm(v.*abs(X'*a))^2/2;
    
    v0 = ones(fc+1,1)/sqrt(fc+1);
    v = v0;
    maxiter = 10001;
    tau = .5; %gradient descent stepsize
    fvals = [];
    cond_in = []; %condition number of inner problem
    cond_varp = []; %condition number of outer problem
    XtX = X'*X;
    err = [];
    for i =1:maxiter
        M = Matv(v);
        %record condition number of inner system
        cond_in(end+1) = cond(M);
        a = -(M\y);
        grad = v - v.*abs(X'*a).^2;
        
        
        %compute Hessian
        u = -v.*(X'*a);
        eta = X'*(X*(u.*v)-y)/lambda;
        A = 1/lambda*(diag(u)*XtX*diag(u))+eye(fc+1);
        D =  1/lambda*(diag(v)*XtX*diag(v))+eye(fc+1);
        B = 1/lambda*(diag(u)*XtX*diag(v) + lambda*diag(eta));
        Hess_varPro = A- B*inv(D)*B';
        cond_varp(end+1) = cond(Hess_varPro);
        
        
        fvals(end+1) = obj_fn(v,a);
        
        %gradient descent update
        v = v-tau*grad;
    end
    beta_varpro = -(v.^2).*(X'*a_fun(v));
    figure(1)
    % clf
    
    mycolor = [1,0,0]*(jmin/j);
    semilogy(cond_in,  '-','color', mycolor,'linewidth', 2)
    hold on
    semilogy([0,maxiter],cond(X)*[1,1],'-.',  'linewidth', 1,'color', [1,1,1]*(jmin/j))
    
    figure(2)
    % clf
    
    mycolor = [1,0,0]*(jmin/j);
    semilogy(cond_varp,  '-','color', mycolor,'linewidth', 2)
    hold on
    
    semilogy([0,maxiter], cond(X)*[1,1],'-.',  'linewidth', 1,'color', [1,1,1]*(jmin/j))
    
    
    figure(3)
    clf
    stem(beta_varpro)
end
%%

figure(1)
xlabel('Iterations', 'fontsize', 18)
% ylabel('Condition Number', 'fontsize', 18)
ax = gca;
ax.FontSize = 18;
xlim([1,8001])
ylim([0,1e20])

set(gca,'fontname','times')  % Set it to times



figure(2)
xlabel('Iterations', 'fontsize', 18)
% ylabel('Condition Number', 'fontsize', 18)
ax = gca;
ax.FontSize = 18;
xlim([1,8001])
ylim([0,1e20])

set(gca,'fontname','times')  % Set it to times


