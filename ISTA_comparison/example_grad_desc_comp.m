%% compare gradient descent in the VarPro, Hadamard formulationw with ISTA
clearvars

resultsdir = 'results/';
mkdir(resultsdir)

rng(1124,'twister')


%fourier
n=  50;
fc = 2;
t = linspace(0,1,n);
Xm =  exp(-2*pi*1i*(-fc:fc)'*t)/sqrt(2*fc+1);
X0 = [real(Xm); imag(Xm)];
m = size(X0,1);
X = X0;
L = norm(X)^2;

beta0 = zeros(n,1);
% beta0(randperm(n,2)) = [.5;-.5];
supp0 = [n/2,ceil(n/4)];
supp0(:);
beta0(supp0) = [1;-1];
y = X*beta0;

lambda = norm(X'*y,'inf')/1000;

mfun = @(beta) 0.5*norm(X*beta - y)^2+lambda*norm(beta,1);
SoftThr = @(beta,tau)sign(beta) .* max(abs(beta)-tau,0);



%%
stepsize = 'fixed';
% stepsize = 'BB';
tol_bb = 1e-20;

niter = 5e5;
niter_ista = niter;
niter_hada = niter;
niter_varp = niter;
niter_hype = niter;

dispint = 1e5;%display interval
%% ground truth

a = (X(:,supp0)'*X(:,supp0))\(X(:,supp0)'*y -lambda*sign(beta0(supp0)));
beta_star = zeros(n,1);
beta_star(supp0)=a;
obj_best = mfun(beta_star);

%% ISTA 
niter = niter_ista;
clf
beta = zeros(n,1);
tau0 = 1/L;
ISTA = [];
b_m = beta;
g =  X'*(X*beta-y);
for it=1:niter
    b_m2 = b_m;
    g_m = g;
    b_m = beta;
    ISTA(end+1) = mfun(beta);
    g = X'*(X*beta-y);
    gnorm = norm(g_m-g)^2;
    
    if strcmp(stepsize,'BB') && it>1% gnorm>tol_bb
        tau = norm((b_m2 - beta))^2/sum((g_m - g).*(b_m2 - beta));
%         tau = sum((g_m - g).*(b_m2 - beta))/gnorm;
        tau = max(0.1/L, tau);
        tau = min(tau,10/L);
    else
        tau = tau0;
    end
    beta = SoftThr( beta - tau * g, tau*lambda );
    
    if mod(it,dispint)==0
        area(linspace(0,1,n),beta, 'FaceColor', 'b');
        title('ista')
        drawnow
    end
    
end
beta_ista = beta;

%% nonconvex hadamard
clf
niter= niter_hada;
niter= 2000;

HADA = [];

gradv = @(uv) uv(:,2) + 1/lambda * uv(:,1) .* X'*( X*(prod(uv,2)) - y );
gradu = @(uv) uv(:,1) + 1/lambda * uv(:,2) .* X'*( X*(prod(uv,2)) - y );
grad = @(uv) [gradu(uv), gradv(uv)  ];
fun = @(uv) mfun(prod(uv,2));

HADA_uvgrad = [];


uv = randn(n,2)/sqrt(n);
% uv(:,1) = uv(:,2)-0.002;

uv_m = uv;
g = grad(uv);
sumv = @(v) sum(v(:));
uv_store = uv(:);
for it=1:niter
    uv_m2 = uv_m;
    g_m = g;
    uv_m = uv;
    
    HADA(end+1) = mfun(prod(uv,2));
    HADA_uvgrad(end+1) = norm(g(:));
    g = grad(uv);
    
    
    gnorm = norm(g_m-g, 'fro')^2;
    if strcmp(stepsize,'BB') && gnorm>tol_bb && it>1
        tau = (sumv((g_m - g).*(uv_m2 - uv)))/gnorm;
        tau = min(max(0.01,tau),5);
    else
        tau = 0.01;
    end
    
    uv = uv - tau*g;
    if mod(it,dispint)==0
        area(linspace(0,1,n),prod(uv,2), 'FaceColor', 'b');
        title('hadamard')
        drawnow
    end
    uv_store(:,end+1) = uv(:);
end


%% % VarPro Flow

clf

niter=niter_varp;
VARP = [];
Matv = @(v) lambda*speye(m)+X*spdiags(v.^2,0,length(v),length(v))*X';
v0 = ones(n,1)/sqrt(n);
v = v0;
v_m = v;
Xta = -X'*(Matv(v)\y);
g = v - v.*(Xta.^2);
for it=1:niter
    v_m2 = v_m;
    g_m = g;
    v_m = v;
    
    alph = -Matv(v)\y;
    Xta = X'*alph;
    beta = -(v.^2).*Xta;
    VARP(end+1) = mfun(beta);
    g = v - v.*(Xta.^2);
        
    gnorm = norm(g_m-g)^2;
    if norm(g)<1e-16
        break
    end
    if gnorm>tol_bb && strcmp(stepsize,'BB') && it>1
        tau = (sum((g_m - g).*(v_m2 - v)))/gnorm;
        tau = min(max(tau,0.1),20);
    else
        tau = .1;
    end
    
    v = v - tau*g;
    if mod(it,dispint)==0
        clf
        area(linspace(0,1,n),beta, 'FaceColor', 'b');
        title('ncvx-pro')
        drawnow
    end
    u = -v.*Xta;
    
    if norm(v-v_m)<1e-8
        break
    end
    
end

beta_pro = beta;

%% hyperbolic entropy

% stepsize = 'fixed';
niter = niter_hype;

% hyperbolic entropy
clf

beta = zeros(n,1);
tau0 = 1/L;
HYPE = [];
b_m = beta;
g =  X'*(X*beta-y);
m_obj = inf;
param =  1/n; %parameter for hyperbolic entropy


for it=1:niter
    b_m2 = b_m;
    g_m = g;
    b_m = beta;
    HYPE(end+1) = mfun(beta);
    g = X'*(X*beta-y);
    
    tau = .1*tau0/param;
    

    beta = param*sinh(SoftThr( asinh(beta/param) - tau * g, tau*lambda ));
    
    if mod(it,dispint)==0
        area(linspace(0,1,n),beta, 'FaceColor', 'b');
        title('hyperbolic')
        drawnow
    end
    
end
beta_hype = beta;



%%
figure(1)

clf
loglog(ISTA-obj_best,'Color',[1 0 1], 'linewidth',3)
hold on
loglog(VARP-obj_best,'b', 'linewidth',3)

loglog(HADA-obj_best,'Color',[0 .8 1],'linewidth',3)

loglog(HYPE-obj_best,'Color',[1 0 0],'linewidth',3)

loglog(.05./(1:niter).^(2/3),'--','Color',[0.3 0 0],'linewidth',1)

loglog(3./(1:niter),'--','Color',[0.6 0.6 0.6],'linewidth',1)


lh = loglog(HADA_uvgrad, '-.','Color',[0,0,0],'linewidth', 4);
lh.Color = [lh.Color 0.4];

lg = legend('ISTA','VarPro',  'Hadamard', 'Hyperbolic', 'k^{-2/3}', 'k^{-1}', 'Hadamard \nabla','NumColumns',2, 'location', 'southwest')

legend boxoff 
ylim([1e-7, 2e1])

% ylim([1e-8, 1e1])

lg.FontSize = 24;
xlabel('k','fontsize',24)
ylabel('Objective value error','fontsize',24)
%
f = gcf;
set(gca,'FontSize',24)


set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 9.], 'PaperUnits', 'Inches', 'PaperSize', [9., 9.])
name = [resultsdir, '/fullcomp_gd_%d_%s.png'];
exportgraphics(f,sprintf(name, n, stepsize),'Resolution',300)



%%
figure(2)
clf
loglog(ISTA-obj_best,'Color',[1 0 1], 'linewidth',2)
hold on
loglog(VARP-obj_best,'b', 'linewidth',4)


loglog(HADA-obj_best,'Color',[0 .8 1],'linewidth',2)

loglog(.05./(1:niter).^(2/3),'--','Color',[0.3 0 0],'linewidth',1)

loglog(3./(1:niter),'--','Color',[0.6 0.6 0.6],'linewidth',1)

lg = legend('ISTA','VarPro',  'Hadamard',  'k^{-2/3}', 'k^{-1}', ...
    'location', 'northeast','NumColumns',2)

legend boxoff 

% ylim([1e-6, 5e1])

ylim([1e-7, 2e1])

lg.FontSize = 24;
xlabel('k','fontsize',24)
ylabel('Objective value error','fontsize',24)

f = gcf;
set(gca,'FontSize',24)


set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 9.], 'PaperUnits', 'Inches', 'PaperSize', [9., 9.])
name = [resultsdir, '/comp_gd_%d_%s.png'];
exportgraphics(f,sprintf(name, n, stepsize),'Resolution',300)

