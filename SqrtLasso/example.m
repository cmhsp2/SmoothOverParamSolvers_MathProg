%compare varpro, coordinate descent and scaled Lasso methods for solving
%  min_x norm(x,1) + norm(X*x - y)/(lambda*sqrt(m));
clearvars
addpath('../')
SaveResults = 0;
type = 'libsvm';
% type = 'Gaussian';
%
switch type
    case 'Gaussian'
        
        n = 2000;
        m = 300;
        X = randn(m,n);
        x0 = zeros(n,1);
        k = 40;
        x0(randperm(n,k)) = randn(k,1);
        y = X*x0+0.1*randn(m,1);
        name = '2000_300_40';
    case 'libsvm'

        name = 'mnist';
%         name = 'leukemia';
        datapath = sprintf('/Users/cmhsp20/Documents/MATLAB/varpro/Neurips2021/LassoDatasets/%s-1.mat',name);
        L = load(datapath);
        X = L.X;
        y = (L.Y)';
        [m,n] = size(X);
end
XtX = X'*X;
Xty = X'*y;

objbest = inf;

lambda0 = norm(X'*y, 'inf')/norm(y)/sqrt(m);
fac = 20;
lambda = lambda0/fac;
nmy = norm(y)^2;
if m<n
    func = @(x) norm(x,1) + norm(X*x - y)/lambda/sqrt(m);
else
    func = @(x) norm(x,1) + sqrt(x'*(XtX*x) + nmy - 2*Xty'*x)/lambda/sqrt(m);
end



OBJ = {};
NAME = {};
TIME = {};
COLOR = {};

%% Varpro for solving min_x norm(x,1) + norm(X*x - y)/lam;
tic
lam = lambda*sqrt(m);
[beta_varpro, fval_varpro] = func_SQLasso_varpro(X,y,lam);
time_varpro = toc
objbest = min(min(fval_varpro),objbest);

OBJ{end+1} = fval_varpro;
NAME{end+1} = 'VarPro';
TIME{end+1} = time_varpro;
COLOR{end+1} = 'b';


%% Coordinate descent
%solves min_beta sqrt( sum((X*beta-y).^2)/m )  +  (lam/m)*norm(beta,1);
lam = lambda*m;
MaxIter = 200;
tic
[betaSQ, sSQ,fvals] = SqrtLassoIterative_WebPage(X, y, lam, ones(n,1), MaxIter);
time_SQ = toc
fvals_SQ = fvals/lambda;
objbest = min(min(fvals_SQ), objbest);

OBJ{end+1} = fvals_SQ;
NAME{end+1} = 'CD';
TIME{end+1} = time_SQ;
COLOR{end+1} = 'r';




%% scaled lasso for solving  min_x norm(x,1) + norm(X*x - y)/lam;

tic
lam = lambda*sqrt(m);

beta = zeros(n,1);
fval_slasso  = func(beta);
time_slasso = 0;
for i=1:6
    beta_old = beta;
    fval_slasso(i) = func(beta);
    
    
    nm = max(1e-12,norm(X*beta - y))
    
    beta = func_l1_varpro(X,y,lam*nm, ones(n,1), XtX, Xty);
    time_slasso(i+1) = toc;
    fval_slasso(i+1) = func(beta);
    if norm(beta-beta_old)/norm(beta)<1e-8
        break
    end
    
end

objbest = min(min(fval_slasso),objbest);

OBJ{end+1} = fval_slasso;
NAME{end+1} = 'Scaled Lasso';
TIME{end+1} = time_slasso;
COLOR{end+1} = [0.25 0.80 0.54];

%%

display_objectives(OBJ, NAME, TIME, COLOR, 'VarPro','obj_min')

f = gcf;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])
if SaveResults
    exportgraphics(f,sprintf('results/%s_%s_%d.png',type,name,fac),'Resolution',300)
end
