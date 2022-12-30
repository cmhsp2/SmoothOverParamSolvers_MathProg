n= 3000;
m = 300;


SaveResults = 0;
addpath('../')

%define I, where I{j} are the indices in the jth group.
I = {};
x0 = zeros(n,1);
s = 1;
over = 5; %overlap between the groups
while s<n
    p = 10;%randi(2)-1; %size of group
    u = min(s:s+p+over-1,n);%sort(randi(n,[10,1]));
    s = s+p;
    I{end+1} = unique(u);
end

%define group truth, with k nonzero groups
k = 3;
for i=1:k
    u = randi(length(I));
    x0(I{u}) = randn(length(I{u}),1);
    
end

%  data matrix
X = randn(m,n)/sqrt(m);
y = X*x0;


% lambda = norm(X'*y, 'inf')/10;
% define the operator D, so that Dx = (y_i), 
% where (y_i)_{i=idx(j), ... idx(j+1)-1} are the elements in the jth group of
% x.
E = speye(n); D = [];
idx = 1;
p=1;
for i=1:length(I)
    D = [D;sqrt(length(I{i}))* E(I{i},:)];
    p = p+length(I{i});
    idx = [idx,p]; %store indices i of where (Dx)_i starts a new group
end
 
proxcalc = 0;% set to 0 if X is not the identity

%%
fac = 20;
lambda_max = 0;
z = X'*y;

for i=1:length(I)
    lambda_max = max(lambda_max,norm(z(I{i}))/sqrt(length(I{i})));
end
lambda = lambda_max/fac;
%%
OBJECTIVE = {};
NAME = {};
TIME = {};
COLOR = {};

% VarPro
[beta_pro,R] = func_OverlappingGroupLasso_varpro(D,X,y,lambda, idx,proxcalc);

warning on;

time = toc
fval_pro = R.err(:,1);

OBJECTIVE{end+1} = fval_pro;
NAME{end+1} = sprintf('VarPro');
TIME{end+1} = time;
COLOR{end+1} = 'b';

clf
stem(x0)
hold on
stem(beta_pro)
objbest = min(fval_pro)

%%

A = X;
niter = 1000;
gamvals = [0.05,1,5,10];
mfunc = @(x) Reg(x,I) + norm(A*x - y)^2/2/lambda;

for i = 1:length(gamvals)
    gamma = gamvals(i);
    proxG = @(z) mygroupthresh(z, lambda/gamma,idx);
    tic
    [x,fvals] = func_ADMM_overlap_group(proxG,D,A,y,gamma,niter,mfunc);
    timedr = toc
    OBJECTIVE{end+1} = fvals;
    NAME{end+1} = sprintf('ADMM \\gamma=%.2g',gamma);
    TIME{end+1} = timedr;
    COLOR{end+1} = [1  0 0]*i/length(gamvals) ;
    
    objbest = min(objbest,min(fvals));

end

%%
display_objectives(OBJECTIVE, NAME, TIME, COLOR, 'VarPro','obj_min')

ylim([1e-10,2e3])

f = gcf;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])
if SaveResults
    exportgraphics(f,sprintf('results/%d_%d_%d_%d.png',over, n,m,fac),'Resolution',300)
end
