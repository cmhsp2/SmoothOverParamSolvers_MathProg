clearvars

SaveResults = 0;
addpath('../')

load('data_cancer.mat')
group_newindex = zeros(size(group));

for i = 1:500
    group_newindex(group == gene(i)) = i;
end

%% Visualize connection of 500 genes
group_newindex = sort(group_newindex, 2);
M = zeros(500);
idx = sub2ind(size(M), group_newindex(:, 1), group_newindex(:, 2));
M(idx) = 1;
M_viso = M + M';
spy(M_viso);
title('correlations of genes')

%% Find not connected genes
range_gene = [1:500];
idx = ismember(range_gene, unique(sort(group_newindex(:))));
node_nconct = range_gene(~idx);

%% Generate group cells

% sigle nodes
G1 = mat2cell(range_gene(~idx)', ones(sum(~idx), 1));
% connected nodes
G2 = mat2cell(group_newindex, ones(size(group_newindex, 1), 1));
% total group
G_idx = [G1; G2];


I = G_idx;


addintercept = 0;
[m,n] = size(expression');


if addintercept
    addones = @(X) [X, ones(size(X,1),1)];
    E = speye(n+1);
    I{end+1} = n+1;
else
    addones = @(X) X;
    E = speye(n);
end


X = addones(expression');
y = labels;


% define the overlapping group operator D, so that Dx = (sqrt(n_j)*y_i),
% where (y_i)_{i=idx(j), ... idx(j+1)-1} are the elements in the jth group of
% x and n_i are the number of elements in group j.
 D = [];
idx = 1;
p=1;
for i=1:length(I)
    D = [D;sqrt(length(I{i}))* E(I{i},:)];
    p = p+length(I{i});
    idx = [idx,p]; %store indices i of where (Dx)_i starts a new group
end
S2 = @(x) Sum2(x,idx);
vec = @(x) x(:);
DDt = D*D';



%% run on full data matrix
lambda_max = 0;
z = X'*y;
fac = 5
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

tic
[beta_pro,R] = func_OverlappingGroupLasso_varpro(D,X,y,lambda, idx,0);
R = R.err(:,1);
time = toc
OBJECTIVE{end+1} = R;
NAME{end+1} = sprintf('VarPro');
TIME{end+1} = time;
COLOR{end+1} = 'b';

clf
stem(beta_pro)



%%

A = X;
niter = 1000;
% gamvals = [.5,1,2,5,10,20];
gamvals = [0.5,1,10,20];
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
end

%%

display_objectives(OBJECTIVE, NAME, TIME, COLOR, 'VarPro','obj_min')

f = gcf;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])
if SaveResults
    exportgraphics(f,sprintf('results/cancer/%d_%d_%d.png', n,m,fac),'Resolution',300)
end
%%