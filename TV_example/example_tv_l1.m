%TV regularisation with L1 loss. 


clearvars
%2D gradient matrix
addpath('../')

SaveImages = 0;
OBJECTIVE = {};
NAME = {};
TIME = {};
COLOR = {};
type = 'colour';
% type = 'phantom';
% type = 'HSData';
switch type
    case 'phantom'
        imgname = 'phantom';
        n1 = 64;
        n2 = 64;
        x0 = phantom(randn(8,6),n1);
        x0 = phantom(n1);
        x0 = x0(:);
        p=n1*n2;
        lambda = .7;
        mode = 'proxcalc';
        
    case 'colour'
        imgname = 'peppers';
        imgname = 'pears';
%         imgname = 'hestain';
        im = im2double(imread(sprintf('%s.png',imgname)));
        [n1,n2,~] = size(im);
        p = n1*n2;
        x0 = reshape(im,n1*n2,[]);
        lambda = .6;
        mode = 'proxcalc';
    case 'HSData'
%         http://www.ehu.eus/ccwintco/index.php/Hyperspectral_Remote_Sensing_Scenes
        imgname = 'HSData';
        L = load('Indian_pines.mat');
        im = L.indian_pines;
        im = im/max(abs(im(:)));
        [n1,n2,~] = size(im);
        p = n1*n2;
        x0 = reshape(im,p,[]);
        lambda = .5;
        mode = 'proxcalc';
        
end

%%
%gradient matrix
[T1,T2] = SparseGrad(n1,n2);
D = [T1;T2];

isotropic= 1;


% mode = 'proxcalc';
% mode = 'inpaint';
% mode = 'sparse';
switch mode
    case 'proxcalc' %denoising, A is identity
        A = speye(p);
        m = size(A,1);
        proxcalc = 1;
    case 'inpaint' %inpainting, A is masking operator
        p = n1*n2;
        A = speye(p);
        indx = randperm(p,ceil(0.3*p));
        A = A(indx, :);
        m = size(A,1);
        proxcalc = 0;
    otherwise %A is random sparse matrix
        proxcalc = 0;
        m = ceil(p);
        k=ceil( m/5);
        
        irand = randperm(m,k);
        jrand = randperm(p,k);
        zrand = sign(randn(k,1));
        A = sparse(irand,jrand,zrand, m,p);
        
        %     m = ceil(p/4);
        %     X = randn(m,p)/sqrt(m);
        
end

%%
y = A*x0;
y = imnoise(y,'salt & pepper',0.25);
[m,n] = size(A);
%%

if isotropic
    t0 = randn(p+m,1);
else
    t0 = randn(2*p+m,1);
end
objbest = inf

%%
Aeye =strcmp(mode, 'proxcalc');
tic
[x_pro,R] = func_TV_L1_pro(A,y,lambda,D,t0,isotropic, Aeye);
time1 = toc
%%
figure(2)
imagesc(reshape(x_pro,n1,[]));
% imagesc(reshape(beta_pro,n1,n2,[]));
f_pro = R.err(:,1);
OBJECTIVE{end+1} = f_pro;
NAME{end+1} = 'VarPro';
TIME{end+1} = time1;
 COLOR{end+1} = 'b';
objbest = min(objbest, min(f_pro));


%% Primal-Dual
[m,n] = size(A);
    p = size(D,1);
    d = size(y,2);
sigvals = [0.1,1,10];
% sigvals = [ 1];
K = [D -speye(p) sparse(p,m); A sparse(m,p) -speye(m) ];

    
for i=1:length(sigvals)
    sig = sigvals(i);
    tau = 0.9/normest(K)^2/sig;
    niter = 1000;
    Aty = A'*y;
    
    
    S2 = @(x) sum(sqrt(sum(x(1:end/2,:).^2+x(end/2+1:end,:).^2,2)));
    mfunc = @(x) S2(D*x(1:n,:)) + sum(sqrt(sum((A*x(1:n,:) - y).^2, 2)))/lambda;
    proxG = @(x) [x(1:n,:); groupthresh_TV(x(1+n:n+p,:), lambda*tau); groupthresh_2(x(n+p+1:end,:), tau)];
    proxFs = @(xi) [xi(1:p,:); xi(p+1:end,:) - sig*y ];

    record_interval = 10;
    
    init.x = [Aty; randn(p+m,d)];
    init.xi = randn(p+m,d);
    %%
    
    
    tic
    [x,fvals] = func_PrimalDual(K,proxG,proxFs, sig,tau,niter,mfunc,record_interval,init);
    time2 = toc
    
    
    OBJECTIVE{i+1} = fvals;
    NAME{i+1} = sprintf('PD \\sigma=%.2g',sig);
    TIME{i+1} = time2;
    COLOR{i+1} = [1,0,0]*i/length(sigvals);
    objbest = min(objbest,min(fvals));
    
end

%%
display_objectives(OBJECTIVE, NAME, TIME, COLOR, 'VarPro','obj_min')
% 
% xlim([0,180])
ylim([1e-4,inf])

f = gcf;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])

if SaveImages
    exportgraphics(f,sprintf('results/%s/tvl1%s_%d_%d_%.2f.png',imgname,mode, n1,n2,lambda),'Resolution',300)
end
%%
%
% ylim([1e-4,inf])
% f = gcf;
% set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])
% exportgraphics(f,sprintf('results/%s/tvl1%s_%d_%d_%.2f.png',imgname,mode, n1,n2,lambda),'Resolution',300)


if SaveImages
    
    if strcmp(type,'colour')
        im_pro =reshape(x_pro,n1,n2,[]);
        yim =reshape(y,n1,n2,[]);
        imwrite(im_pro, sprintf('results/%s/tvl1Reconstr_%s_%d_%d_%.2f.png', imgname,mode,n1,n2,lambda));
        imwrite(yim, sprintf('results/%s/tvl1Input_%s_%d_%d_%.2f.png', imgname,mode,n1,n2,lambda));
        
    end
end
