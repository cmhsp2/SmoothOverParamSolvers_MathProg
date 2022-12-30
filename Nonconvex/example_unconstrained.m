clc
clearvars
%Figure 11 in paper. Uncontrained Lq minimisation.
clf

for runs = 1:20
m =45;
% m =85;
n=256;
T = 100;
s = 40;

TIME= {};
OBJ = {};
NAME = {};


q= 2/3;%should be at least 0.5
lambda = 0.1;
% normq = @(x,q) sum(sqrt(sum(abs(x).^2,2)).^q);


X = randn(m,n);
x0 = randn(n,T);
x0(randperm(n,n-s),:)=0;
y = X*x0;

myfunc = @(x) 3/2*sum(sqrt(sum(abs(x).^2,2)).^q) + 1/(2*lambda)*norm(X*x-y,'fro')^2;


%varpro lq (2 factors on outside, linear system on inside)
disp(['run: ',num2str(i)]);
disp('Running varpro with 3 factors, 1 inside ...');


tic
[b,f] = func_lq_varpro(X,y,lambda,q);
toc
err = norm(b-x0)/norm(x0)

t=toc;

objbest = min(f);

OBJ{end+1} = f;
NAME{end+1} = 'Varpro (1 inside)';
TIME{end+1} = t;

%%
disp('Running varpro with 3 factors, 2 inside (solved by l1-varpro) ...');

tic
[x_u, f_u] = func_lq_nonlin_varpro(X,y,lambda, q);
t=toc;

err_u = norm(x_u-x0)/norm(x0)
% myfunc(x_u)

objbest = min(objbest,min(f_u));


OBJ{end+1} = f_u;
NAME{end+1}= 'Varpro (2 inside)';
TIME{end+1} = t

%% IRLS
lambda0 = lambda*3/2;
disp('Running IRLS ...')
tic
[irls_x, f_irls ]= IRLS(X,y,lambda0,q);
t=toc;

err_irls = norm(irls_x-x0)/norm(x0)

f_irls = f_irls/lambda;

objbest = min(objbest,min(f_irls));


OBJ{end+1} = f_irls;
NAME{end+1} = 'IRLS';

TIME{end+1} = t;
%% reweighted l1
disp('runnning Reweighted l1 ...')
tic
[x_rwl1, f_rwl1] = rwl1(X,y,lambda0,q);
t=toc;
err_rwl1 = norm(x_rwl1-x0)/norm(x0)


f_rwl1 = f_rwl1/lambda;

objbest = min(objbest,min(f_rwl1));


OBJ{end+1} = f_rwl1;
NAME{end+1} = 'RWl1';
TIME{end+1} = t;


%%
newcolors = [0.83 0.14 0.14 
             1.00 0.54 0.00 
             0.47 0.25 0.80 
             0.25 0.80 0.54 ];
         
colororder(newcolors)

name = [];
for i = 1:length(TIME)
    
    time =  TIME{i};
    f = OBJ{i};
    
    
    lh = semilogy(linspace(0,time, length(f)),(f-objbest)./f, 'linewidth',2);
    lh.Color = [lh.Color 0.6];

%     name(end+1) = NAME{i};
    hold on
end
% legend(NAME)
if runs==1
    lgd = legend(NAME, 'AutoUpdate','off');
    lgd.FontSize =24;
end
    
end
ylabel('Objective value error', 'fontsize',24);

xlabel('Time in seconds', 'fontsize',24);
%%
ylim([1e-4,inf])
xlim([0,2])
f = gcf;
exportgraphics(f,sprintf('results/T%d.png',T),'Resolution',300)