% Gradient descent VarPro for TV denoising

n = 512; %% increasing this slow down the convergence

x = linspace(0,1,n)';
f = abs(x-.5)<=.3;

grad = @(x)x-x([2:end 1]);
div  = @(x)x([end 1:end-1])-x;
L = 4; % norm of laplacian



%% projected gradient descent on the dual
Proj = @(z)z ./ max(abs(z),1);
tau = 1.8/L;

lambda = 10;

mfun = @(x) lambda*norm(grad(x),1)+ norm(x-f)^2/2;

z = zeros(n,1);
niter = 10000;
E = [];
fval = [];

saveits = [1,500,1000,4000,9000];
for i=1:niter
    u = div(z) + f/lambda;
    E(i) = norm(u, 'fro')^2;
    z = Proj( z + tau*grad( u ) );
    % denoising result.
    f1 = f+lambda*div(z);
    
    
    %savefigure
    if any(ismember(saveits,i))
        clf; hold on;
        plot(x, f, 'k--');
        plot(x, f1, 'r', 'LineWidth', 2);
        axis([0 1 0 1]);
        axis off
        fghandle = gcf;
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4., 4.], 'PaperUnits', 'Inches', 'PaperSize', [4.25, 4.])

        exportgraphics(fghandle,sprintf('results/Chambolle_denoising_%d.png',i),'Resolution',300)
    end
    
    
    if mod(i,50)==1
        clf; hold on;
        plot(x, f, 'k--');
        plot(x, f1, 'r', 'LineWidth', 2);
        axis([0 1 0 1]);
        drawnow;
    end
    fval(i) = mfun(f1);
end



%% VarPro
% niter =100;
%sparse gradient matrix
e = ones(n,1);
D =  spdiags([ e -1*e ],0:1,n,n);
D(n,1)=-1;
% tau = 0.01;.1/L;
test = randn(n,1);
disp(['should be zero:', num2str(norm(D*test-grad(test)))])


v = ones(n,1);
fval_v = [];
saveits = [1,5,10,20,30];
for i=1:niter
    
    a = (lambda*(D*D')+spdiags(v.^2,0,n,n))\(D*f);
    g = v - a.^2.*v;
    v = v - tau*g;
    
    f1v = -lambda*D'*a + f;
    
    
    %savefigure
    if any(ismember(saveits,i))
        clf; hold on;
        plot(x, f, 'k--');
        plot(x, f1v, 'r', 'LineWidth', 2);
        axis([0 1 0 1]);
        axis off
        fghandle = gcf;
        set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4., 4.], 'PaperUnits', 'Inches', 'PaperSize', [16.25, 3.])

        exportgraphics(fghandle,sprintf('results/VarPro_denoising_%d.png',i),'Resolution',300)
    end

    if mod(i,50)==1
        clf; hold on;
        plot(x, f, 'k--');
        plot(x, f1v, 'r', 'LineWidth', 2);
        axis([0 1 0 1]);
        drawnow;
    end
   
    fval_v(i) = mfun(f1v);
end
a = (lambda*(D*D')+spdiags(v.^2,0,n,n))\(D*f);
f1v = -lambda*D'*a + f;

%% plot objective
figure(3)
clf
objmin = min(min(fval_v), min(fval));
loglog(fval - objmin)
hold on
loglog(fval_v - objmin, 'r')
legend('PCG', 'VarPro')

%%
C = diag(1-a.^2) + diag(a.*v)*inv(lambda*(D*D')+spdiags(v.^2,0,n,n)) *diag(a.*v);

