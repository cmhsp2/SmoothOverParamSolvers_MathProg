 %% plot results
m = linspace(80,113,12);
varpro1 = [0,1,3,6,8,11,16,31,42,48,64,76];
varpro1b = [0,3,4,8,12,17,33,44,67,80,83,93];
varpro2 = [3,5,13,25,34,60,76,91,93,96,99,100];
irls = [16,32,50,70,81,91,95,99,99,100,100,100];
rwl1 = [6,13,21,38,50,71,87,95,98,98,100,100];
l1 = [0,0,1,1,2,3,3,7,19,33,49,62];
clf


plot(m,varpro1, '-+', 'linewidth',2, 'Color', [0,0,1], 'Markersize', 15);
hold on
plot(m,varpro1b,'-^', 'linewidth',2, 'Color', [0,0,0.7], 'Markersize', 15);
plot(m,varpro2,'->', 'linewidth',2, 'Color', [0,0,0.5], 'Markersize', 15);
plot(m,irls,'-o', 'linewidth',2, 'Color', [1,0,0], 'Markersize', 15);
plot(m,rwl1,'-d', 'linewidth',2, 'Color', [0.75,0,0], 'Markersize', 15);


plot(m,l1, '-x', 'linewidth',2, 'Color', [0.5,0.5,0], 'Markersize', 15);

% lgd = legend('VarPro (1 inside)', 'VarPro v2  (1 inside)', 'VarPro (2 inside)', 'IRLS', 'Rw-l1','l1')
% lgd.FontSize = 24;
% lgd.Location = 'southeast';
xlabel('m','Fontsize', 24)
ylabel('# Successful Recovery','Fontsize', 24)

f = gcf;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])
exportgraphics(f,'results/T1.png','Resolution',300)

 

%%
figure


m = linspace(40,55,6);
varpro1 = [0,11,73,93,99,100];
varpro1b = [0,17,82,99,100,100];
varpro2 = [0,5,87,98,100,100];
irls = [0,5,80,96,100,100];
rwl1 = [0,0,85,99,99,100];
clf

% plot(m,varpro1, 'linewidth',2);
% hold on
% plot(m,varpro1b, 'linewidth',2);
% plot(m,varpro2, 'linewidth',2);
% plot(m,irls, 'linewidth',2);
% plot(m,rwl1, 'linewidth',2);

plot(m,varpro1, '-+', 'linewidth',2, 'Color', [0,0,1], 'Markersize', 15);
hold on
plot(m,varpro1b,'-^', 'linewidth',2, 'Color', [0,0,0.7], 'Markersize', 15);
plot(m,varpro2,'->', 'linewidth',2, 'Color', [0,0,0.5], 'Markersize', 15);
plot(m,irls,'-o', 'linewidth',2, 'Color', [1,0,0], 'Markersize', 15);
plot(m,rwl1,'-d', 'linewidth',2, 'Color', [0.75,0,0], 'Markersize', 15);

plot(m,zeros(size(m)), '-x', 'linewidth',2, 'Color', [0.5,0.5,0], 'Markersize', 15);


lgd = legend('VarPro (1 inside)', 'VarPro v2  (1 inside)', 'VarPro (2 inside)', 'IRLS', 'Rw-l1','l1')
lgd.FontSize = 24;
lgd.Location = 'southeast';
% lgd = legend('VarPro (1 inside)', 'VarPro v2  (1 inside)', 'VarPro (2 inside)', 'IRLS', 'Rw-l1')
xlabel('m','Fontsize', 24)
ylabel('# Successful Recovery','Fontsize', 24)
lgd.FontSize = 24;
lgd.Location = 'southeast';
f = gcf;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])

exportgraphics(f,'results/T50.png','Resolution',300)


%%
figure


m = linspace(40,55,6);
varpro1 = [0,39,94,94,99,100];
varpro1b = [0,46,97,99,99,100];
varpro2 = [0,10,96,99,99,100];
irls = [0,22,85,96,99,100];
rwl1 = [0,7,90,99,99,100];
clf

% plot(m,varpro1 , 'linewidth',2);
% hold on
% plot(m,varpro1b, 'linewidth',2);
% plot(m,varpro2 , 'linewidth',2);
% plot(m,irls, 'linewidth',2);
% plot(m,rwl1, 'linewidth',2);


plot(m,varpro1, '-+', 'linewidth',2, 'Color', [0,0,1], 'Markersize', 15);
hold on
plot(m,varpro1b,'-^', 'linewidth',2, 'Color', [0,0,0.7], 'Markersize', 15);
plot(m,varpro2,'->', 'linewidth',2, 'Color', [0,0,0.5], 'Markersize', 15);
plot(m,irls,'-o', 'linewidth',2, 'Color', [1,0,0], 'Markersize', 15);
plot(m,rwl1,'-d', 'linewidth',2, 'Color', [0.75,0,0], 'Markersize', 15);

plot(m,zeros(size(m)), '-x', 'linewidth',2, 'Color', [0.5,0.5,0], 'Markersize', 15);



lgd = legend('VarPro (1 inside)', 'VarPro v2  (1 inside)', 'VarPro (2 inside)', 'IRLS', 'Rw-l1','l1')
lgd.FontSize = 24;
lgd.Location = 'southeast';
% 
% lgd = legend('VarPro (1 inside)', 'VarPro v2  (1 inside)', 'VarPro (2 inside)', 'IRLS', 'Rw-l1')
% % lgd.FontSize = 24;
% lgd.Location = 'southeast';
xlabel('m','Fontsize', 24)
ylabel('# Successful Recovery','Fontsize', 24)
f = gcf;
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 9., 7.], 'PaperUnits', 'Inches', 'PaperSize', [9., 7.])

exportgraphics(f,'results/T100.png','Resolution',300)
