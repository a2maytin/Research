% This script was used to make the figures from my presentation 
% on 7/20.
% The effect of the boundary on simulated trajectories was plotted.
% Andrew Maytin

close all
D = (0:0.01:2);
x = [0.1,0.3,0.5,0.7,0.9,1.5,2];
y = [.095, .275, .45, .61, .76, 1.2, 1.5];
y2 = x-y;

figure
plot(x,y,'.','MarkerSize',12)
%plot(x,y)
errorbar(x,y,.1*sqrt(y))
title('D_{model} vs. D_{sim}')
legend('Closest match 10 simulations','Location','northwest')
xlabel('Dsim (um^2/s)')
ylabel('Dmodel (um^2/s)')

figure
plot(x,x-y,'.','MarkerSize',12)
hold on
plot(x, .161*x.^(1.53))
errorbar(x,x-y,.1*sqrt(x-y),'LineStyle', 'none')
xlim([0,2.1])
ylim([0,0.55])
title('(D_{sim} - D_{model}) vs. D_{sim}')
legend('Closest match 10 simulations', 'Dsim-Dmodel = .161*Dsim.^{1.53}','Location','northwest')
xlabel('Dsim (um^2/s)')
ylabel('Dsim-Dmodel (um^2/s)')

Bp = polyfit(log10(x), log10(y2), 1);
Yp = polyval(Bp,log10(x));
figure
plot(log10(x), log10(y2), 'pg')
hold on
plot(log10(x), Yp, '-r')
hold off
grid
title('log(D_{sim} - D_{model}) vs. log(D_{sim})')
legend('Closest match 10 simulations', 'log(Dsim-Dmodel) = 1.53*log(Dsim)-.793','Location','northwest')
xlabel('log(Dsim)')
ylabel('log(Dsim-Dmodel)')
