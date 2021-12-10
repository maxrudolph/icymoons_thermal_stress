% test the goldsby kohlstedt rheology.
% I used this code to reproduce the curves in Figure 7
% of Goldsby and Kohlstedt (2001).
clear
close all

addpath ../core/

sigma = logspace(-4,2)*1e6;
T = 258* ones(size(sigma));
d=1e-3;
P=1.3*1000*1e4;

eta = goldsby_kohlstedt(sigma,T,d,P);

figure()
plot(sigma/1e6,sigma./(2*eta));
xlabel('stress (MPa)');
ylabel('strain rate (1/s)');
set(gca,'YLim',[1e-12 1e0]);
set(gca,'XLim',[1e-4 1e2]);
set(gca,'XScale','log');
set(gca,'YScale','log');

