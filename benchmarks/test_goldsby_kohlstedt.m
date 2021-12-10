% test the goldsby kohlstedt rheology.
clear
close all

addpath ../core/

% T = linspace(100,273,1000);
T=273;
sigma = logspace(-4,2)*1e6;
d=100e-3;
P=1.3*1000*1e4;

for i=1:length(sigma)
    eta(i) = goldsby_kohlstedt(sigma(i),T,d,P);
end

figure()
plot(sigma/1e6,sigma./(2*eta));
xlabel('stress (MPa)');
ylabel('strain rate (1/s)');
set(gca,'YLim',[1e-12 1e0]);
set(gca,'XLim',[1e-4 1e2]);
set(gca,'XScale','log');
set(gca,'YScale','log');
% 
% hold on
% sigma = 1e7;
% for iT=1:length(T)
%     eta(iT) = goldsby_kohlstedt(sigma,T(iT),d,P);
% end
% plot(eta,T,'r');
