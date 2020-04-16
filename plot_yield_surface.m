clear
close all

ns = 200;
smin = -3e6;
smax = 3e6;
sigma = linspace(smin,smax,ns);

f = zeros(ns,ns);
for i=1:ns
   for j=1:ns
        sig_r = sigma(i);
        sig_t = sigma(j);
        smean = 1/3*(sig_r + 2*sig_t);
        sigma_rD = sig_r - smean;
        sigma_tD = sig_t - smean;
        J2 = 1/2*(sigma_rD^2 + 2*sigma_tD^2);
        sig_eff = sqrt(3*J2);
        
        p(i,j) = -smean;
        f(i,j) = sig_eff - (1.5e5 + 0.1*p(i,j));
   end
end


figure()
contourf(sigma,sigma,f,256,'color','none'); colorbar;
hold on
contour(sigma,sigma,f,[0 0],'k');
axis equal;

figure()
pcolor(sigma,sigma,p); shading flat; colorbar;
axis equal;