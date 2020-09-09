clear
close all

nx = 100;
H=1;
Tl=0;
Tr=1;
rho=1000;
cp=5;
kappa=1e-6;
k=kappa*rho*cp;

x = linspace(0,1,nx);

L = zeros(nx,nx);
R = zeros(nx,1);
for i=1:nx
    dx=x(2)-x(1);
    coef_minus = k/dx/dx;
    coef_center = -k/dx/dx - k/dx/dx;
    coef_plus = k/dx/dx;
    if i==1
        L(i,i) = coef_center;
        L(i,i+1) = coef_plus - coef_minus;
        R(i) = H -2*coef_minus*Tl;
    elseif i==nx
        L(i,i-1) = coef_minus - coef_plus;
        L(i,i) = coef_center;
        R(i) = H -2*coef_plus*Tr;
    else
        L(i,i-1) = coef_minus;
        L(i,i) = coef_center;
        L(i,i+1) = coef_plus;
        R(i) = H;
    end
end

T = L\R;
figure, plot(x,T);

dtdx = zeros(nx,1);
for i=2:nx-1
    dtdx(i) = (T(i+1)-T(i-1))/dx/2;
end
dtdx(1) = (T(2)-T(1))/dx;
dtdx(nx) = (T(nx)-T(nx-1))/dx;
figure, plot(x,dtdx)

% 2nd method - note that even though T is not correct at the boundary, the
% approximation dt/dx is consistent.
dtdx = zeros(nx,1);
for i=2:nx-1
    dtdx(i) = (T(i+1)-T(i-1))/dx/2;
end
Tg = Tl - (T(2)-Tl);
dtdx(1) = (T(2)-Tg)/2/dx;
dtdx(nx) = (T(nx)-T(nx-1))/dx;
figure, plot(x,dtdx)
