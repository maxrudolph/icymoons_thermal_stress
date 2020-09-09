function [T,dTdotdr] = solve_temperature_shell(grid_r,T_last,Tb,Ts,k,rho_i,Cp,H,dt,delta_rb)
% Solve the energy equation in spherical coordinates
% VECTOR inputs
% grid_r - node locations
% T_last - last timestep solution
%
% SCALAR inputs
% Tb - basal temperature
% Ts - surface temperature
% k - thermal conductivity
% rho - density
% Cp - heat capacity
% H - (constant) internal heating
% dt - timestep

nr = length(grid_r);

% L = sparse(nr,nr);
R = zeros(nr,1);
row=zeros(3*nr,1);
col=zeros(3*nr,1);
val=zeros(3*nr,1);
ind=1;
for i=1:nr
    r = grid_r(i);
    if i==1
        drm = grid_r(i+1)-grid_r(i);
    else
        drm = grid_r(i)-grid_r(i-1);
    end
    if i==nr
        drp = drm;
    else
        drp = grid_r(i+1)-grid_r(i);
    end
    rA = r + drp/2;
    rB = r - drm/2;
    kA = k;% thermal conductivities
    kB = k;
    dr = rA-rB;
    coef_plus   = -kA*rA^2/r^2/drp/dr;
    coef_center =  rho_i*Cp/dt + kA*rA^2/r^2/drp/dr + kB*rB^2/r^2/drm/dr;
    coef_minus  = -kB*rB^2/r^2/drm/dr;
    
    if( i==1 )
        row(ind) = i; col(ind) = i; val(ind) =  coef_center; ind = ind + 1;
        row(ind) = i; col(ind) = i+1; val(ind) =  coef_plus-coef_minus; ind = ind + 1;
        R(i) = rho_i*Cp/dt*T_last(i) + H - 2*Tb*coef_minus;
        %             L(i,i+1) = coef_plus-coef_minus;
        %             R(i) = R(i) - 2*Tb*coef_minus;
        %         R(i) = coef_center*Tb;
    elseif i==nr
        row(ind) = i; col(ind) = i;   val(ind) = coef_center; ind = ind+1;
        row(ind) = i; col(ind) = i-1; val(ind) = coef_minus-coef_plus; ind = ind+1;
        %             L(i,i-1) = coef_minus-coef_plus;
        R(i) = rho_i*Cp/dt*T_last(i) + H - 2*Ts*coef_plus;
        %         R(i) = coef_center*Ts;
    else
        row(ind) = i; col(ind) = i-1; val(ind) = coef_minus;  ind = ind + 1;
        row(ind) = i; col(ind) = i;   val(ind) = coef_center; ind = ind + 1;
        row(ind) = i; col(ind) = i+1; val(ind) = coef_plus;   ind = ind + 1;
        
        R(i) = rho_i*Cp/dt*T_last(i) + H;
    end
end
row = row(1:ind-1);
col = col(1:ind-1);
val = val(1:ind-1);
L = sparse(row,col,val,nr,nr);
T = L\R;

dTdr_b = (T(2)-Tb)/(grid_r(2)-grid_r(1));
Tdot = (T-T_last)/dt;
Tdot(1) = dTdr_b*delta_rb/dt;
dTdotdr = zeros(nr,1);
for i=2:nr-1
    dTdotdr(i) = (Tdot(i+1)-Tdot(i-1))/(grid_r(i+1)-grid_r(i-1));
end
dTdotdr(nr) = (0-Tdot(nr-1))/(grid_r(nr)-grid_r(nr-1));
%         dTdotdr(1)  = (Tdot(2)-Tdot(1))/(grid_r(2)-grid_r(1)); % First-order approximation at i=1
dr1 = grid_r(2)-grid_r(1);
dTdotdr(1) = 1/dr1*(2*(Tdot(2)-Tdot(1)) - (Tdot(3)-Tdot(1))/2); % 2nd order approximation to dTdr using a right-weighted stencil.
end