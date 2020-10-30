function [sigma_r,sigma_t,sigma_rD,sigma_tD] = solve_stress_viscoelastic_shell(grid_r,mu,sigma_r_last,alpha_dTdotdr,Pex,E,nu,dt)
% Solver for viscoelastic stresses
% Subject to basal pressure boundary condition given by Pex
% surface boundary condition is sigma_r = 0 (free surface)
%
% VECTOR inputs
% grid_r contains nodal coordinates
% mu contains viscosities at each node
% sigma_r_last contains previous stress (sigma_r) values
% alpha_dTdotdr is the alpha_l * d/dt(dT/dr)
%
% SCALAR inputs
% Pex excess pressure
% E shear modulus
% nu poisson ratio
% dt timestep
%
%

nr = length(grid_r);
m1row = zeros(3*nr,1);
m1col = zeros(3*nr,1);
m1val = zeros(3*nr,1);
m2row = zeros(3*nr,1);
m2col = zeros(3*nr,1);
m2val = zeros(3*nr,1);
% M1 = sparse(nr,nr); % coefficients on (dsigma/dt)
% M2 = sparse(nr,nr); % coefficients on (sigma_r)
R = zeros(nr,1);
ind1=1;
ind2=1;
for i=1:nr
    if i==1
        drm = grid_r(i+1)-grid_r(i);
    else
        drm = grid_r(i) - grid_r(i-1);
    end
    if i==nr
        drp = grid_r(i)   - grid_r(i-1);
    else
        drp = grid_r(i+1) - grid_r(i);
    end
    rA = grid_r(i)+drp/2; % half a cell +
    rB = grid_r(i)-drm/2; % half a cell -
    drc = rA-rB;
    this_mu = mu( i ); % viscosity
    
    % M1 - coefficients of dsigma_r/dt
    const1 = (3-3*nu)/(2*E);
    const2 = (1-nu)/E; % coefficient on d(r/2 d/dr)/dr term
    coef_a = rA/2/drp/drc;
    coef_b = rB/2/drm/drc;
    coef_plus   =1/dt*(  const1/(drp+drm) + const2*coef_a);
    coef_center =1/dt*(                   - const2*coef_a - const2*coef_b);
    coef_minus  =1/dt*( -const1/(drp+drm) + const2*coef_b);
    if i==1
        m1row(ind1) = i; m1col(ind1) = i;   m1val(ind1) = coef_center;          ind1 = ind1+1;
        m1row(ind1) = i; m1col(ind1) = i+1; m1val(ind1) = coef_plus-coef_minus; ind1 = ind1+1;
        R(i) = R(i) - 2*coef_minus*Pex + 2*coef_minus*sigma_r_last(1); % sigma_r(i-1) + sigma_r(i+1) = 2*Pex
        %                 M1(i,i)   = coef_center;
        %                 M1(i,i+1) = coef_plus-coef_minus;
        %                 R(i) = R(i) - 2*coef_minus*P;
    elseif i==nr
        m1row(ind1) = i; m1col(ind1) = i-1; m1val(ind1) = coef_minus-coef_plus; ind1 = ind1 + 1;
        m1row(ind1) = i; m1col(ind1) = i;   m1val(ind1) = coef_center; ind1 = ind1 + 1;
        %         M1(i,i-1) = coef_minus-coef_plus;
        %         M1(i,i)   = coef_center;
        R(i) = R(i) - 2*coef_plus*0 + 2*coef_plus*sigma_r_last(nr);
    else
        m1row(ind1) = i; m1col(ind1) = i-1; m1val(ind1) = coef_minus;  ind1 = ind1+1;
        m1row(ind1) = i; m1col(ind1) = i;   m1val(ind1) = coef_center; ind1 = ind1+1;
        m1row(ind1) = i; m1col(ind1) = i+1; m1val(ind1) = coef_plus;  ind1 = ind1+1;
        %         M1(i,i-1) = coef_minus;
        %         M1(i,i)   = coef_center;
        %         M1(i,i+1) = coef_plus;
    end
    
    % M2 - coefficients of sigma_r
    if i == 1
        % extrapolate the log-viscosity gradient
        mu_ghost = exp( log(mu(i)) - ( log(mu(i+1))-log(mu(i)) ) );
        mu_B = exp(mean(log([mu_ghost mu(i)])));
    else
        mu_B = exp(mean(log([mu(i-1) mu(i)]))); % viscosity halfway between nodes i,i-1
    end
    if i == nr
        mu_ghost = exp( log(mu(i)) + ( log(mu(i))-log(mu(i-1)) ) );
        %         mu_ghost = exp(mean(log([mu(i-1) mu(i)])));
        mu_A = exp(mean(log([mu_ghost mu(i)])));
    else
        mu_A = exp(mean(log([mu(i) mu(i+1)]))); % viscosity halfway between nodes i,i+1
    end
    coef_plus   = 1/(4*this_mu)/(drp+drm) + rA/12/mu_A/drp/drc;
    coef_center =                          -rA/12/mu_A/drp/drc - rB/12/mu_B/drm/drc;
    coef_minus  =-1/(4*this_mu)/(drp+drm) + rB/12/mu_B/drm/drc;
    if i==1
        m2row(ind2) = i; m2col(ind2) = i;   m2val(ind2) = coef_center;          ind2=ind2+1;
        m2row(ind2) = i; m2col(ind2) = i+1; m2val(ind2) = coef_plus-coef_minus; ind2=ind2+1;
        %         M2(i,i)   = coef_center;
        %         M2(i,i+1) = coef_plus-coef_minus;
        R(i) = R(i) - 2*coef_minus*Pex;
    elseif i==nr
        m2row(ind2) = i; m2col(ind2) = i-1; m2val(ind2) = coef_minus-coef_plus; ind2=ind2+1;
        m2row(ind2) = i; m2col(ind2) = i;   m2val(ind2) = coef_center; ind2=ind2+1;
        %         M2(i,i-1) = coef_minus-coef_plus;
        %         M2(i,i)   = coef_center;
        R(i) = R(i) - 2*coef_plus*0; % surface sigma_r = 0
    else
        m2row(ind2) = i; m2col(ind2) = i-1; m2val(ind2) = coef_minus; ind2=ind2+1;
        m2row(ind2) = i; m2col(ind2) = i;  m2val(ind2) = coef_center; ind2=ind2+1;
        m2row(ind2) = i; m2col(ind2) = i+1; m2val(ind2) = coef_plus;  ind2=ind2+1;
        
        %         M2(i,i-1) = coef_minus;
        %         M2(i,i)   = coef_center;
        %         M2(i,i+1) = coef_plus;
    end
    if i==1
        %         R(i) = R(i) - coef_minus*2*Pex;
        %         m2row(ind2) = i; m2col(ind2) = i+1; m2val(ind2) = -coef_minus; ind2=ind2+1;
        %         M2(i,i+1) = M2(i,i+1) - coef_minus;
    elseif i==nr
%         m2row(ind2) = i; m2col(ind2) = i-1; m2val(ind2) = -coef_plus; ind2=ind2+1;
        %         M2(i,i-1) = M2(i,i-1) - coef_plus;
%         R(i) = R(i) - coef_plus*0; % no change because sigma_r = 0 at surface
    end
    R(i) = R(i) - alpha_dTdotdr(i);
    %         R(i) = R(i)+alpha_l*(Tdot(i+1)-Tdot(i))/2/drc; % this term
    %         includes the coupling to the energy equation - Tdot needs
    %         to be updated
end
m1row = m1row(1:ind1-1);
m1col = m1col(1:ind1-1);
m1val = m1val(1:ind1-1);
m2row = m2row(1:ind2-1);
m2col = m2col(1:ind2-1);
m2val = m2val(1:ind2-1);

M1 = sparse(m1row,m1col,m1val,nr,nr);
M2 = sparse(m2row,m2col,m2val,nr,nr);
LHS = (M1+M2);
R1term = M1*sigma_r_last; % this represents terms involving dsigma/dr at previous timestep
% R1term(1) = 0;
RHS = (R+R1term);

% LHS(1,:) = 0;
% LHS(1,1) = abs(LHS(2,2));
% RHS(1)   = Pex*LHS(1,1);
% LHS(nr,:) = 0;
% LHS(nr,nr) = abs(LHS(nr-1,nr-1));
% RHS(nr) = LHS(nr,nr)*0;
LHS = sparse(LHS);
sigma_r = LHS\RHS;

% 4. calculate the tangential stress sigma_t
% first, calculate dsr/dr
dsrdr = zeros(size(sigma_r));
for i=2:nr-1
    dr = grid_r(i+1)-grid_r(i-1);
    dsrdr(i) = (sigma_r(i+1)-sigma_r(i-1))/dr;
end
sigma_g = Pex - (sigma_r(2) - Pex);
dsrdr(1) =  (sigma_r(2)-sigma_g)/2/(grid_r(2)-grid_r(1)); % special formula using ghost value
sigma_g = 0 - (sigma_r(nr-1) - 0);
dsrdr(nr) = (sigma_g-sigma_r(nr-1))/2/(grid_r(nr)-grid_r(nr-1));

sigma_t = sigma_r+(grid_r'/2).*dsrdr;
sigma_t(end) = 0 + grid_r(end)/2 * dsrdr(end);

%deviatoric stresses, from Hillier and Squyres (1991) equations A8-9
sigma_tD =  grid_r'/6.*dsrdr;
sigma_rD = -grid_r'/3.*dsrdr;
end