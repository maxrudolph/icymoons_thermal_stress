function [T_last,sigma_r_last,sigma_t_last,er_last,et_last] = interpolate_solution(new_grid_r,grid_r,T_last,sigma_r_last,sigma_t_last,er_last,et_last,Tb)

interp_r    = grid_r;
interp_r(1) = new_grid_r(1);% account for thickening by extending first cell.

tmp = T_last;
tmp(1) = Tb;
T_last = interp1(interp_r,tmp,new_grid_r)';
sigma_r_last = interp1(interp_r,sigma_r_last,new_grid_r)';% note imposes sigma_r=sigma_t at base
sigma_t_last = interp1(interp_r,sigma_t_last,new_grid_r)';
er_last = interp1(interp_r,er_last,new_grid_r)';
et_last = interp1(interp_r,et_last,new_grid_r)';

% re-mesh onto new grid
% new_grid_r = linspace(Ri-z,Ro,nr);
% interp_r = [new_grid_r(1) grid_r];
% T_last = interp1(interp_r,[Tb; T_last],new_grid_r)'; % very important to set this equal to Tb to prevent drift due to numerical inaccuracy.
% sigma_r_last = interp1(interp_r,[sigma_r_last(1); sigma_r_last],new_grid_r)';
% sigma_t_last = interp1(interp_r,[sigma_t_last(1); sigma_t_last],new_grid_r)';
% er_last = interp1(interp_r,[er_last(1); er_last],new_grid_r)';
% et_last = interp1(interp_r,[et_last(1); et_last],new_grid_r)';
