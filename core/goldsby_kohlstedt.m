function [etaeff] = goldsby_kholstedt(s,T,d,P)
% s is stress in Pascal
% T is temperature in K
% d is grain size in m
% P is pressure in Pascal
% etaeff is effective viscosity in Pa-s
assert(all(size(s) == size(T))); % ensure that sizes of s and T are compatible.
%
%
% PetscScalar goldsbyKohlstedt(PetscScalar s, PetscScalar T, PetscScalar d, PetscScalar P){
%   /* returns strain rate predicted by goldsby and kohlstedt 2001 in units of s^-1*/
R = 8.3144621; % ideal gas constant, J/mol/K */
%   /* diffusion creep */

if( P < 1.0 )
    P = 1.0; %/* prevent spurious pressures from affecting viscosity */
end
%   {
%     /* quantities from table 6: Diffusion Creep Parameters */
b = 4.52e-10; %/* burgers vector */
Vm = 1.97e-5; %/* m^3 */
Qv = 59.4*1000.0;%/* Activation energy for volume diffusion, J/mol */
delta = 2.0*b; %/* grain boundary width, m */
Qb = 49*1000.0;
D0v = 9.10e-4; %/* m^2/s */
D0b = 5.8e-4; %/* m^2/s - this is lower of two upper bounde presented in text */
%     /* diffusivities */
Dv = D0v*exp(-Qv./(R*T));
Db = D0b*exp(-Qb./(R*T));
%     /* note that all units above are SI, so strain rate will be in 1/s for stress in Pa */
ediff0 = 42.0*Vm./(R*T*d*d).*(Dv + pi*delta/d*Db );

ediff = ediff0 .* s;
%   {
%     PetscScalar A,n,Q;
egbs = zeros(size(T));
V = -13e-6; %/* m^3/mol, activation volume from text */
p = 1.4; %/* grain size dependence */
mask = T<255;

A=3.9e-3;% /* MPa^-1.8 m^1.4 s^-1 */
n=1.8;
Q=49*1000.0; %/* J/mol */
egbs(mask) = A*(s(mask)/1e6).^n/d^p.*exp( (-Q + P*V)./(R*T(mask)) );
%     }else{
A = 3.0e26; %/* MPa^-1.8 m^1.4 S^-1 */
n = 1.8;
Q = 192*1000.0;
egbs(~mask) = A*(s(~mask)/1e6).^n/d^p.*exp( (-Q + P*V)./(R*T(~mask)) );

%   PetscScalar edisl;
%   {
%     PetscScalar A,n,Q;
if( T < 258 )
    A=4.0e5;% /* MPa^-4.0 s^-1 */
    n=4.0;
    Q=60*1000.0; %/* J/mol */
else
    A = 6.0e28; %/* MPa^-4.0 S^-1 */
    n = 4.0;
    Q = 180*1000.0; %/* Note that there is a mistake in G&K table 5. Should be 180 kJ/mol per text */
end
V = -13e-6; %/* m^3/mol, activation volume from text */
edisl = A*(s/1e6).^n.*exp( (-Q + P*V)./(R*T) );
%   }
%   PetscScalar ebs;
%   {
A = 5.5e7; %/* MPa^-2.4 s^-1 */
n = 2.4;
Q = 60.0*1000.0; %/* J/mol */
V = -13e-6; %/* m^3/mol, activation volume from text */
ebs = A*(s/1e6).^n.*exp( (-Q + P*V)./(R*T) );
%   }
etot = ediff + 1.0./(1.0./ebs + 1.0./egbs) + edisl;
%
etaeff = s./(2.0*etot);

mask = s < 1e0; %{/* if stress is really small (i.e. for first timestep), we cannot compute dislocation creep or gbs strain rates */
etaeff(mask) = 0.5./ediff0(mask);
