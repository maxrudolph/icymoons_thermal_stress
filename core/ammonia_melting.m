function [tm,dtdx] = ammonia_melting(X)
% melting point of ammonian-water mixtures
% X is the mass fraction of ammonia
% equation (3a) from Croft et al. (1988)
% valid for X<0.329
%
assert(X<0.329);
tm=273.15-53.07*X-1651.4*X^2+11842.0*X^3-46269.0*X^4+56695.0*X^5;
dtdx=0.0 - 53.07 - 2*1651.4*X + 3*11842.0*X^2 - 4*46269.0*X^3 + 5*56695.0*X^4;
end