function [ jacb ] = jacob( bigz, domain_var )
%JACOB  Jacobian factor |F'(Z)| of the 'heart' conformal map.
%
%   jacb = jacob(bigz, domain_var)
%
%     |F'(Z)| = | a (1 + q^2 Z^2 - 2 i g q^2 Z^3) / ((1 - q^2 Z^2)(1 - i g Z))^2 |
%
%   (the same expression as DF1.m). Lengths in the physical plane are |F'(Z)|
%   times lengths in the image plane. Used in energy_bounded_MZS.m for the
%   term (kappa^2 / 4 pi) log|F'(Z)|.
%
%   bigz        image-plane position(s), any array size (evaluated elementwise)
%   domain_var  [q, g], the heart shape parameters

q=domain_var(1); g=domain_var(2); a=alpha_heart(q,g,pi);
jacb=abs(a*(1+q^2*bigz.^2-2*1i*q^2*g*bigz.^3)   ...
     ./((1-q^2*bigz.^2).^2.*(1-1i*g*bigz).^2));

end