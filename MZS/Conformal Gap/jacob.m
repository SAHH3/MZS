function jacb = jacob(Z, domain_var)
%JACOB  Jacobian factor |F'(Z)| of the 'gap' conformal map z = F(Z) = 2iZ / (Z^2 - 1).
%
%   jacb = jacob(Z, domain_var)
%
%   Lengths in the physical plane are |F'(Z)| times lengths in the image
%   plane. Used in energy_bounded_MZS.m for the term (kappa^2 / 4 pi) log|F'(Z)|.
%
%   Z           image-plane position(s), any array size (evaluated elementwise)
%   domain_var  passed through to DF1.m (not used by the gap map)

jacb = abs(DF1(Z, domain_var));