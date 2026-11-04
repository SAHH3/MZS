function C1 = C1(z, domain_var)
%C1  Conformal-map correction to the strain a vortex feels from itself.
%
%   C1 = C1(z, domain_var)
%
%     C1(Z) = lim_{Z' -> Z} d^2/dw'^2 [ log(Z - Z') - log(w - w') ]
%           = -F'''/(3 F'^3) + 3 F''^2/(4 F'^4),   w = F(Z), w' = F(Z')
%
%   The free-space part of the Green's function in the image plane differs
%   from the physical-plane log|w - w'|; this is the leftover second
%   derivative at the vortex itself. Zero for the unit disc. Gives the h7
%   term of energy_bounded_MZS.m.
%
%   z           image-plane position(s) Z
%   domain_var  conformal map parameters, passed to DF1/DF2/DF3

df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var);

C1 = -df3 ./ (3 .* (df1.^3)) + 3 .* (df2.^2) ./ (4 .* (df1.^4));