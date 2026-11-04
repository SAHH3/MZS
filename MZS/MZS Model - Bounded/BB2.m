function BB2 = BB2(z, zp, domain_var)
%BB2  Second derivative of the image part of the Green's function with
%     respect to the PHYSICAL source position w' = F(Z').
%
%   BB2 = BB2(z, zp, domain_var)
%
%     BB2 = d^2/dw'^2 [ -log(1 - conj(Z) Z') ]
%         = ( conj(Z)/(1 - conj(Z) Z') - F''(Z')/F'(Z') ) conj(Z)
%           / ( (1 - conj(Z) Z') F'(Z')^2 )
%
%   The image counterpart of BB1.m (Re[-log(1 - conj(Z) Z')] is the image
%   term of the disc Green's function, written so it is analytic in Z').
%   BB2(Z, Z) is the strain a vortex feels from its own image (gQX0.m).
%
%   z           image-plane field point(s) Z
%   zp          image-plane source point(s) Z'
%   domain_var  conformal map parameters, passed to DF1/DF2

zbar = conj(z);
df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

BB2 = (zbar./(1 - zbar.*zp) - df2./df1) .* zbar ./ ((1 - zbar.*zp) .* (df1.^2));