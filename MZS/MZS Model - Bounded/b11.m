function b11 = b11(z, zp, domain_var)
%B11  Derivative of BB1(Z, Z') with respect to the field point Z.
%
%   b11 = b11(z, zp, domain_var)
%
%     b11 = d/dZ BB1(Z, Z') = ( 2/(Z - Z') - F''(Z')/F'(Z') ) / ( (Z - Z') F'(Z') )^2
%
%   Used in uvUV_bounded_MZS.m for the velocity of point vortex i (at Z)
%   due to the ellipticity of vortex j (at Z').
%
%   z           image-plane field point(s) Z
%   zp          image-plane source point(s) Z'
%   domain_var  conformal map parameters, passed to DF1/DF2

df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

b11 = (2./(z - zp) - df2./df1) ./ (((z - zp) .* df1).^2);