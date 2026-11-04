function BB1 = BB1(z, zp, domain_var)
%BB1  Second derivative of the free-space part of the Green's function with
%     respect to the PHYSICAL source position w' = F(Z').
%
%   BB1 = BB1(z, zp, domain_var)
%
%     BB1 = d^2/dw'^2 log(Z - Z')
%         = -( 1/(Z - Z') - F''(Z')/F'(Z') ) / ( (Z - Z') F'(Z')^2 )
%
%   using d/dw' = (1/F'(Z')) d/dZ'. Together with BB2.m this gives the strain
%   that the vortex at Z' feels from a point vortex at Z, which couples to the
%   ellipticity of the vortex at Z' (see GQX.m and uvUV_bounded_MZS.m).
%
%   z           image-plane field point(s) Z
%   zp          image-plane source point(s) Z'
%   domain_var  conformal map parameters, passed to DF1/DF2

df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

BB1 = -(1./(z - zp) - df2./df1) ./ ((z - zp) .* (df1.^2));