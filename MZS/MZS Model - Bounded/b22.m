function b22 = b22(z, zp, domain_var)
%B22  Derivative of BB1 with respect to its SOURCE point, with the source at z.
%
%   b22 = b22(z, zp, domain_var)
%
%     b22(Z, Z') = d/dZ [ BB1(Z', Z) ]
%       = ( 2/(Z - Z')^2 + 3 F''/(F' (Z - Z')) + 3 (F''/F')^2 - F'''/F' ) / ( (Z - Z') F'^2 )
%
%   with F', F'', F''' evaluated at Z. Used in uvUV_bounded_MZS.m for the
%   velocity of vortex i (at Z) due to its own ellipticity interacting with
%   point vortex j (at Z').
%
%   z           image-plane position(s) Z (the source of BB1, being moved)
%   zp          image-plane position(s) Z'
%   domain_var  conformal map parameters, passed to DF1/DF2/DF3

df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var);

b22 = (2./(z - zp).^2 + 3.*df2 ./ (df1 .* (z - zp)) + 3.*(df2./df1).^2 - df3./df1) ./ ((z - zp) .* (df1.^2));