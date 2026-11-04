function b44 = b44(z, zp, domain_var)
%B44  Derivative of BB2 with respect to its SOURCE point, with the source at z.
%
%   b44 = b44(z, zp, domain_var)
%
%     b44(Z, Z') = d/dZ [ BB2(Z', Z) ],  with w = conj(Z'), D = 1 - Z w:
%       = ( 2 (w/D)^2 - 3 w F''/(F' D) + 3 (F''/F')^2 - F'''/F' ) w / ( D F'^2 )
%
%   with F', F'', F''' evaluated at Z. The image counterpart of b22.m.
%   b44(Z, Z) is used for the self term.
%
%   z           image-plane position(s) Z (the source of BB2, being moved)
%   zp          image-plane position(s) Z'
%   domain_var  conformal map parameters, passed to DF1/DF2/DF3

zpbar = conj(zp);
df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var);

b44 = (2 .* (zpbar ./ (1 - z.*zpbar)).^2 - 3.*zpbar.*df2 ./ (df1 .* (1 - z.*zpbar)) + 3.*(df2./df1).^2 - df3./df1) .* zpbar ./ ((1 - z.*zpbar) .* (df1.^2));