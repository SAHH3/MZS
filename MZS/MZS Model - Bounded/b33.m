function b33 = b33(z, zp, domain_var)
%B33  Derivative of BB2(Z, Z') with respect to conj(Z) (BB2 depends on the
%     field point only through conj(Z)).
%
%   b33 = b33(z, zp, domain_var)
%
%     b33 = d/dconj(Z) BB2(Z, Z')
%         = ( 2 conj(Z)/(1 - conj(Z) Z') - F''(Z')/F'(Z') ) / ( (1 - conj(Z) Z') F'(Z') )^2
%
%   The image counterpart of b11.m. b33(Z, Z) is used for the self term.
%
%   z           image-plane field point(s) Z
%   zp          image-plane source point(s) Z'
%   domain_var  conformal map parameters, passed to DF1/DF2

zbar = conj(z);
df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

b33 = (2.*zbar ./ (1 - zbar.*zp) - df2./df1) ./ (((1 - zbar.*zp) .* df1).^2);