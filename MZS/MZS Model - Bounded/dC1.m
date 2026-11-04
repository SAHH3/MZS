function dC1 = dC1(z, domain_var)
%DC1  Derivative dC1/dZ of C1.m.
%
%   dC1 = dC1(z, domain_var)
%
%     dC1/dZ = -F''''/(3 F'^3) + 5 F''' F''/(2 F'^4) - 3 F''^3/F'^5
%
%   Used in uvUV_bounded_MZS.m for the velocity a vortex gets from the h7 term.
%
%   z           image-plane position(s) Z
%   domain_var  conformal map parameters, passed to DF1/DF2/DF3/DF4

df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var); df4 = DF4(z, domain_var);

dC1 = -df4 ./ (3 .* (df1.^3)) + 5.*df3.*df2 ./ (2 .* (df1.^4)) - 3 .* (df2.^3) ./ (df1.^5);