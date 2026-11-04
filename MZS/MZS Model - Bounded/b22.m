function b22 = b22(z, zp, domain_var)

df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var);

b22 = (2./(z - zp).^2 + 3.*df2 ./ (df1 .* (z - zp)) + 3.*(df2./df1).^2 - df3./df1) ./ ((z - zp) .* (df1.^2));