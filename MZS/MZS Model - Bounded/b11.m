function b11 = b11(z, zp, domain_var)

df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

b11 = (2./(z - zp) - df2./df1) ./ (((z - zp) .* df1).^2);