function BB1 = BB1(z, zp, domain_var)

df1 = DF1(zp, domain_var); df2 = DF2(zp, domain_var);

BB1 = -(1./(z - zp) - df2./df1) ./ ((z - zp) .* (df1.^2));