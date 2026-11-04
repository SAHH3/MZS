function C1 = C1(z, domain_var)

df1 = DF1(z, domain_var); df2 = DF2(z, domain_var); df3 = DF3(z, domain_var);

C1 = -df3 ./ (3 .* (df1.^3)) + 3 .* (df2.^2) ./ (4 .* (df1.^4));