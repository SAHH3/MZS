function jacb = jacob(Z, domain_var)

% JACOB: 
% returns the value of the jacobian |F'(Z)| 
% of the 'gap' mapping z = F(Z) = 2iZ / (Z^2 - 1)
% note that Z will typically be vector valued

jacb = abs(DF1(Z, domain_var));