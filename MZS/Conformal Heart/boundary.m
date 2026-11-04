function[xs, ys] = boundary(NP,domain_var)
%BOUNDARY  Points on the boundary of the heart-shaped domain.
%
%   [xs, ys] = boundary(NP, domain_var)
%
%   Returns xs + i*ys = F(exp(i*t)) for t = 0 : 2*pi/NP : 2*pi, i.e. the image
%   of the boundary |Z| = 1 under the heart map F (see heart_map.m). That is NP+1
%   points, with the first equal to the last so the curve closes when plotted.
%   Used for plotting and by xyvors.m to find a box enclosing the domain.
%
%   NP          number of boundary intervals
%   domain_var  [q, g], the heart shape parameters

q=domain_var(1); g=domain_var(2); a=alpha_heart(q,g,pi);
t = (0: (2*pi)/NP: 2*pi);                    % angle around the boundary |Z| = 1
% real and imaginary parts of F(exp(i*t)), written out explicitly
xs = -cos(t).*a.*(2*sin(t)*g*q^2-1+q^2)./ ...
	 ((-2*q*cos(t)+1+q^2).*(2*q*cos(t)+1+q^2).*(g^2+2*sin(t)*g+1));
ys = -a*(-sin(t)-sin(t)*q^2-g+2*g*q^2*cos(t).^2-g*q^2)./ ...
	 ((-2*q*cos(t)+1+q^2).*(2*q*cos(t)+1+q^2).*(g^2+2*sin(t)*g+1));

end