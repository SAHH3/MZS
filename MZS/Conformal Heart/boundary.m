
function[xs, ys] = boundary(NP,domain_var)

% BOUNDARY:
% returns NP points on the boundary of the heart_shaped domain
% for use with plotting, set_ic_bounded etc. 

q=domain_var(1); g=domain_var(2); a=alpha_heart(q,g,pi);
t = (0: (2*pi)/NP: 2*pi);
xs = -cos(t).*a.*(2*sin(t)*g*q^2-1+q^2)./ ...
	 ((-2*q*cos(t)+1+q^2).*(2*q*cos(t)+1+q^2).*(g^2+2*sin(t)*g+1));
ys = -a*(-sin(t)-sin(t)*q^2-g+2*g*q^2*cos(t).^2-g*q^2)./ ...
	 ((-2*q*cos(t)+1+q^2).*(2*q*cos(t)+1+q^2).*(g^2+2*sin(t)*g+1));

end