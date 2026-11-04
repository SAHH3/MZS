
function[alpha] = alpha_heart(q, g, A)

alpha = sqrt(A/pi)*(1-g^2)*(1+q^2*g^2)*(1-q^4)                   ...
 /sqrt(1+q^4+2*g^6*q^8-g^4*q^2*(q^6+4*q^4+q^2-4)+g^2*(4*q^6-2*q^4-4*q^2));