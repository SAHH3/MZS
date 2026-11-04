
function [bigz]=xyvors(nvor,domain_var) 
       
%XYVORS:
%place nvor vortices in a bounded domain with params domain_var
%will fail sometimes if n is not reasonably large
%etavor returns complex coord in image domain (circle)

     [x_bound, y_bound] = boundary(200,domain_var);
     ymax=1.1*max(abs(y_bound));
     xmax=1.1*max(abs(x_bound));
      
     x=-xmax+2*xmax*rand(1,fix(4*nvor));
     y=-ymax+2*ymax*rand(1,fix(4*nvor));
     z=x+1i*y;
     
     bigz=inverse_map(z,domain_var);          
     pl1=randperm(length(bigz));  
     bigz=bigz(pl1(1:nvor));
    
end