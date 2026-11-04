function [bigz]=xyvors(nvor,domain_var)
%XYVORS  nvor random positions, uniformly distributed over the physical
%        heart domain, returned in image-plane (unit-disc) coordinates.
%
%   bigz = xyvors(nvor, domain_var)
%
%   Draws 4*nvor points uniformly in a box slightly larger than the heart,
%   maps them back to the unit disc with inverse_map.m (which drops the
%   points outside the heart), and keeps nvor of them at random. Uniform in
%   the physical plane means NOT uniform in the image plane.
%   Errors if fewer than nvor of the 4*nvor points land inside the heart.
%
%   nvor        number of positions wanted
%   domain_var  [q, g], the heart shape parameters
%   bigz        1 x nvor complex image-plane positions
%   Needs boundary.m, so works with the Conformal Heart folder only.

     % bounding box of the heart, enlarged by 10%
     [x_bound, y_bound] = boundary(200,domain_var);
     ymax=1.1*max(abs(y_bound));
     xmax=1.1*max(abs(x_bound));

     x=-xmax+2*xmax*rand(1,fix(4*nvor));
     y=-ymax+2*ymax*rand(1,fix(4*nvor));
     z=x+1i*y;

     bigz=inverse_map(z,domain_var);     % only points inside the heart are returned
     pl1=randperm(length(bigz));
     bigz=bigz(pl1(1:nvor));

end