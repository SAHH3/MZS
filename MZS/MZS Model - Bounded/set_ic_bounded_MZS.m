function [ xyvor ] = set_ic_bounded_MZS(n, nvor, eval, eshell, XY, eps, parvec)
%SET_IC_BOUNDED_MZS  Random initial positions in the heart with a given energy.
%
%   xyvor = set_ic_bounded_MZS(n, nvor, eval, eshell, XY, eps, parvec)
%
%   Finds positions satisfying eval - eshell/2 < E < eval + eshell/2.
%   First makes nvor candidate positions, uniformly distributed over the
%   physical domain (xyvors.m). It then shuffles them, splits them into
%   floor(nvor/n) groups of n, and returns the first group whose energy is in
%   the shell, reshuffling as many times as needed. The printed NGO counts
%   the reshuffles.
%
%   n             number of vortices
%   nvor          number of candidate positions to generate
%   eval, eshell  centre and width of the energy shell
%   XY, eps       vortex shapes [X; Y] (2n x 1) and sizes (n x 1)
%   parvec        [kappas; domain_var]
%   xyvor         [x; y], 2n x 1, image-plane positions
%   Needs boundary.m, so works with the Conformal Heart folder only.

    xyvor=[]; domain_vec=parvec(n+1:length(parvec))';
    bigzs=xyvors(nvor,domain_vec);       % candidate image-plane positions
    nrec=floor(nvor/n); NGO=1;           % number of groups of n per shuffle

    while (size(xyvor,1)==0)
        NGO
        i1=randperm(nvor);

        for ii=1:nrec

            % positions of the ii-th group of n
            xys = [real(bigzs(i1((ii-1)*n+1:ii*n))), ...
                 imag(bigzs(i1((ii-1)*n+1:ii*n)))]';
            hvec = energy_bounded_MZS(0, xys, XY, eps, parvec);

            if (hvec>eval-eshell/2 && hvec<eval+eshell/2)
                xyvor=xys;
                return
            end
        end

	NGO=NGO+1;
    end
end