function [ xyvor ] = set_ic_MZS_chain_punc(N, eval, eshell, XY, eps, factor, varargin)
%SET_IC_MZS_CHAIN_PUNC  Random initial positions in the doubly periodic domain
%                       with a given energy and a minimum separation.
%
%   xyvor = set_ic_MZS_chain_punc(N, eval, eshell, XY, eps, factor)
%   xyvor = set_ic_MZS_chain_punc(N, eval, eshell, XY, eps, factor, flag)
%
%   Finds positions with eval - eshell/2 < E < eval + eshell/2 and with every
%   pair (selected by flag, see minpm_doubly.m) further apart than
%   factor * sqrt(8*eps(1)), i.e. factor times the radius of a circular
%   vortex of size eps(1).
%
%   Markov chain method: start from N random positions, then repeatedly move
%   one randomly chosen vortex to a new random position, accepting the move
%   only if it takes the energy towards eval. Stop as soon as both conditions
%   hold; if the supply of random positions runs out first, reshuffle and
%   start again (the printed NGO counts the restarts).
%
%   N             number of vortices (even); the first N/2 have k = +1, the rest -1
%   eval, eshell  centre and width of the energy shell
%   XY, eps       vortex shapes [X; Y] (2N x 1) and sizes (N x 1)
%   factor        minimum separation, in units of the vortex radius
%   flag          which pairs to check (default 2: all pairs)
%   xyvor         [x; y], 2N x 1
%   NB. NI = 6 and the circulations are set here, not taken from the main script.

flag = 2;
if nargin == 7
    flag = varargin{1};
end

kappas=ones(1,N); kappas(N/2+1:N)=-1; NI=6;
parvec=[NI,kappas]';                      %parameters for energy_doubly_MZS

nrec=100000; nvor=N*nrec;
xyvor=[];
xyposs=2*pi*rand(2*nvor,1);               % supply of random coordinates in [0, 2pi)

NGO=1;
while (size(xyvor,1)==0)
    NGO
    i1=randperm(2*nvor);
    xyposs=xyposs(i1);
    xys=xyposs(1:2*N);                    % starting positions
    hvec0 = energy_doubly_MZS(0, xys, XY, eps, parvec);

    for ii=2*N+1:2:nrec
        xys1=xys;
        jj=ceil(N*rand()); % index of vortex to replace
        xys1(jj)=xyposs(ii); xys1(jj+N)=xyposs(ii+1);
        hvec = energy_doubly_MZS(0, xys1, XY, eps,  parvec);

        % accept if the energy moved towards eval
        if (hvec-hvec0)*(eval-hvec0)>0
            xys=xys1; hvec0=hvec;
            [ii,hvec];
            [minpm,~,~] = minpm_doubly(xys, parvec, flag);

            if (hvec > eval - eshell/2 && hvec < eval + eshell/2 && minpm > factor * sqrt(8*eps(1)))
                xyvor=xys;
                return
            end
        end
    end
    NGO=NGO+1;
end