function [ xyvor ] = set_ic_MZS_chain_punc(N, eval, eshell, XY, eps, factor, varargin)

% Function to generate random i.c.s for main_doubly
% satisfying eval - eshell/2 < E < eval + eshell/2
% using a Markov chain method
% ensures min dist > factor * radius = factor * sqrt(8 * eps)

flag = 2;
if nargin == 7
    flag = varargin{1};
end

kappas=ones(1,N); kappas(N/2+1:N)=-1; NI=6;
parvec=[NI,kappas]';                      %parameters for RK4
% generate sets of vortices with (Px,Py) in range (epsilon) of 
% lattice points

nrec=100000; nvor=N*nrec;
xyvor=[];
xyposs=2*pi*rand(2*nvor,1);

NGO=1;
while (size(xyvor,1)==0)
    NGO
    i1=randperm(2*nvor); 
    xyposs=xyposs(i1);
    xys=xyposs(1:2*N);
    hvec0 = energy_doubly_MZS(0, xys, XY, eps, parvec);

    for ii=2*N+1:2:nrec                                     
        xys1=xys;
        jj=ceil(N*rand()); % index of vortex to replace
        xys1(jj)=xyposs(ii); xys1(jj+N)=xyposs(ii+1);
        hvec = energy_doubly_MZS(0, xys1, XY, eps,  parvec);

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