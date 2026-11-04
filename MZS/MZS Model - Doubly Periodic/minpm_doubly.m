

function [ minpm, i1, i2 ] = minpm_doubly( xyvor, parvec , flag )
%
% function to find minimum value of distances between opposite or same 
% signed vortices. i1 and i2 are indices of the vortices
%
N=size(xyvor,1)/2;                   %no vortices  
kappas=parvec(2:N+1);                %extract params.
mat1=repmat(1:N,N,1); mat2=mat1';
if flag==0 % opposite signed
    utripm=find((mat1>mat2) & (kappas*kappas'==-1)); %p/m distances
elseif flag==1
    utripm=find((mat1>mat2) & (kappas*kappas'==1)); %p/p,m/m distances    
else
    utripm=find((mat1>mat2));   % flag>1 all distances (for gp)
end

xdists=repmat(xyvor(1:N),1,N)-repmat(xyvor(1:N),1,N)'; 
ydists=repmat(xyvor(N+1:2*N),1,N)-repmat(xyvor(N+1:2*N),1,N)'; 
xdists=mod(xdists(utripm)+pi ,2*pi)-pi; %map to [-pi,pi)
ydists=mod(ydists(utripm)+pi ,2*pi)-pi;

[minpm, ind]=min(sqrt(xdists.^2+ydists.^2));
i1=mat1(utripm(ind)); i2=mat2(utripm(ind));

end