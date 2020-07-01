function [NORM]=PCF_normalisation_rectangular(Lx,Ly)
%USES ADJACENCY MATRIX

%This function calculates the normalisation part of the PCF with the
%taxi-cab metric using adjacency matrices.

%L is a rectangular lattice

%Create space to store A^m
storeA=NaN(Lx*Ly,Lx*Ly,2*max(Lx,Ly));

%A^0=Id
storeA(:,:,1)=eye(Lx*Ly);

%Generates the adjacency matrix
A=zeros(Lx*Ly);
for i=1:Lx*Ly
    if mod(i,Lx)~=0
        A(i,i+1)=1;
    end    
    if i<=(Ly-1)*Lx
        A(i,i+Lx)=1;
    end
end
size(A);
storeA(:,:,2)=A+A';

%Fills storeA with A^m
for m=3:min(Lx,Ly)+1
   storeA(:,:,m)=storeA(:,:,m-1)*storeA(:,:,2);
end

%%Fill norm(m)=A^m-A^(m-2)

%We only care about if there is a single route, and not how many there are
%in total
storeA=storeA>0;
NORM=zeros(1,min(Lx,Ly));
for m=2:min(Lx,Ly)
        NORM(m)=(nnz(storeA(:,:,m+1))-nnz(storeA(:,:,m-1)))/2;
end
NORM(1)=nnz(storeA(:,:,2))/2;
end