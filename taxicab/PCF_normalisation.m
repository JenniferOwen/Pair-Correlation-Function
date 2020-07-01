function [NORM]=PCF_normalisation(L)
%USES ADJACENCY MATRIX

%This function calculates the normalisation part of the PCF with the
%taxi-cab metric using adjacency matrices.

%Require L to be a square matrix

%Place to store A^m
%storeA(:,:,m)=A^(m+1)

storeA=NaN(L^2,L^2,2*L);

%Store the identity matrix
storeA(:,:,1)=eye(L^2);


%%Produce the adjacency matrix for a square lattice of size L by L
A=zeros(L^2);
for i=1:L^2
    if mod(i,L)~=0
        A(i,i+1,1)=1;
    end 
    if i<=L^2-L
        A(i,i+L,1)=1;
    end
end
%A is symmetric since if a is connected to b, b is connected to a.
storeA(:,:,2)=A+A';

%Fill storeA matrix with A^m
for m=3:2*L
    storeA(:,:,m)=storeA(:,:,m-1)*storeA(:,:,2);
end

storeA=storeA>0;
NORM=zeros(1,2*L);
%Calculate the norm to be norm(m)=A^m-A^(m-2)
for m=2:2*L-1
    NORM(m)=(nnz(storeA(:,:,m+1))-nnz(storeA(:,:,m-1)))/2;
end
NORM(1)=nnz(storeA(:,:,2))/2;
end