function [NORM]=PCF_normalisation_rectangular_sparse_periodic(Lx,Ly)
%USES ADJACENCY MATRIX in sparse format

%Much faster than without sparse format

%This function calculates the normalisation part of the PCF with the
%taxi-cab metric using adjacency matrices.

%L is a rectangular lattice

%%ADJACENCY MATRIX BUILD
%Generates the adjacency matrix for lattice of size Lx by Ly in sparse format

%Include all pairs connected on the right 
right_pairs_r=zeros(Lx*Ly-Lx,1);
right_pairs_c=zeros(Lx*Ly-Lx,1);
k=1;
right_pairs_periodic_r=zeros(Ly,1);
right_pairs_periodic_c=zeros(Ly,1);

for i=0:Ly-1
    for j=1:Lx-1 %Every element of each row except for the last
    right_pairs_r(k)=i*Lx+j;
    right_pairs_c(k)=i*Lx+j+1;
    k=k+1;
    end
end

%Periodic connections
for i=1:Ly
    right_pairs_periodic_r(i)=i*Lx;
    right_pairs_periodic_c(i)=(i-1)*Lx+1;
end

%Include all pairs connected below
down_pairs_r=zeros((Lx-1)*Ly,1);
down_pairs_c=zeros((Lx-1)*Ly,1);

down_pairs_periodic_r=zeros(Lx,1);
down_pairs_periodic_c=zeros(Lx,1);


for i=1:(Lx-1)*Ly  %Up until last row
    down_pairs_r(i)=i;
    down_pairs_c(i)=i+Lx;
end

%Periodicconnections
for i=1:Lx
    down_pairs_periodic_r(i)=Ly*(Lx-1)+i;
    down_pairs_periodic_c(i)=i;
end




%Include all pairs and incorporate symmetry.
R=[right_pairs_r; right_pairs_periodic_r; down_pairs_r; down_pairs_periodic_r; right_pairs_c; right_pairs_periodic_c; down_pairs_c; down_pairs_periodic_c];
C=[right_pairs_c; right_pairs_periodic_c; down_pairs_c; down_pairs_periodic_c; right_pairs_r; right_pairs_periodic_r; down_pairs_r; down_pairs_periodic_r];

A=sparse(R,C,ones(size(C,1),1),Lx*Ly,Lx*Ly);
spy(A)
%%Fill norm(m)=A^m-A^(m-2)
%We only care about if there is a single route, and not how many there are
%in total
NORM=zeros(1,min(Lx,Ly));
storeNNZ=zeros(1,min(Lx,Ly));

%%Save the number of non zero elements of A^r for all r=1:min(Lx,Ly)
%storeNNZ(r)=nnz(A^r)
B=A;
for r=1:min(Lx,Ly)
    storeNNZ(r)=nnz(B);
    B=A*B;
end

for m=3:min(Lx,Ly)
        NORM(m)=(storeNNZ(m)-storeNNZ(m-2))/2;
end

if min(Lx,Ly)>=2
NORM(2)=(storeNNZ(2)-Lx*Ly)/2;
end
NORM(1)=storeNNZ(1)/2;
end