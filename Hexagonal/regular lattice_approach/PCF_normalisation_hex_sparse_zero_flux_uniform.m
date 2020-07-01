function [NORM]=PCF_normalisation_hex_sparse_zero_flux_uniform(Lx,Ly)
%USES ADJACENCY MATRIX in sparse format

%Much faster than without sparse format

%This function calculates the normalisation part of the PCF with the
%taxi-cab metric using adjacency matrices.

%L is a rectangular lattice

%%ADJACENCY MATRIX BUILD
%Generates the adjacency matrix for lattice of size Lx by Ly in sparse format

%Include all pairs connected on the right 
right_pairs_r=zeros(Ly*(Lx-1),1);
right_pairs_c=zeros(Ly*(Lx-1),1);
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
down_pairs_r=zeros((Ly-1)*Lx,1);
down_pairs_c=zeros((Ly-1)*Lx,1);

down_pairs_periodic_r=zeros(Lx,1);
down_pairs_periodic_c=zeros(Lx,1);


for i=1:Lx*(Ly-1)  %Up until last row
    down_pairs_r(i)=i;
    down_pairs_c(i)=i+Lx;
end

%Periodicconnections
for i=1:Lx
    down_pairs_periodic_r(i)=Ly*(Lx-1)+i;
    down_pairs_periodic_c(i)=i;
end

%Include all diagonal pairs right
 diag_right_pairs_r=zeros((Ly-1)*(Lx-1)-1-Lx/2*Ly,1);
 diag_right_pairs_c=zeros((Ly-1)*(Lx-1)-1-Lx/2*Ly,1);

k=1;
for i=0:Ly-2
    for j=1:Lx-1
       if mod(j,2)==1 %only even numbers
        diag_right_pairs_r(k,1)=i*Lx+j;
        diag_right_pairs_c(k,1)=i*Lx+j+Lx+1;
        k=k+1;
       end
    end
end

%Include all diagonal pairs left
 diag_left_pairs_r=zeros((Ly-1)*(Lx-1)-1-Lx/2*Ly,1);
 diag_left_pairs_c=zeros((Ly-1)*(Lx-1)-1-Lx/2*Ly,1);

k=1;
for i=0:Ly-2
    for j=2:Lx
        if mod(j,2)==1
        diag_left_pairs_r(k,1)=i*Lx+j;
        diag_left_pairs_c(k,1)=i*Lx+j+Lx-1;
        k=k+1;
        end
    end
end


%Include all pairs and incorporate symmetry.
Rows=[right_pairs_r; down_pairs_r; diag_left_pairs_r; diag_right_pairs_r];
Columns=[right_pairs_c; down_pairs_c; diag_left_pairs_c; diag_right_pairs_c];

R=[Rows;Columns];
C=[Columns; Rows];

A=sparse(R,C,ones(size(C,1),1),Lx*Ly,Lx*Ly);
%spy(A)
%%Fill norm(m)=A^m-A^(m-2)
%We only care about if there is a single route, and not how many there are
%in total
NORM=zeros(1,min(Lx,Ly)-1);
storeNNZ=zeros(1,min(Lx,Ly)-1);

%%Save the number of non zero elements of A^r for all r=1:min(Lx,Ly)
%storeNNZ(r)=nnz(A^r)
B=A;
for r=1:min(Lx,Ly)-1
    storeNNZ(r)=nnz(B);
    B=A*B;
end

for m=3:min(Lx,Ly)-1
        NORM(m)=(storeNNZ(m)-storeNNZ(m-1))/2;
end

if min(Lx,Ly)-1>=2
NORM(2)=(storeNNZ(2)-storeNNZ(1)-Lx*Ly)/2;
end
NORM(1)=storeNNZ(1)/2;
end
%
%Y_c=[];
%for i=1:Ly
%    Y_c=[Y_c;i*ones(Lx,1)];
%end
%X_c=repmat(1:Lx,1,Ly)';
%gplot(A,[X_c,Y_c])
%axis([0 Lx+1 0 Ly+1])
