function [NORM]=PCF_normalisation_sparse(Lx,Ly)
%USES ADJACENCY MATRIX

%This function calculates the normalisation part of the PCF with the
%taxi-cab metric using adjacency matrices in sparse format.

%Require L to be a square matrix

c=[];
r=[];
for i=1:Lx*Ly
    if mod(i,Lx)~=0
        r=[r; i];
        c=[c; i+1];
    end
    if i<(Lx-1)*Ly+1
        r=[r; i];
        c=[c; i+Lx];
    end
end


R=[r; c];
C=[c; r];


A=sparse(R,C,ones(size(C)));
id=sparse(1:Lx*Ly,1:Lx*Ly,ones(Lx*Ly,1));

NORM=zeros(1,max(Lx,Ly));
NORM(1)=nnz(A)/2;
NORM(2)=(nnz(A^2)-nnz(id))/2;
for m=3:max(Lx,Ly)
    
    NORM(m)=(nnz(A^m)-nnz(A^(m-2)))/2;
    
end

end