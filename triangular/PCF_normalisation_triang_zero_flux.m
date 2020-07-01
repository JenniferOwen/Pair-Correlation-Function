function [NORM]=PCF_normalisation_triang_zero_flux(L)

Rows=1:6*L^2-1;
Columns=2:(6*L^2);
i=(2*L+1);
j=1;
while j<=L
    Rows(i)=0;
    Columns(i)=0;
    Rows(6*L^2-i)=0;
    Columns(6*L^2-i)=0;
    %Cumulative row position
    i=i+2*L+1+2*j;
    %Stores row of the hexagon
    j=j+1;
end %for

R=1:L;
R=1+(R-1).*(2*L+R-1);

k=1;
%For j=1 to middle of hexagonal
for j=1:L-1
    %From start of row to the end of the row
    for i=R(j):2:R(j)+2*L+(j-1)*2;
        Rows_down(k)=i;
        Columns_down(k)=i+2*L+1+(j-1)*2+1;
        k=k+1;
        Rows_down(k)=6*L^2-i+1;
        Columns_down(k)=(6*L^2+1-i)-2*L-1-(j-1)*2-1;
        k=k+1;
        
    end
end

for i=R(L):2:3*L^2;
    Rows_down(k)=i;
    Columns_down(k)=i+2*L+1+(L-1)*2;
    k=k+1;
end


[trash, trash, Rows]=find(Rows);
[trash, trash, Columns]=find(Columns);

[trash, trash, Rows_down]=find(Rows_down);
[trash, trash, Columns_down]=find(Columns_down);


A=sparse([Rows, Rows_down],[Columns, Columns_down], ones(1,size(Rows,2)+size(Rows_down,2)),6*L^2,6*L^2);
A=A+A';
%imagesc(A)
B=A;

storennz_A=zeros(1,L);
for m=1:L
    storennz_A(m)=nnz(B);
    B=B*A;   
end

NORM=zeros(1,L);
for m=3:L
        NORM(m)=storennz_A(m)/2-storennz_A(m-2)/2;
end
NORM(1)=(nnz(A))/2;
NORM(2)=(storennz_A(2)-6*L^2)/2;


end