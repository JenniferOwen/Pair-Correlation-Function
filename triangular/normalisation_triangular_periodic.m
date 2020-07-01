function [NORM]=normalisation_triangular_periodic(L);


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

%Periodic boundary conditions from top to bottom of hexagon
rows_periodic_down=2:2:2*L+1;
columns_periodic_down=6*L^2-(2*L+1)+2:2:6*L^2;
k=1;
%Periodic boundary conditions from top left to bottom right of hexagon
for j=1:L
    rows_periodic_left_down(k)=R(j);
    columns_periodic_left_down(k)=6*L^2+1-R(j);
    k=k+1;
end
columns_periodic_left_down=flip(columns_periodic_left_down);
k=1;
%Periodic boundary conditions from top right to bottom left of hexagon
for j=1:L
    rows_periodic_right_down(k)=R(j)+2*L+(j-1)*2;
    columns_periodic_right_down(k)=6*L^2-2*L-((j-1)*2)-R(j)+1;
    k=k+1;
end
 columns_periodic_right_down=flip(columns_periodic_right_down);

periodic_rows=[rows_periodic_right_down,rows_periodic_left_down, rows_periodic_down];
periodic_columns=[columns_periodic_right_down,columns_periodic_left_down, columns_periodic_down];

A=sparse([Rows, Rows_down,periodic_rows],[Columns, Columns_down,periodic_columns], ones(1,size(Rows,2)+size(Rows_down,2)+size(periodic_rows,2)),6*L^2,6*L^2);
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