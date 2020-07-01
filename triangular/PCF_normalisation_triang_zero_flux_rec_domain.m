function [NORM]=PCF_normalisation_triang_zero_flux_rec_domain(Lx,Ly)
if mod(Lx,2)~=0
    error('Lx has to be even')
end
Rows=zeros(Lx*Ly,1);
Columns=Rows;
k=1;
for i=1:Lx
    for j=1:Ly
        if i<Lx
            Rows(k)=i+Lx*(j-1);
            Columns(k)=Rows(k)+1;
            k=k+1;
        end
        if j<Ly  %if not in last row
             if mod(j,2)~=0 %if on an odd row
                if mod(i,2)~=0 %all odd entries are connected
                    Rows(k)=i+Lx*(j-1);
                    Columns(k)=Rows(k)+Lx;
                    k=k+1;
                end
            else %if on an even row
                if mod(i,2)==0  %all even entries are connected
                    Rows(k)=i+Lx*(j-1);
                    Columns(k)=Rows(k)+Lx;
                    k=k+1;
                end
            end
        end
        
    end
end




A=sparse([Rows, Columns],[Columns, Rows], ones(1,2*size(Rows,1)),Lx*Ly,Lx*Ly);
%imagesc(A)
B=A;
spy(A)
A
storennz_A=zeros(1,min(Lx,Ly));
for m=1:min(Lx,Ly)
    storennz_A(m)=nnz(B);
    B=B*A;
end

NORM=zeros(1,min(Lx,Ly));
for m=3:min(Lx,Ly)
    NORM(m)=storennz_A(m)/2-storennz_A(m-2)/2;
end
NORM(1)=(nnz(A))/2;
NORM(2)=(storennz_A(2)-Lx*Ly)/2;


end