function [NORM]=PCF_normalisation_3D_sparse_zero_flux_uniform(Lx,Ly,Lz)

Lx=Lx+2;
Ly=Ly+2;
Lz=Lz+2;
Columns=zeros(1,26*(Lx-2)*(Ly-2)*(Lz-2));
Rows=Columns;
select=zeros(1,(Lx-2)*(Ly-2)*(Lz-2));
h=1;
hh=1;
for i=2:Lx-1
    for j=2:Ly-1
        for k=2:Lz-1
            select(hh)=(i+(j-1)*Lx+(k-1)*Lx*Ly);
            hh=hh+1;
            Columns(h:h+25)=(i+(j-1)*Lx+(k-1)*Lx*Ly)*ones(1,26);
            Rows(h:h+25)=Columns(h:h+25)+[1,-1,Lx,-Lx,Lx*Ly,-Lx*Ly,Lx+1,Lx-1,-Lx+1,-Lx-1,Lx*Ly+1,Lx*Ly-1,-Lx*Ly+1,-Lx*Ly-1,Lx*Ly+Lx,Lx*Ly-Lx,-Lx*Ly+Lx,-Lx*Ly-Lx,Lx*Ly+Lx+1,Lx*Ly+Lx-1,Lx*Ly-Lx+1,Lx*Ly-Lx-1,-Lx*Ly+Lx+1,-Lx*Ly+Lx-1,-Lx*Ly-Lx+1,-Lx*Ly-Lx-1];
            h=h+26;
        end
    end
end
      
R=[Rows,Columns];
C=[Columns,Rows];
A=sparse(R,C,ones(size(C,1),1),Lx*Ly*Lz,Lx*Ly*Lz);
Lx=Lx-2;
Ly=Ly-2;
Lz=Lz-2;
A=A(select,select);
A=double(A>0);
spy(A)
%%
% the formula for the periodic taxicam is 3*(2*m-1)*Lx*Ly*Lz
% spy(A)
% %Fill norm(m)=A^m-A^(m-2)
% We only care about if there is a single route, and not how many there are
% in total
NORM=zeros(1,floor((min(Lx,min(Ly,Lz))-1)/2));
storeNNZ=zeros(1,floor((min(Lx,min(Ly,Lz))-1)/2)-1);

%Save the number of non zero elements of A^r for all r=1:min(Lx,Ly)
%storeNNZ(r)=nnz(A^r)
B=A;
for r=1:floor((min(Lx,min(Ly,Lz))-1))
    storeNNZ(r)=nnz(B);
    B=A*B;
end

for m=3:floor((min(Lx,min(Ly,Lz))-1))
        NORM(m)=(storeNNZ(m)-storeNNZ(m-1))/2;
end

if floor((min(Lx,Ly)-1))>=2
NORM(2)=(storeNNZ(2)-storeNNZ(1)-Lx*Ly*Lz)/2;
end
NORM(1)=storeNNZ(1)/2;


%Y_c=[];
%for i=1:Ly
%    Y_c=[Y_c;i*ones(Lx,1)];
%end
%X_c=repmat(1:Lx,1,Ly)';
%gplot(A,[X_c,Y_c])
%axis([0 Lx+1 0 Ly+1])

end