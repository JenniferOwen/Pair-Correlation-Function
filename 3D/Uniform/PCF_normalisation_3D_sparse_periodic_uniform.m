function [NORM]=PCF_normalisation_3D_sparse_zero_flux_uniform(Lx,Ly,Lz)

Columns=zeros(1,3*(Lx-1)*(Ly-1)*(Lz-1));
Rows=Columns;
h=1;
for i=1:Lx
    for j=1:Ly
        for k=1:Lz
            if i<Lx
                Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                Rows(h)=Columns(h)+1;
                h=h+1;
                if j<Ly
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)+Lx+1;
                    h=h+1;
                end
            else
               Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                Rows(h)=Columns(h)+1-Lx;
                h=h+1; 
                if j<Ly
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)+1;
                    h=h+1;
                else
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)+1-Lx*Ly;
                    h=h+1;
                end
            end
            if j<Ly
                Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                Rows(h)=Columns(h)+Lx;
                h=h+1;
                if k<Lz
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)+Lx*(1+Ly);
                    h=h+1;
                    if i<Lx
                        Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                        Rows(h)=Columns(h)+Lx*(1+Ly)+1;
                        h=h+1;
                    end
                end
            else
                Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                Rows(h)=Columns(h)-Lx*(Ly-1);
                h=h+1;
                if k<Lz
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)+Lx;
                    h=h+1;
                    if i<Lx
                        Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                        Rows(h)=Columns(h)+Lx+1;
                        h=h+1;
                    else
                        Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                        Rows(h)=Columns(h)+1;
                        h=h+1;
                    end
                else
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)-Lx*Ly*(Lz-1)-Lx*(Ly-1);
                    h=h+1;
                    if i<Lx
                        Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                        Rows(h)=Columns(h)+-Lx*Ly*(Lz-1)-Lx*(Ly-1)+1;
                        h=h+1;
                    else
                        Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                        Rows(h)=1;
                        h=h+1;
                    end
                end
            end
            if k<Lz
                Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                Rows(h)=Columns(h)+Lx*Ly;
                h=h+1;
                if i<Lx
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)+Lx*Ly+1;
                    h=h+1;
                else
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)+Lx*(Ly-1);
                    h=h+1;
                end
            else
                Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                Rows(h)=Columns(h)-Lx*Ly*(Lz-1);
                h=h+1;
                if i<Lx
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=Columns(h)-Lx*Ly*(Lz-1)+1;
                    h=h+1;
                else
                    Columns(h)=i+(j-1)*Lx+(k-1)*Lx*Ly;
                    Rows(h)=1+(j-1)*Lx;
                    h=h+1;
                end
            end
        end
    end
end
      


R=[Rows,Columns];
C=[Columns, Rows];
%

A=sparse(R,C,ones(size(C,1),1),Lx*Ly*Lz,Lx*Ly*Lz);
%%
% the formula for the periodic taxicam is 3*(2*m-1)*Lx*Ly*Lz
% spy(A)
% %Fill norm(m)=A^m-A^(m-2)
% We only care about if there is a single route, and not how many there are
% in total
NORM=zeros(1,floor((min(Lx,Ly)-1)/2));
storeNNZ=zeros(1,floor((min(Lx,Ly)-1)/2)-1);

%Save the number of non zero elements of A^r for all r=1:min(Lx,Ly)
%storeNNZ(r)=nnz(A^r)
B=A;
for r=1:floor((min(Lx,Ly)-1)/2)
    storeNNZ(r)=nnz(B);
    B=A*B;
end

for m=3:floor((min(Lx,Ly)-1)/2)
        NORM(m)=(storeNNZ(m)-storeNNZ(m-1))/2;
end

if floor((min(Lx,Ly)-1)/2)>=2
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