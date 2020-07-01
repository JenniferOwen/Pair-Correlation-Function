function [NORM]=PCF_normalisation_3D_sparse_zero_flux_taxcab(Lx,Ly,Lz)


%USES ADJACENCY MATRIX in sparse format
%Much faster than without sparse format

%This function calculates the normalisation part of the PCF with the
%taxi-cab metric using adjacency matrices.

%L is a rectangular lattice

%%ADJACENCY MATRIX BUILD
%Generates the adjacency matrix for lattice of size Lx by Ly in sparse format

%Include all pairs connected on the right


% right_pairs_r=zeros(Lz*Ly*(Lx-1)+Lz*Lx*(Ly-1)+Lx*Ly*(Lz-1),1);
% right_pairs_c=zeros(Lz*Ly*(Lx-1)+Lz*Lx*(Ly-1)+Lx*Ly*(Lz-1),1);
k=1;
right_pairs_periodic_r=zeros(Lz*Ly+Lx*Ly+Lx*Lz,1);
right_pairs_periodic_c=zeros(Lz*Ly+Lx*Ly+Lx*Lz,1);

%x-y plane
for h=0:Lz-1
    for i=0:Ly-1
        for j=1:Lx-1 %Every element of each row except for the last
            right_pairs_r(k)=i*Lx+j+h*Lx*Ly;
            right_pairs_c(k)=i*Lx+j+h*Lx*Ly+1;
            k=k+1;
        end
    end
end


%y-z plane
for h=1:Lz-1
    for i=1:Ly*Lx
        right_pairs_r(k)=i+h*Ly*Lx;
        right_pairs_c(k)=i+(h-1)*Ly*Lx;
        k=k+1;
    end
end

%x-z plane
for h=0:Lz-1
    for i=1:Lx
        for j=0:Ly-2 %Every element of each row except for the last
            right_pairs_r(k)=i+Lx*j+Lx*Ly*h;
            right_pairs_c(k)=i+Lx*j+Lx*Ly*h+Lx;
            k=k+1;
        end
    end
end
% 
% %Periodic x-y
% for i=1:Lx*Ly
%             right_pairs_r(k)=i;
%             right_pairs_c(k)=i+Lx*Ly*(Lz-1);
%             k=k+1;
% end
% 
% %Periodic x-z
% for i=1:Lx
%     for h=0:Lz-1
%             right_pairs_r(k)=i+Lx*Ly*h;
%             right_pairs_c(k)=i+Lx*Ly*h+Lx*(Ly-1);
%             k=k+1;
%     end
% end
% 
% % %Periodic y-z
%     for j=0:Ly-1
%         for h=0:Lz-1
%                 right_pairs_r(k)=1+Lx*j+Lx*Ly*h;
%                 right_pairs_c(k)=1+Lx*j+Lx*Ly*h+Lx-1;
%                 k=k+1;
%                 
%         end
%     end

%
%
% %x-z plane
% for h=0:Ly-1
%     for i=0:Lz-1
%         for j=1:Lx-1 %Every element of each row except for the last
%             right_pairs_r(k)=i*Lx+j+h*Lx*Lz;
%             right_pairs_c(k)=i*Lx+j+h*Lx*Lz+1;
%             k=k+1;
%         end
%     end
% end
%
%
% %y-z plane
% for h=0:Lx-1
%     for i=0:Ly-1
%         for j=1:Lz-1 %Every element of each row except for the last
%             right_pairs_r(k)=i*Lz+j+h*Lz*Lx;
%             right_pairs_c(k)=i*Lz+j+h*Lz*Lx+1;
%             k=k+1;
%         end
%     end
% end


%
% Periodic connections
% for i=1:Ly
%     right_pairs_periodic_r(i)=i*Lx;
%     right_pairs_periodic_c(i)=(i-1)*Lx+1;
% end
%
% Include all pairs connected below
% down_pairs_r=zeros((Ly-1)*Lx,1);
% down_pairs_c=zeros((Ly-1)*Lx,1);
%
% down_pairs_periodic_r=zeros(Lx,1);
% down_pairs_periodic_c=zeros(Lx,1);
%
%
% for i=1:Lx*(Ly-1)  %Up until last row
%     down_pairs_r(i)=i;
%     down_pairs_c(i)=i+Lx;
% end
%
% Periodic connections
% for i=1:Lx
%     down_pairs_periodic_r(i)=(Ly-1)*Lx+i;
%     down_pairs_periodic_c(i)=i;
% end
%
% Include all diagonal pairs right
% diag_right_pairs_r=zeros((Ly-1)*(Lx-1)-1,1);
% diag_right_pairs_c=zeros((Ly-1)*(Lx-1)-1,1);
%
% k=1;
% for i=0:Ly-2
%     for j=1:Lx-1
%         diag_right_pairs_r(k)=i*Lx+j;
%         diag_right_pairs_c(k)=i*Lx+j+Lx+1;
%         k=k+1;
%     end
% end
% diag_right_pairs_r(k)
%
% diag_right_pairs_periodic_r=zeros(Ly,1);
% diag_right_pairs_periodic_c=zeros(Ly,1);
%
% Periodic connections
% for i=1:Ly-1
%     diag_right_pairs_periodic_r(i)=Lx*i;
%     diag_right_pairs_periodic_c(i)=Lx*i+1;
% end
%
% for i=1:Lx-1
%     diag_right_pairs_periodic_r(Ly+i)=Lx*(Ly-1)+i;
%     diag_right_pairs_periodic_c(Ly+i)=i+1;
% end
% diag_right_pairs_periodic_r(Ly)=Lx*Ly;
% diag_right_pairs_periodic_c(Ly)=1;
%
% Include all diagonal pairs left
% diag_left_pairs_r=zeros((Ly-1)*(Lx-1)-1,1);
% diag_left_pairs_c=zeros((Ly-1)*(Lx-1)-1,1);
%
% k=1;
% for i=0:Ly-2
%     for j=2:Lx
%         diag_left_pairs_r(k)=i*Lx+j;
%         diag_left_pairs_c(k)=i*Lx+j+Lx-1;
%         k=k+1;
%     end
% end
%
%
% diag_left_pairs_periodic_r=zeros(Ly+Lx-1,1);
% diag_left_pairs_periodic_c=zeros(Ly+Lx-1,1);
%
% Periodic connections
% for i=0:Ly-2
%     diag_left_pairs_periodic_r(i+1)=Lx*i+1;
%     diag_left_pairs_periodic_c(i+1)=Lx*(i+2);
% end
%
% for i=1:Lx-1
%     diag_left_pairs_periodic_r(Ly+i)=Lx*(Ly-1)+1+i;
%     diag_left_pairs_periodic_c(Ly+i)=i;
% end
% diag_left_pairs_periodic_r(Ly)=Lx*(Ly-1)+1;
% diag_left_pairs_periodic_c(Ly)=Lx;
%
% Include all pairs and incorporate symmetry.
Rows=[right_pairs_r];
Columns=[right_pairs_c];
%
% Periodic conditions
% Rows=[Rows; diag_right_pairs_periodic_r;diag_left_pairs_periodic_r;right_pairs_periodic_r;down_pairs_periodic_r];
% Columns=[Columns; diag_right_pairs_periodic_c;diag_left_pairs_periodic_c;right_pairs_periodic_c;down_pairs_periodic_c];
%
%
R=[Rows;Columns];
C=[Columns; Rows];
%
A=sparse(R,C,ones(2*size(C,2),1),Lx*Ly*Lz,Lx*Ly*Lz);
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
        NORM(m)=(storeNNZ(m)-storeNNZ(m-2))/2;
end

if floor((min(Lx,Ly)-1)/2)>=2
NORM(2)=(storeNNZ(2)-Lx*Ly*Lz)/2;
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