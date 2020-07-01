comp=[];
Ly=30;
Lx=25;
for m=1:min(Lx,Ly)-1
    k=mod(m,2);
    NORM=PCF_normalisation_hex_sparse_zero_flux_uniform(Lx,Ly);
    
    %comp(m,:)=[m NORM(m) 3*Lx*Ly*m-1/4*(7*m^2+m-2*k)*Lx-2*m^2*Ly+1/3*(m^3+3*m^2-m-4*k^3-3*k^2*m+6*k*m^2-6*k*m+k)];
        comp(m,:)=[m NORM(m) 3*Lx*Ly*m-1/4*(7*m^2+k)*Lx-2*m^2*Ly+11/12*m^3-(2-3*k)*m/12];
  

end
G=[comp comp(:,3)-comp(:,2)]

% A=[];
% dim=4;
% for i=1:dim
%     for j=1:dim
%         A(i,j)=i^(dim-j);
%     end
% end
% 
% A\(-G(1:2:2*dim,end))