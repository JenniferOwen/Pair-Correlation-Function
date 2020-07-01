comp=[];
%m=
Lx=20;
Lz=20;
Ly=20;
for m=1:8
   % Lx=Ly;
    NORM=PCF_normalisation_3D_sparse_zero_flux_uniform(Lx,Ly,Lz);
    %comp(Ly,:)=[Ly NORM(m) Lx*Ly*Lz*(2*m^2+1)-(m+2/3*m*(m-1)*(m+1))*(Lx*Ly+Ly*Lz+Lz*Lx)+1/6*(m-1)*m^2*(1+m)*(Lx+Ly+Lz)-1/30*m*(m+1)*(m-1)*(6+5*m+m^2)+1/6*m*(m-1)*(m+1)*(m+2)];
    comp(m,:)=[m NORM(m) Lx*Ly*Lz*(12*m^2+1)-m*(8*m^2+1)*(Lx*Ly+Ly*Lz+Lz*Lx)+(m^2*(5*m^2+1))*(Lx+Ly+Lz)-3*m^5-m^3];
    %-6*(3/4*m^2*(m+1)^2-3/2*m*(m-1))
end
G=[comp (comp(:,3)-comp(:,2))]
%plot(comp(:,3)-comp(:,2),'*')
% hold on
% plot(comp(:,3),'o')
% hold off


% A=[];
% dim=8;
% for i=1:dim
%     for j=1:dim
%         A(i,j)=i^(dim-j);
%     end
% end
% 
% A\(-G(1:dim,end))