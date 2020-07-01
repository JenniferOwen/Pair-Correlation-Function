comp=[];
%m=;
Lx=15;
Lz=16;
Ly=17
for m=1:7
    NORM=PCF_normalisation_3D_sparse_zero_flux_taxcab(Lx,Ly,Lz);
    %comp(m,:)=[Ly NORM(m) Lx*Ly*Lz*(2*m^2+1)-1/3*(2*m^3+m)*(Lx*Ly+Ly*Lz+Lz*Lx)+1/6*(m-1)*m^2*(1+m)*(Lx+Ly+Lz)-1/30*m*(m+1)*(m-1)*(6+5*m+m^2)+1/6*m*(m-1)*(m+1)*(m+2)];
    comp(m,:)=[Ly NORM(m) Lx*Ly*Lz*(2*m^2+1)-1/3*(2*m^3+m)*(Lx*Ly+Ly*Lz+Lz*Lx)+m^2/6*(m^2-1)*(Lx+Ly+Lz)-m^5/30+m^3/6-2*m/15];
    
end
[comp comp(:,3)-comp(:,2)]
