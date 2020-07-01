comp=[];
m=6;
Lx=15;
Lz=15;
for Ly=15:20
    NORM=PCF_normalisation_3D_sparse_zero_periodic_taxcab(Lx,Ly,Lz);
    comp(Ly,:)=[Ly NORM(m) Lx*Ly*Lz*(2*m^2+1)];
end
[comp comp(:,3)-comp(:,2)]
