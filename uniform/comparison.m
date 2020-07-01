comp=[];
Ly=30;
Lx=25;
for m=1:24
    NORM=PCF_normalisation_rectangular_sparse_zero_flux_uniform(Lx,Ly);
    comp(m,:)=[Ly NORM(m) 4*m*Lx*Ly-m^2*(3*Lx+3*Ly)+2*m^3];%-(m^2*(Lx+Ly)-1/3*(m^3-m))];
end
[comp comp(:,3)-comp(:,2)]
