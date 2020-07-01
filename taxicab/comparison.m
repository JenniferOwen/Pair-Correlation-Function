comp=[];
Ly=30
Lx=25;
for m=1:25
    NORM=PCF_normalisation_rectangular(Lx,Ly);
    comp(m,:)=[NORM(m) 2*m*Lx*Ly-(m^2*(Lx+Ly)-1/3*(m^3-m))];
end
[comp comp(:,1)-comp(:,2)]
