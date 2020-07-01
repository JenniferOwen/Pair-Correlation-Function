function [PCF ] = PCF_lattice_unif(D)
%Function calculates the PCF for a given rectangular domain D where D is 
%size Lx by Ly

agt=nnz(D); 
%Number of agents on lattice

% Size domain
Lx=size(D,2);
Ly=size(D,1);

% Coordinates agents
C=zeros(agt,2);

[C(:,1),C(:,2)]=find(D);

% Initialise PCF
PCF=zeros(1,min(Lx,Ly));

for r=1:min(Lx,Ly)
    for i=1:agt-1
        for j=i+1:agt
            %Calculate the number of agents of distance r from each other
            if abs(C(i,1)-C(j,1))+abs(C(i,2)-C(j,2))==r
                PCF(r)=PCF(r)+1;
            end %if
        end %for
    end %for 
end %for

PCF_NORM=1:min(Lx,Ly);
%Normalise
PCF_NORM=4*PCF_NORM*Lx*Ly-(PCF_NORM.^2*(3*Lx+3*Ly)+2*PCF_NORM.^3);
rho=agt/(Lx*Ly); %Choose one element on lattice
rho_hat=(agt-1)/(Lx*Ly-1); %Choose second element on lattice

PCF=PCF./(PCF_NORM*rho*rho_hat);

subplot(1,2,1)
spy(D)
title('D')
subplot(1,2,2)
plot(PCF(1:min(Lx,Ly)))
axis([0,min(Lx,Ly),0,2])
title('PCF(D)')
end

