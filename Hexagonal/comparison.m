comp=[];
Ly=12;
Lx=22;
for m=1:8
    NORM=PCF_normalisation_triang_zero_flux_rec_domain(Lx,Ly);
%     if m==1
%     comp(Ly/2,:)=[Ly NORM(m) 3*m*Lx*Ly/2-Ly*m-(1.5*m-1)*Lx];    
%     elseif m/2==floor(m/2) %even
%     comp(Ly/2,:)=[Ly NORM(m) 3*m*Lx*Ly/2-Ly*(m+(m-1)*m-(m/2-1)*m/2)-(1.5*m-1)*Lx];
%     else
%     comp(Ly/2,:)=[Ly NORM(m) 3*m*Lx*Ly/2-Ly*(m+(m-1)*m-ceil(m/2)*floor(m/2)+1+2)-(1.5*m-1)*Lx ];
%     end
if m==1
    comp(m,:)=[m NORM(m) 3*m*Lx*Ly/2-m^2/2*Lx-Ly];
else
    %comp(Ly,:)=[Ly NORM(m) 3*m*Lx*Ly/2-m^2/2*Lx-((m-1)^2+(m+1))*Ly-2*Ly*(floor((m-6)/4)+1)*floor((m-6)/4)-Ly*mod(m-6,4)*floor((m-3)/4)];
     M=max(-1,floor((m-7)/4));
    K=mod(max(0,m-7),4);
    corner=1/3*(m-1)*(m^2-2*m+6)-10/3*M*(M+1)*(2*M+1)+2*(M+1)*M-(5*(M+1)^2-M-1)-6*(M+1)*M-(5*(M+1)^2-M-1)*K-(M+1)*(K+1)*K;
     comp(m,:)=[m NORM(m) 3*m*Lx*Ly/2-m^2/2*Lx-((m-1)^2+(m+1))*Ly+Ly*(2*(floor((m-6)/4)+1)*floor((m-6)/4)+mod(m-6,4)*floor((m-3)/4))+corner];
    
end
end
G=[comp ((comp(:,end)-comp(:,end-1)))]
GG=G(2:end,end)-G(1:end-1,end);
GGG=GG(2:end,end)-GG(1:end-1,end);