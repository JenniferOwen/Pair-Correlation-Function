comp=[];
Ly=25;
Lx=30;
for m=1:25
    NORM=PCF_normalisation_triang_zero_flux_rec_domain(Lx,Ly);
%     if m==1
%     comp(Ly/2,:)=[Ly NORM(m) 3*m*Lx*Ly/2-Ly*m-(1.5*m-1)*Lx];    
%     elseif m/2==floor(m/2) %even
%     comp(Ly/2,:)=[Ly NORM(m) 3*m*Lx*Ly/2-Ly*(m+(m-1)*m-(m/2-1)*m/2)-(1.5*m-1)*Lx];
%     else
%     comp(Ly/2,:)=[Ly NORM(m) 3*m*Lx*Ly/2-Ly*(m+(m-1)*m-ceil(m/2)*floor(m/2)+1+2)-(1.5*m-1)*Lx ];
%     end
if m==1
    comp(m,:)=[1 NORM(1) 3*Lx*Ly/2-Lx/2-Ly];
elseif m==2
    comp(m,:)=[2 NORM(2) 6*Lx*Ly/2-2*Lx-4*Ly+2];
    
else
    %comp(Ly,:)=[Ly NORM(m) 3*m*Lx*Ly/2-m^2/2*Lx-((m-1)^2+(m+1))*Ly-2*Ly*(floor((m-6)/4)+1)*floor((m-6)/4)-Ly*mod(m-6,4)*floor((m-3)/4)];
    
     a=floor((m-6)/4);
     b=floor((m-3)/4);
     c=floor((m-7)/4);
     polish=3*m*Lx*Ly/2-m^2/2*Lx+Ly*(2*a*(a-2*b+1)+b*(m-6)-m^2+m-2)+1/3*(m-1)*(m^2-2*m+6)-1/3*(c+1)*(20*c^2+37*c+12)-(m-7-4*c)*(c+1)*(m+c-2);
    % corner=1/3*(m-1)*(m^2-2*m+6)-1/3*(c+1)*(20*c^2+37*c+12)-(m-7-4*c)*(c+1)*(m+c-2);
     comp(m,:)=[m NORM(m) 3*m*Lx*Ly/2-m^2/2*Lx+Ly*(2*a*(a-2*b+1)+b*(m-6)-m^2+m-2)+1/3*(m-1)*(m^2-2*m+6)-1/3*(c+1)*(20*c^2+37*c+12)-(m-7-4*c)*(c+1)*(m+c-2)];
    comp(m,:)=[m NORM(m) polish];
    
end
end
G=[comp ((comp(:,end)-comp(:,end-1)))]
GG=G(2:end,end)-G(1:end-1,end);
GGG=GG(2:end,end)-GG(1:end-1,end);