for m=1:20
    
    M=max(-1,floor((m-8)/4));
    K=mod(max(0,m-8),4);
    ser(m)=1/3*(m-2)*((m-1)^2-2*(m-1)+6)-2*(10/6*M*(M+1)*(2*M+1)-(M+1)*M+1/2*(5*(M+1)^2-M-1)+3*(M+1)*M+1/2*(5*(M+1)^2-M-1)*K+(M+1)*(K+1)*K/2);
    
end