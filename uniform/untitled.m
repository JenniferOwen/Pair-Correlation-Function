sum=0;
m=5;
for i=1:m-1
    sum=sum+i^2;
end
for i=1:m-2
    sum=sum+i*(2*m-i);
end