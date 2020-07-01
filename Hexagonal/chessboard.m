A=zeros(45,55);
for i1=3:10:43
    for i2=3:10:53
        A(i1-2:i1+2,i2-2:i2+2)=ones(5,5);
    end
end
for i1=8:10:43
    for i2=8:10:58
        A(i1-2:i1+2,i2-2:i2+2)=ones(5,5);
    end
    
end