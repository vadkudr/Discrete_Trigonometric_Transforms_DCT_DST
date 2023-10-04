function x=TschebyshevU(n);

if n==0
    x=1; return
elseif n==1
    x=[2 0]; return
else
    a=1; b=[2 0];
    for i=2:n
        x=2*[b 0]-[zeros(1,length(b)+1-length(a)) a];
        a=b; b=x;
    end
end