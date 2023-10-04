function x=TschebyshevT(n);

if n==0
    x=1; return
elseif n==1
    x=[1 0]; return
else
    a=1; b=[1 0];
    for i=2:n
        x=2*[b 0]-[zeros(1,length(b)+1-length(a)) a];
        a=b; b=x;
    end
end