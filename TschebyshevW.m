function x=TschebyshevW(n);

if n==0
    x=1; return
else
    a=TschebyshevU(n-1);
    b=TschebyshevU(n);
    x=b+[zeros(1, length(b)-length(a)) a];
end