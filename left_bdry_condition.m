function b=left_bdry_condition(a,k)

[n,~]=size(a);
b=a;

switch k,
    case 1, % even symmetry relative to a(0)
        b(2:n,1)=a(2:n,1)+a(1,2:n)'; b(2:n-1,2)=a(2:n-1,2)+a(1,3:n)';
    case 2, % odd symmetry relative to a(-1)
        b(1:n-2,1)=a(1:n-2,1)-a(1,3:n)';
    case 3, % even symmetry relative to a(-1/2)
        b(1:n-1,1)=a(1:n-1,1)+a(1,2:n)'; b(1:n-2,2)=a(1:n-2,2)+a(1,3:n)';
    case 4, % odd symmetry relative to a(-1/2)
        b(1:n-1,1)=a(1:n-1,1)-a(1,2:n)'; b(1:n-2,2)=a(1:n-2,2)-a(1,3:n)';
end
