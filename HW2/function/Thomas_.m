function x = Thomas(A, d)
n = size(A,1);
a = [0; diag(A,-1)];
b = diag(A);
c = diag(A,1);
x = zeros(n,1);

for i = 2:n
    b(i) = b(i)-a(i)*c(i-1)/b(i-1);
    d(i) = d(i)-a(i)*d(i-1)/b(i-1);
end
    
x(n) = d(n)/b(n);

for i = n-1:-1:1
    x(i) = (d(i)-c(i)*x(i+1))/b(i);
end