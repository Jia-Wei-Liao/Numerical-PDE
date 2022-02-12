function U = FiniteDifference(nu, x, f, alpha, beta)
h = 1/nu;
U = ones(nu+1,1);
U(1) = alpha;  U(end) = beta;
A = diag(-2*ones(nu-1,1)) + diag(ones(nu-2,1), 1) + diag(ones(nu-2,1), -1);

F = f(x(2:end-1))';
F(1) = F(1)-alpha/h^2;
F(end) = F(end)-beta/h^2;

U(2:end-1) = A\F*h^2;
%U(2:end-1) = Thomas(A, F)*h^2;

end