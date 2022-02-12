function U = FDM(m, x, f, alpha, beta, method)
h = (x(end)-x(1))/(m+1);

U = ones(m+2,1);
U(1) = alpha;  U(end) = beta;

A = diag(-2*ones(m,1)) + diag(ones(m-1,1), 1) + diag(ones(m-1,1), -1);
Ah = A/h^2;

F = [0; f(x(2:end-1))'; 0];

% 4th order convergence
if nargin == 6
    F(2:end-1) = (F(1:end-2) + 10*F(2:end-1) + F(3:end))/12;
    F(2) = F(2)+alpha/h^2;
    F(end-1) = F(end-1)-beta/h^2;
end

F = F(2:end-1);
U(2:end-1) = Thomas(Ah, F);

end