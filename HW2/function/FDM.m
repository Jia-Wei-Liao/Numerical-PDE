function U = FDM(m, x, f, alpha, beta, um_1, method)
h = (x(end)-x(1))/(m+1);

if nargin==5
    U = ones(m+2,1);
    U(1) = alpha;  U(end) = beta;
    A = diag(-2*ones(m,1)) + diag(ones(m-1,1), 1) + diag(ones(m-1,1), -1);

    F = f(x(2:end-1))';
    F(1) = F(1)-alpha/h^2;
    F(end) = F(end)-beta/h^2;

    %U(2:end-1) = A\F*h^2;
    U(2:end-1) = Thomas(A, F)*h^2;

else
    F = f(x)';

    if strcmp(method, 'OneSideDiff')  % First order convergence
        A = diag([-h; -2*ones(m,1); h]) ...
            + diag([ones(m,1); -h], -1) ...
            + diag([h; ones(m,1)], 1);
        F(1) = alpha;  F(end) = beta;

    elseif strcmp(method, 'CenterDiff')  % Second order convergence
        A = diag([-1; -2*ones(m,1); -1]) ...
            + diag(ones(m+1,1), -1) ...
            + diag(ones(m+1,1), 1);
        F(1) = F(1)/2 + alpha/h;
        F(end) = F(end)/2 - beta/h;

    end

    Ah = A/h^2;
    U = Thomas(Ah, F, um_1);

end
end
