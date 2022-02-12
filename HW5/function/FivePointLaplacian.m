function U = FivePointLaplacian(M, N, F, U)
dx = 1/(M+1); dy = 1/(N+1);

F = F(2:end-1, 2:end-1);

F(1,:) = F(1,:) - U(1,2:end-1)/(dx*dy);
F(N,:) = F(N,:) - U(N+2,2:end-1)/(dx*dy);
F(:,1) = F(:,1) - U(2:end-1,1)/(dx*dy);
F(:,M) = F(:,M) - U(2:end-1,M+2)/(dx*dy);

Lambda = TridiagED(-4, 1, 1, M);
b = F*dx*dy;

%bbar = b*Q;
bbar = zeros(N,M);
for i=1:M
    bbar(i, :) = dst(b(i,:))*sqrt(2/(N+1));
end

Ubar = zeros(N,M);
for k = 1:M
    A = diag(Lambda(k)*ones(N,1)) + diag(ones(N-1,1), 1) + diag(ones(N-1,1), -1);
    Ubar(:,k) = Thomas(A, bbar(:,k));
end

%U(2:end-1, 2:end-1) = Ubar*Q;
for i=1:M
    U(1+i, 2:end-1) = dst(Ubar(i,:))*sqrt(2/(N+1));
end


end