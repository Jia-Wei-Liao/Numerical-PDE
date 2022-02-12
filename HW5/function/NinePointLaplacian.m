function U = NinePointLaplacian(M, N, F, U)
dx = 1/(M+1); dy = 1/(N+1);
Lambda = TridiagED(-20, 4, 4, M);
Gamma  = TridiagED(4, 1, 1, M);

Fbar = FivePointF(F, U, dx, dy);
b = Fbar*6*dx*dy;

%bbar = b*Q;
bbar = zeros(N,M);
for i=1:M
    bbar(i, :) = dst(b(i,:))*sqrt(2/(N+1));
end

Ubar = zeros(N,M);
for k = 1:M
    A = diag(Lambda(k)*ones(N,1)) + diag(Gamma(k)*ones(N-1,1), 1) + diag(Gamma(k)*ones(N-1,1), -1);
    Ubar(:,k) = Thomas(A, bbar(:,k));
end

%U(2:end-1, 2:end-1) = Ubar*Q;
for i=1:M
    U(1+i, 2:end-1) = dst(Ubar(i,:))*sqrt(2/(N+1));
end

end