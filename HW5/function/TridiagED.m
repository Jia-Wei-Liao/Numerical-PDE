function [Lambda, Q] = TridiagED(a, b, c, N)
dx = 1/(N+1);
[Qx, Qy] = meshgrid(1:N,1:N);
Q = sin(Qx.*Qy*pi*dx)*sqrt(2/(N+1));
Lambda = a+2*sqrt(b*c)*cos((1:N)'*pi*dx);

end