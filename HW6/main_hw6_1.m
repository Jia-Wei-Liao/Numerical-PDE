clc; clear; close all;
addpath('function/');

c = 1/3;  % center
m = 15;  % mesh
x0 = 0;     xm_1 = 1;
alpha = 0;  beta = 0;
h = (xm_1-x0)/(m+1);  x = x0:h:xm_1;

ErrorList = zeros(m+2,4);

% Exact solution
Uhat = GreenFunc(x', c);

%% FDM
for ii =1:4
    % Numerical solution
    f = delta(h, ii, c);
    U = FDM(m, x, f, alpha, beta);
    ErrorList(:,ii) = abs(U-Uhat);
    
    subplot(2,2,ii)
    plot(x, U, 'ro', x, Uhat, 'b*');
    title(['$d_h($' int2str(ii) '$)$'], 'interpreter', 'latex');
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('$u(x)$', 'interpreter', 'latex');
    legend('Numerical', 'Exact', 'Location', 'best');
end
