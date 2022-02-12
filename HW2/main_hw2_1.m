clc; clear; close all;
addpath('function/');

G = @(x, c) (heaviside(x)-heaviside(x-c)).*(c-1).*x + ...
    + (heaviside(x-c)-heaviside(x-1)).*c.*(x-1);
alpha = 0;  beta = 0;  c = 0.31;

x0 = 0; xm_1 = 1;
mList = [15 ,31, 63, 127];
mListLength = length(mList);
ErrorList = zeros(mListLength,1);
RatioList = zeros(mListLength-1,1);

figure(1);
for i = 1:mListLength
    m = mList(i);
    h = (xm_1-x0)/(m+1);  x = x0:h:xm_1;
    f = delta_eps(h, c);
 
    U = FDM(m, x, f, alpha, beta);  % numerical solution
    U_hat = G(x', c);  % exact solution
    ErrorList(i) = max(abs(U-U_hat));  % infinity norm
    
    % Plot the numerical solution and exact solution
    subplot(2,mListLength/2,i);
    plot(x, U, 'bo', x, U_hat, 'r.');
    title(['$h=1/$', int2str(m+1)], 'interpreter', 'latex');
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('$u(x)$', 'interpreter', 'latex');
    legend('Numerical', 'Exact', 'Location', 'best');
    
    if i > 1
        RatioList(i) = log2(ErrorList(i-1)/ErrorList(i));
    end
    
end
hold off;
%%
% Plot error with loglog
figure(2);
PlotLogError(mList, ErrorList)

