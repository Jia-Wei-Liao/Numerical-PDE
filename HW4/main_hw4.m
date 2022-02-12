clc; clear; close all;
addpath('function/');

u = @(x) -sin(pi*x)/pi^2;
f = @(x) sin(pi*x);
x0 = 0;     xm_1 = 1;
alpha = 0;  beta = 0;

%% FDM
mList = [15 ,31, 63, 127];
mListLength = length(mList);
ErrorList = zeros(mListLength,2);
RatioList = zeros(mListLength-1,2);

figure(1);
for i = 1:mListLength
    m = mList(i);
    h = (xm_1-x0)/(m+1);  x = x0:h:xm_1;
    
    % Exact solution
    U_hat = u(x');
    
    % 2nd order convergence
    U1 = FDM(m, x, f, alpha, beta);
    
    % 4th order convergence
    U2 = FDM(m, x, f, alpha, beta, 4);
    
    % Compute error with infinity norm
    ErrorList(i,:) = max(abs([U1-U_hat, U2-U_hat]));
    
    % Plot the numerical solution and exact solution
    subplot(2,mListLength/2,i);
    plot(x, U1, 'b*',x , U2, 'bo', x, U_hat, 'r.');
    title(['$h=1/$', int2str(m+1)], 'interpreter', 'latex');
    xlabel('$x$', 'interpreter', 'latex');
    ylabel('$u(x)$', 'interpreter', 'latex');
    legend('2-nd', '4-th', 'Exact', 'Location', 'north');

    if i > 1
        RatioList(i,:) = log2(ErrorList(i-1,:)./ErrorList(i,:));
    end

end
hold off;

%% Plot error with loglog
figure(2);
subplot(121);  PlotLogError(mList, ErrorList(:,1));
subplot(122);  PlotLogError(mList, ErrorList(:,2));

%%
figure(3);
PlotLogError(mList, ErrorList(:,1));
figure(4);
PlotLogError(mList, ErrorList(:,2));
