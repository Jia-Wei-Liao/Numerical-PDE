clc; clear; close all;
addpath('function/');

%% Test 1
u = @(x,y) x.*(1-x).*y.*(1-y);
f = @(x, y) -2*y.*(1-y)-2*x.*(1-x);

%% Test 2
u = @(x,y) exp(x).*sin(y);
f = @(x, y) 0 + x.*0 + y.*0;

%% Exercise 1
u = @(x,y) x.^2+y.^2;
f = @(x, y) 4 + x.*0 + y.*0;

%%  Exercise 2
u = @(x, y) exp(1).^(x+y);
f = @(x, y) 2*exp(1).^(x+y);

%% FDM
figure(1)
mList = [7, 15 ,31, 63];
mListLength = length(mList);
ErrorList = zeros(mListLength,2);
RatioList = zeros(mListLength-1,2);

for ii = 1:mListLength
    M = mList(ii);  N = mList(ii);
    dx = 1/(M+1); dy = 1/(N+1);
    x0 = 0; xM_1 = 1; x = x0:dx:xM_1;
    y0 = 0; yN_1 = 1; y = y0:dy:yN_1;
    [X, Y] = meshgrid(x,y);
    F = f(X,Y);
    
    % Exact solution
    Uhat = u(X,Y);
    
    % Boundary condiction
    U = zeros(N+2, M+2);
    U(1,:) = u(X(1,:), Y(1,:));
    U(N+2,:) = u(X(N+2,:), Y(N+2,:));
    U(:,1) = u(X(:,1), Y(:,1));
    U(:,M+2) = u(X(:,M+2), Y(:,M+2));
    
    % Five-point Laplacian method
    U1 = FivePointLaplacian(M, N, F, U);
    
    % Nine-point Laplacian method
    U2 = NinePointLaplacian(M, N, F, U);
    
    % Compute error
    ErrorList(ii,1) = max(abs(U1-Uhat), [], 'all');
    ErrorList(ii,2) = max(abs(U2-Uhat), [], 'all');
    if ii > 1
        RatioList(ii,:) = log2(ErrorList(ii-1,:)./ErrorList(ii,:));
    end
    
    subplot(4, 3, 1+3*(ii-1));
    surf(X,Y,Uhat, 'edgecolor', 'none');
    view([0,90]);
    title('exact', 'interpreter', 'latex');
    colorbar; axis square;

    subplot(4, 3, 2+3*(ii-1));
    surf(X,Y,U1, 'edgecolor', 'none');
    view([0,90]);
    title('5-point', 'interpreter', 'latex');
    colorbar; axis square;

    subplot(4, 3, 3+3*(ii-1));
    surf(X,Y,U2, 'edgecolor', 'none');
    view([0,90]);
    title('9-point', 'interpreter', 'latex');
    colorbar; axis square;
end

%%
figure(2);
PlotLogError(mList, ErrorList(:,1))
title('Five-point');

figure(3);
PlotLogError(mList, ErrorList(:,2))
title('Nine-point');

