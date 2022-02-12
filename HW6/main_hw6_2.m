clc; clear; close all;
addpath('function/');

c = 1/3;  % center

mList = [15 ,31, 63];
mListLength = length(mList);
ErrorList = zeros(mListLength,4);
RatioList = zeros(mListLength-1,4);

%% FDM
for i = 1:3
    m = mList(i);
    x0 = 0;     xm_1 = 1;
    alpha = 0;  beta = 0;
    h = (xm_1-x0)/(m+1);  x = x0:h:xm_1;
    
    % Exact solution
    Uhat = ((x.^2)/2)' + GreenFunc(x', c) - x'/2;
    
    for j =1:4
        % Numerical solution
        d = delta(h, j, c);
        f =@(x) d(x)+1;
        
        U = FDM(m, x, f, alpha, beta);
        
        ErrorList(i,j) = max(abs(U-Uhat));
        if i>1
            RatioList(i-1,j) = log2(ErrorList(i-1,j)/ErrorList(i,j));
        end
        
        subplot(4,3,(j-1)*3+i)
        plot(x, U, 'b*', x, Uhat, 'ro');
        daspect([10 5 1])
        title(['$d_h($' int2str(j) '$)$ ' 'with $h=1/$' int2str(m+1)], 'interpreter', 'latex');
        xlabel('$x$', 'interpreter', 'latex');
        ylabel('$u(x)$', 'interpreter', 'latex');
    end
end
