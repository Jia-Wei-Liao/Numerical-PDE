clc; clear; close all;
addpath('function/');
format long e

f1 =@(x) sin(pi*x);
u1 =@(x) -sin(pi*x)/pi^2;
alpha1 = 0; beta1 = 0;

f2 =@(x) 0*x+1;
u2 =@(x) x.^2/2;
alpha2 = 0; beta2 = 1/2;

f3 =@(x) (heaviside(x-1/2)-heaviside(x-1)).*x.^2 ...
         + (heaviside(x)-heaviside(x-1/2)).*x/2;
u3 =@(x) (heaviside(x-1/2)-heaviside(x-1)).*(x.^4/12-5*x/64-1/192) ...
         +(heaviside(x)-heaviside(x-1/2)).*(x.^3/12 -19*x/192);
alpha3 = 0; beta3 = 0;
     
%% Exercise 1
close all;

NList = [8 ,16, 32, 64];
ErrorList = zeros(4,1);
RatioList = zeros(3,1);

figure(1);
for i = 1:length(NList)
    nu = NList(i);
    x0 = 0;  xn = 1;  dx = (xn-x0)/nu;  x = x0:dx:xn;
    
    U = FiniteDifference(nu, x, f1, alpha1, beta1);  % numerical solution
    U_hat = u1(x');  % exact solution
    ErrorList(i) = max(abs(U-U_hat));  % infinity norm
    
    subplot(2,2,i);
    plot(x, U, 'r*', x, U_hat, 'bo');
    title(['h=1/', int2str(nu)]);
    xlabel('x');  ylabel('u(x)');
    legend('Numerical', 'Exact', 'Location', 'N');
    hold on;
    
    if i > 1
        RatioList(i) = log2(ErrorList(i-1)/ErrorList(i));
    end
    
end
hold off;

figure(2);
loglog(NList, ErrorList, '-ro');
axis([min(NList), max(NList), ...
    min(ErrorList), max(ErrorList)]);
xlabel('log(h)');
ylabel('log(error)');

%% Exercise 2
close all;

NList = [8 ,16, 32, 64];
ErrorList = zeros(1,4);
RatioList = zeros(1,3);

figure(1);
for i = 1:length(NList)
    nu = NList(i);
    x0 = 0;  xn = 1;  dx = (xn-x0)/nu;  x = x0:dx:xn;
    
    U = FiniteDifference(nu, x, f2, alpha2, beta2);  % numerical solution
    U_hat = u2(x');  % exact solution
    ErrorList(i) = max(abs(U-U_hat));  % infinity norm
    
    subplot(2,2,i)
    plot(x, U, 'r*', x, U_hat, 'bo');
    title(['h=1/', int2str(nu)]);
    xlabel('x');  ylabel('u(x)');
    legend('Numerical', 'Exact', 'Location', 'NW');
    hold on;
    
    if i > 1
        RatioList(i) = log2(ErrorList(i-1)/ErrorList(i));
    end
    
end

figure(2);
loglog(NList, ErrorList, '-ro');
axis([min(NList), max(NList), ...
    min(ErrorList), max(ErrorList)]);
xlabel('log(h)');
ylabel('log(error)');

%% Exercise 3
close all;

NList = [8 ,16, 32, 64];
ErrorList = zeros(1,4);
RatioList = zeros(1,3);

figure(1);
for i = 1:length(NList)
    nu = NList(i);
    x0 = 0;  xn = 1;  dx = (xn-x0)/nu;  x = x0:dx:xn;
    
    U = FiniteDifference(nu, x, f3, alpha3, beta3);  % numerical solution
    U_hat = u3(x');  % exact solution
    ErrorList(i) = max(abs(U-U_hat));  % infinity norm
    
    subplot(2,2,i)
    plot(x, U, 'r*', x, U_hat, 'bo');
    title(['h=1/', int2str(nu)]);
    xlabel('x');  ylabel('u(x)');
    legend('Numerical', 'Exact', 'Location', 'N');
    hold on;
    
    if i > 1
        RatioList(i) = log2(ErrorList(i-1)/ErrorList(i));
    end
    
end

figure(2);
loglog(NList, ErrorList, '-ro');
axis([min(NList), max(NList), ...
    min(ErrorList), max(ErrorList)]);
xlabel('log(h)');
ylabel('log(error)');



