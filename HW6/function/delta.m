function  f = delta(h, idx, c)

if idx==1
    f =@(x) (heaviside(x+h)-heaviside(x-h)).*(h-abs(x))/h^2;
end

if idx==2
    f =@(x) (heaviside(x+2*h)-heaviside(x-2*h)).*(2*h-abs(x))/(4*h^2);
end

if idx==3
    f =@(x) (heaviside(x+2*h)-heaviside(x-2*h)).*(cos(pi*x/(2*h))+1)/(4*h);
end

if idx==4
    f =@(x) ((heaviside(abs(x))-heaviside(abs(x)-h)).*(1-(abs(x)/h).^2)...
    +(heaviside(abs(x)-h)-heaviside(abs(x)-2*h)).*(2 - 3*abs(x)/h + (abs(x)/h).^2))/h;
end

if nargin == 3
    f = @(x) f(x-c);
end

end