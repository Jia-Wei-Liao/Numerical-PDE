function G = GreenFunc(x,c)
G = (heaviside(x)-heaviside(x-c)).*(c-1).*x + ...
    + (heaviside(x-c)-heaviside(x-1)).*c.*(x-1);
end