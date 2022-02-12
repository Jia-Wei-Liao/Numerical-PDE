function  f = delta_eps(eps, c)
f =@(x) (heaviside(x+eps)-heaviside(x)).*(eps+x)/eps^2 ...
    + (heaviside(x)-heaviside(x-eps)).*(eps-x)/eps^2;

if nargin == 2
    f = @(x) f(x-c);
end

end