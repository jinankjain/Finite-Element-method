function uv=ue(z,t)
x=z(1); y=z(2);
%uv=x*(1-x)*y*(1-y);
%uv = .1*sin(20*pi*t).*exp(-10*(x^2 + y^2));
%uv = exp(-2*pi^2*t).*sin(pi*x).*sin(pi*y);
if(x<1)
uv =  exp(sin(t))*sin(pi*x)*sin(pi*y);
else
uv = -exp(sin(t))*sin(2*pi*x)*sin(pi*y);
end
