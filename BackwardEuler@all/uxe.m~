function [uxv,uyv]=uxe(z,t)
x=z(1); y=z(2);
if(x<1)
uxv = exp(sin(t))*pi*sin(pi*y)*cos(pi*x);
uyv = exp(sin(t))*pi*sin(pi*x)*cos(pi*y);
else
uxv = -exp(sin(t))*2*pi*sin(pi*y)*cos(2*pi*x);
uyv = -exp(sin(t))*pi*sin(2*pi*x)*cos(pi*y);
end
 %uxv = -2*x*sin(20*pi*t)*exp(-10*(x^2 + y^2));
 %uyv = -2*y*sin(20*pi*t)*exp(-10*(x^2 + y^2));
%  1/10sin(20*pi*t)exp(-10*(x^2 + y^2))
%uxv = pi*exp(-2*pi^2*t).*cos(pi*x).*sin(pi*y);
%uyv = pi*exp(-2*pi^2*t).*sin(pi*x).*cos(pi*y);
