function [uxv,uyv]=uxe(z,t)
x=z(1); y=z(2);
if(x<1)
uxv = exp(sin(t))*pi*sin(pi*y)*cos(pi*x);
uyv = exp(sin(t))*pi*sin(pi*x)*cos(pi*y);
else
uxv = -exp(sin(t))*2*pi*sin(pi*y)*cos(2*pi*x);
uyv = -exp(sin(t))*pi*sin(2*pi*x)*cos(pi*y);
end

