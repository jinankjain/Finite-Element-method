function uv=ue(z,t)
x=z(1); y=z(2);
if(x<1)
uv =  exp(sin(t))*sin(pi*x)*sin(pi*y);
else
uv = -exp(sin(t))*sin(2*pi*x)*sin(pi*y);
end
