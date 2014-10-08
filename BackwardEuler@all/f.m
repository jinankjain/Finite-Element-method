function VolumeForce = f(x,t)
%VolumeForce = 2*pi*cos(20*pi*t).*exp(-10*(x(1)^2 + x(2)^2)) - 2*sin(20*pi*t).*exp(-10*(x(1)^2 + x(2)^2)).*(20*(x(1)^2 + x(2)^2) - 2);
%VolumeForce = pi*cos(pi*t).*exp(-10*(x(1)^2 + x(2)^2)) - 20*sin(pi*t).*exp(-10*(x(1)^2 + x(2)^2)).*(20*(x(1)^2 + x(2)^2) - 2);
if x(1)<1
	VolumeForce = exp(sin(t))*sin(pi*x(1))*sin(pi*x(2))*(cos(t) + 2*pi*pi);
else
	VolumeForce = -exp(sin(t))*sin(2*pi*x(1))*sin(pi*x(2))*(cos(t) + 5*pi*pi*0.5);
end