function VolumeForce = f(x,t)
if x(1)<1
	VolumeForce = exp(sin(t))*sin(pi*x(1))*sin(pi*x(2))*(cos(t) + 2*pi*pi);
else
	VolumeForce = -exp(sin(t))*sin(2*pi*x(1))*sin(pi*x(2))*(cos(t) + 5*pi*pi*0.5);
end
