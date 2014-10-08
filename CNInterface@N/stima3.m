function M = stima3(vertices)

a = vertices(:,1);

% condition could be modified according to interface here the interface was at (1,0)
if((a(1)<1 && a(2)<1) || (a(2)<1 && a(3)<1) || (a(3)<1 && a(1)<1) || (a(1)==1&&a(2)==1&&a(3)<1)|| (a(3)==1&&a(2)==1&&a(1)<1)|| (a(1)==1&&a(3)==1&&a(2)<1))
	d = size(vertices,2);
	G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
	M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
else
	d = size(vertices,2);
	G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
	M = 0.5*det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
end
