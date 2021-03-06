function exact(c4n,elements3,elements4,t);
u = zeros(size(c4n,1),1);
for i = 1:size(c4n,1)
	a = c4n(i,:);
	x = a(1);
	y = a(2);	
	if(a(1)<1)
		u(i) =  exp(sin(t))*(sin(pi*x))*sin(pi*y)
	else
		u(i) =  -exp(sin(t))*sin(2*pi*x)*sin(pi*y)
	end
end
trisurf(elements3,c4n(:,1),c4n(:,2),u,'facecolor','interp')
hold on
trisurf(elements4,c4n(:,1),c4n(:,2),u,'facecolor','interp')
hold off
view(10,40);
title('Exact Solution of the Problem')
