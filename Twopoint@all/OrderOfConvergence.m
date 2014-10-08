clc;
clear all;
format long
c4n = [0 0;1 0;2 0;2 1;1 1;0 1;0.5 0.5;1.5 0.5];
n4sDb = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];
n4e = [1 2 7;2 5 7;5 6 7;6 1 7;2 3 8;3 4 8;4 5 8;5 2 8];
%load c4n.mat;
%load n4e.mat;
%load n4sDb.mat;
T = .1; 
c = .1;
NoIteration = 3;
maxl2reconsesti = zeros(NoIteration,5);
l1spaceesti = zeros(NoIteration,5);
spaceestimate1 = zeros(NoIteration,5);
timereconesti1 = zeros(NoIteration,5);
timereconesti2 = zeros(NoIteration,5);
dt = zeros(NoIteration,1);
N = zeros(NoIteration,1);
h = zeros(NoIteration,1); 
for k = 1:NoIteration
    h(k) = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
    dt(k) = c*h(k);
    N(k) = T/dt(k);
	figure
    [U,A,B,ndof] = TwoPoint(c4n,n4e,unique(n4sDb),N(k,1),dt(k,1));
	for l = 1:(N(k)+1)
		spaceestimate1(k,l) = dt(k,1)*timelevelspaceestimator1(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U);
		timereconesti1(k,l) = dt(k,1)*timeleveltimerecestimator1(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U);
		timereconesti2(k,l) = max(timeleveltimerecestimator2(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U));
		maxl2reconsesti(k,l) = max(timelevell2spaceestimate(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U));
		l1spaceesti(k,l) = dt(k,1)*timelevelspaceestimator(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U);
	end
	%[maxL2e(k),L2H1e(k)] = ocparafem(N(k,1),c4n,n4e,n4sDb,T,dt(k,1),U);
    n4sMarked = markUniform(n4e);
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,[],n4sMarked);
end

for j=1:NoIteration-1
	for i=1:N(j)+1
		eocspace1(j,i) = log(spaceestimate1(j+1,2*i-1)/spaceestimate1(j,i))/log(0.5);
		eoctimerec1(j,i) = log(timereconesti1(j+1,2*i-1)/timereconesti1(j,i))/log(0.5);
		eoctimerec2(j,i) = log(timereconesti2(j+1,2*i-1)/timereconesti2(j,i))/log(0.5);
		eocl1space(j,i) = log(l1spaceesti(j+1,2*i-1)/l1spaceesti(j,i))/log(0.5);
		eocmaxl2(j,i) = log(maxl2reconsesti(j+1,2*i-1)/maxl2reconsesti(j,i))/log(0.5);
	end
end
