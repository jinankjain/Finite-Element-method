clc;
clear all;
%c4n = [0 0;1 0;2 0;2 1;1 1;0 1;0.5 0.5;1.5 0.5];
%n4sDb = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];
%n4e = [1 2 7;2 5 7;5 6 7;6 1 7;2 3 8;3 4 8;4 5 8;5 2 8];
%triplot(n4e,c4n(:,1),c4n(:,2));
format long;
T = .1;  
c = .1;
load c4n.mat;
load n4e.mat;
load n4sDb.mat;
NoIteration = 4;
maxl2reconsesti = zeros(NoIteration,1);
l1timeesti = zeros(NoIteration,1);
l2timeh1spaceesti = zeros(NoIteration,1);
l1spaceesti = zeros(NoIteration,1);
dt = zeros(NoIteration,1);
N = zeros(NoIteration,1);
h = zeros(NoIteration,1);
for k = 1:NoIteration
    h(k) = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
    dt(k) = c*h(k)^2;
    N(k) = T/dt(k);
    figure
    [U,A,B,ndof] = FEMPARABOLIC(c4n,n4e,unique(n4sDb),N(k,1),dt(k,1));
		l1spaceesti(k) = dt(k,1)*timelevelspaceestimator(c4n,n4e,n4sDb,T,N(k,1),dt(k,1),h(k),U);
		maxl2reconsesti(k) = max(timelevell2spaceestimate(c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U));
		l1timeesti(k) = dt(k,1)*timeleveltimeestimator(c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U);
		l2timeh1spaceesti(k) = sqrt(dt(k,1)*timelevelH1spaceestimate(c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U,k));
    %[maxL2e(k),L2H1e(k)] = ocparafem(N(k,1),c4n,n4e,n4sDb,T,dt(k,1),U);
    n4sMarked = markUniform(n4e);
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,[],n4sMarked);
end
%for z = 1:k-1
%	for l=1:(N(k)+1)
%		eocmaxl2(z,l) = log(maxl2reconsesti(z+1,l)/maxl2reconsesti(z,l))/log(h(z+1)/h(z))
%    	eocl1time(z,l) = log(l1timeesti(z+1,l)/l1timeesti(z,l))/log(h(z+1,1)/h(z,1))
%		eocl2h1(z,l) = log(l2timeh1spaceesti(z+1,l)/l2timeh1spaceesti(z,l))/log(h(z+1,1)/h(z,1))
%		eocl1space(z,l) = log(l1spaceesti(z+1,l)/l1spaceesti(z,l))/log(h(z+1,1)/h(z,1))
%	end
%end
%eocl1time
