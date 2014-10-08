clc;
clear all;
%c4n = [-1,-1;1,-1;1,1;-1,1;0,0];
%n4e = [1 2 5;2 3 5; 3 4 5; 4 1 5]; 
%n4sDb = [1 2;2 3;3 4;1 4];
%triplot(n4e,c4n(:,1),c4n(:,2));


% c4n = [0 0;2 0;2 1;0 1;0.5 0.5; 1.5 1.5;1 0;1 1]
% n4e = [1 7 5; 7 8 5; 8 4 5; 5 4 1; 7 2 6; 2 3 6; 3 8 6; 6 8 7]

 c4n = [0 0;1 0;2 0;2 1;1 1;0 1;0.5 0.5;1.5 0.5];
 n4sDb = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];
 n4e = [1 2 7;2 5 7;5 6 7;6 1 7;2 3 8;3 4 8;4 5 8;5 2 8];
T = .1;  
c = .1;
NoIteration = 4;
dt = zeros(NoIteration,1);
%dt(1) = 0.04;
N = zeros(NoIteration,1);
h = zeros(NoIteration,1);
for k = 1:NoIteration
    h(k) = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2)
    dt(k) = c*h(k)^2;
    N(k) = T/dt(k)
    figure
    [U1,A,B,ndof] = FEMPARABOLIC(c4n,n4e,unique(n4sDb),N(k,1),dt(k,1));
    %ndof
	%n4sDb
	[U2,A,B,ndof] = second(c4n,n4e,unique(n4sDb),N(k,1),dt(k,1),U1);
    [maxL2e(k),L2H1e(k)] = ocparafem(N(k,1),c4n,n4e,n4sDb,T,dt(k,1),U2)
    n4sMarked = markUniform(n4e);
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,[],n4sMarked);
end
 for Iter = 1:(k-1)
      ocmaxL2(Iter) = log(maxL2e(Iter)/maxL2e(Iter+1))/log(h(Iter)/h(Iter+1));
      ocL2H1(Iter) = log(L2H1e(Iter)/L2H1e(Iter+1))/log(h(Iter)/h(Iter+1));
  end
 ocmaxL2
 ocL2H1
