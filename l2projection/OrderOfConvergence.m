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
NoIteration = 8;
dt = zeros(NoIteration,1);
N = zeros(NoIteration,1);
h = zeros(NoIteration,1); 
for k = 1:NoIteration
    h(k) = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
	%figure
    [L2f(k),ndof] = L2Projection(c4n,n4e,unique(n4sDb));
    n4sMarked = markUniform(n4e);
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,[],n4sMarked);
    %L2f - f(c4n(:,1),c4n(:,2))
end
L2f
for Iter = 1:(k-1)
      ocL2pro(Iter) = log(L2f(Iter)/L2f(Iter+1))/log(h(Iter)/h(Iter+1));
  end
 ocL2pro


