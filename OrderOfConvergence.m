clc;
clear all;
% c4n represents coordinates for node. In this portion we have declared the initial mess structure that users need to modify according to their initial mess representation.
c4n = [0 0;1 0;2 0;2 1;1 1;0 1;0.5 0.5;1.5 0.5];

% n4sDb represents nodes of dirchlet boundary condition
n4sDb = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];

% n4e represents nodes for elements each node numbering which is done according to indexing their in c4n such as (0,0) represents nodal number 1 and similarly. Using n4e we try mark down number of triangulations formed on each refinement and numbering is done always in anti-clockwise fashion  
n4e = [1 2 7;2 5 7;5 6 7;6 1 7;2 3 8;3 4 8;4 5 8;5 2 8];

% Time T upto which we want to iteartions from t = 0 to t = T
T = .1; 
% c is constant used for obtaining the order of convergence 
c = .1;
NoIteration = 4; % Number iterations can be changed by users according to need
dt = zeros(NoIteration,1); % Represents the time interval in which we want to break 0 to T
N = zeros(NoIteration,1);
h = zeros(NoIteration,1); % Mesh size parameter
for k = 1:NoIteration
    h(k) = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2)
    dt(k) = c*h(k)^2;
    N(k) = T/dt(k);
    [U,A,B,ndof] = femParabolic(c4n,n4e,unique(n4sDb),N(k,1),dt(k,1)); % This subroutine helps in caluclating in finite element solution which gives for paramater degree of U(finite element solution), ndof(number of free nodes), Mass Matrix(A) and Volume Matrix(B)
    [maxL2e(k),L2H1e(k)] = ocparafem(N(k,1),c4n,n4e,n4sDb,T,dt(k,1),U)% This is a subrountine helps in calculating in order of convergence
    n4sMarked = markUniform(n4e);
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,[],n4sMarked); % This subroutines helps in refinement of current mesh.
end
 for Iter = 1:(k-1)
      ocmaxL2(Iter) = log(maxL2e(Iter)/maxL2e(Iter+1))/log(h(Iter)/h(Iter+1));
      ocL2H1(Iter) = log(L2H1e(Iter)/L2H1e(Iter+1))/log(h(Iter)/h(Iter+1));
  end
 ocmaxL2 % L2 order of convergence
 ocL2H1 % H1 order of convergence
