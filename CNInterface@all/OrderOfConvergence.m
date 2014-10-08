clc;
clear all;
format long
%c4n = [0 0;1 0;2 0;2 1;1 1;0 1;0.5 0.5;1.5 0.5];

 %n4sDb represents nodes of dirchlet boundary condition
%n4sDb = [1 2; 2 3; 3 4; 4 5; 5 6; 6 1];

% n4e represents nodes for elements each node numbering which is done according to indexing their in c4n such as (0,0) represents nodal number 1 and similarly. Using n4e we try mark down number of triangulations formed on each refinement and numbering is done always in anti-clockwise fashion  
%n4e = [1 2 7;2 5 7;5 6 7;6 1 7;2 3 8;3 4 8;4 5 8;5 2 8];
load c4n.mat;
load n4e.mat;
load n4sDb.mat;
% Time T upto which we want to iteartions from t = 0 to t = T

T = .1; 
% c is constant used for obtaining the order of convergence 
c = .1;
NoIteration = 3; % Number iterations can be changed by users according to need
maxl2reconsesti = zeros(NoIteration,33);
l1spaceesti = zeros(NoIteration,33);
spaceestimate1 = zeros(NoIteration,33);
timereconesti1 = zeros(NoIteration,33);
timereconesti2 = zeros(NoIteration,33);
dt = zeros(NoIteration,1); % Represents the time interval in which we want to break 0 to T
N = zeros(NoIteration,1);
h = zeros(NoIteration,1); % Mesh size parameter
for k = 1:NoIteration
    h(k) = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
    dt(k) = c*h(k)
    N(k) = T/dt(k)
	figure
    [U,A,B,ndof] = CNFEMInterface(c4n,n4e,unique(n4sDb),N(k,1),dt(k,1)); 
	%WH1Norm = Wvalue(c4n,n4e,n4sDb,T,N(k,1),dt(k,1))  %%% Maxvalue = 55.5
	%WL2Norm = Wvalue1(c4n,n4e,n4sDb,T,N(k,1),dt(k,1)) %%% Maxvalue = 9.1
    %GNorm = g(c4n,n4e,U,T,dt(k,1),N(k,1))%% Max value 744.5
    %[maxL2e(k),L2H1e(k)] = ocparafem(N(k,1),c4n,n4e,n4sDb,T,dt(k,1),U);% This is a subrountine helps in calculating in order of convergence

for l = 1:(N(k)+1)
	spaceestimate1(k,l) = dt(k,1)*timelevelspaceestimator1(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),h(k),U)
	timereconesti1(k,l) = dt(k,1)*timeleveltimerecestimator1(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),h(k),U)
	timereconesti2(k,l) = max(timeleveltimerecestimator2(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),h(k),U))
	maxl2reconsesti(k,l) = max(timelevell2spaceestimate(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),U))
	l1spaceesti(k,l) = dt(k,1)*timelevelspaceestimator(l,c4n,n4e,n4sDb,T,N(k,1),dt(k,1),h(k),U)
end
    n4sMarked = markUniform(n4e);
    [c4n,n4e,n4sDb,n4sNb] = refineRGB(c4n,n4e,n4sDb,[],n4sMarked); % This subroutines helps in refinement of current mesh.
end
 %for z = 1:k-1
	  %eocspace1(z) = log(spaceestimate1(z+1,1)/spaceestimate1(z,1))/log(h(z+1)/h(z));
	  %eoctimerec1(z) = log(timereconesti1(z+1,1)/timereconesti1(z,1))/log(h(z+1,1)/h(z,1));
      %eoctimerec2(z) = log(timereconesti2(z+1,1)/timereconesti2(z,1))/log(h(z+1,1)/h(z,1));
      %ocL2H1(Iter) = log(L2H1e(Iter)/L2H1e(Iter+1))/log(h(Iter)/h(Iter+1));
  %end
%%%%% For tau proportional to h
for j=1:NoIteration-1
	for i=1:N(j)+1
		eocspace1(j,i) = log(spaceestimate1(j+1,2*i-1)/spaceestimate1(j,i))/log(0.5);
		eoctimerec1(j,i) = log(timereconesti1(j+1,2*i-1)/timereconesti1(j,i))/log(0.5);
		eoctimerec2(j,i) = log(timereconesti2(j+1,2*i-1)/timereconesti2(j,i))/log(0.5);
		eocl1space(j,i) = log(timelevelspaceestimator(j+1,2*i-1)/timelevelspaceestimator(j,i))/log(0.5);
		eocmaxl2(j,i) = log(timelevell2spaceestimate(j+1,2*i-1)/timelevell2spaceestimate(j,i))/log(0.5);
end
end
%%%%% For tau proportional to h*h
for j=1:NoIteration-1
	for i=1:N(j)+1
		eocl1time(j,i) = log(l1timeesti(j+1,4*i-3)/l1timeesti(j,i))/log(0.5);
		eocl2h1(j,i) = log(l2timeh1spaceesti(j+1,4*i-3)/l2timeh1spaceesti(j,i))/log(0.5);
		eocl1space(j,i) = log(timelevelspaceestimator(j+1,4*i-3)/timelevelspaceestimator(j,i))/log(0.5);
		eocmaxl2(j,i) = log(timelevell2spaceestimate(j+1,4*i-3)/timelevell2spaceestimate(j,i))/log(0.5);
end
end		
%eocspace1
%eoctimerec1
%eoctimerec2
%eocl1space
%eocmaxl2
