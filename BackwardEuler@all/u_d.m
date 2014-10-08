function DirichletBoundaryValue = u_d(x,t)
DirichletBoundaryValue =  zeros(size(x,1),1);                                                                                                                                                                                                                                                                                                                                                               
DirichletBoundaryValue(find(x(:,1)==0)) = 0; %sin(pi*t).*exp(-10*(1 + x(find(x(:,1)== -1),2).^2));  %% Along x = -1
DirichletBoundaryValue(find(x(:,1)== 2)) = 0; %sin(pi*t).*exp(-10*(x(find(x(:,1)==1),2).^2 + 1));  %% Along x = 1
DirichletBoundaryValue(find(x(:,2)== 0)) = 0; %sin(pi*t).*exp(-10*(1 + x(find(x(:,2)==-1),1).^2));  %% Along y = -1
DirichletBoundaryValue(find(x(:,2) == 1)) = 0; %sin(pi*t).*exp(-10*(1 + x(find(x(:,2)== 1),1).^2));  %% Along y = 1

% DirichletBoundaryValue(find(x(:,1)==-1)) = .1 * sin(20*pi*t).*exp(-10*(1 + x(find(x(:,1)== -1),2).^2));  %% Along x = -1
% DirichletBoundaryValue(find(x(:,1)== 1)) = .1 * sin(20*pi*t).*exp(-10*(x(find(x(:,1)==1),2).^2 + 1));  %% Along x = 1
% DirichletBoundaryValue(find(x(:,2)== -1)) = .1 *sin(20*pi*t).*exp(-10*(1 + x(find(x(:,2)==-1),1).^2));  %% Along y = -1
% DirichletBoundaryValue(find(x(:,2) == 1)) = .1*sin(20*pi*t).*exp(-10*(1 + x(find(x(:,2)== 1),1).^2));  %% Along y = 1
