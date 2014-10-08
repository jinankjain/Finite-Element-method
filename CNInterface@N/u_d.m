function DirichletBoundaryValue = u_d(x,t)
DirichletBoundaryValue =  zeros(size(x,1),1);                                                                                                                                                                                                                                                                                                                                                               
DirichletBoundaryValue(find(x(:,1)==0)) = 0; % Along x = 0
DirichletBoundaryValue(find(x(:,1)== 2)) = 0; % Along x = 2
DirichletBoundaryValue(find(x(:,2)== 0)) = 0; % Along y = 0
DirichletBoundaryValue(find(x(:,2) == 1)) = 0; % Along y = 1

