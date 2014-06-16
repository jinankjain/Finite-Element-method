function[U, A, B, ndof] = femParabolic(c4n,n4e,dirichlet,N,dt)
	% Intialize various paramaters	
	M = size(c4n,1);
	A = sparse(M,M);
	B = sparse(M,M);
	x = zeros(M,1); 
	U = sparse(M,N+1);
	ndof = length(FreeNodes);
	
	% Segreagating FreeNodes with Interface Nodes
	interFaceNodes = find(c4n(:,1)==1);
	FreeNodes = setdiff(1:M,unique(dirichlet));
	FreeNodes = setdiff(FreeNodes,interFaceNodes);
	
	% This loops helps in calculating Stifness Matrix using one subroutine to calculate stiffness matrix
	for j = 1:size(n4e,1)
  		A(n4e(j,:),n4e(j,:)) = A(n4e(j,:),n4e(j,:))...
      		+ stima3(c4n(n4e(j,:),:));
	end
	
	% This loops helps in calculating Volume force Matrix
	for j = 1:size(n4e,1)
  		B(n4e(j,:),n4e(j,:)) = B(n4e(j,:),n4e(j,:)) ...
      		+ det([1,1,1;c4n(n4e(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
	end
	
	% Initializing U matrix at t=0
	x = find(c4n(:,1)<1); 
	U(x,1) = sin(pi*c4n(x,1)).*sin(pi*c4n(x,2));
	x = find(c4n(:,1)>1);
	U(x,1) = -sin(2*pi*c4n(x,1)).*sin(pi*c4n(x,2));
	U(find(c4n(:,1)==1),1) = 0; 
	
	% Assembly Method 
	for n = 2:N+1
  		b = sparse(M,1);
  		for j = 1:size(n4e,1)
    		b(n4e(j,:)) = b(n4e(j,:)) + det([1,1,1; c4n(n4e(j,:),:)']) * dt*f(sum(c4n(n4e(j,:),:))/3,n*dt)/6;
  		end
  	b = b + B * U(:,n-1);
  	u = sparse(size(c4n,1),1);
  	u(unique(dirichlet)) = u_d(c4n(unique(dirichlet),:),n*dt);
  	u(interFaceNodes) = 0;
  	b = b - (dt*A + B) * u;
  	u(FreeNodes) = (dt*A(FreeNodes,FreeNodes)+ ...
      	B(FreeNodes,FreeNodes))\b(FreeNodes);
  	U(:,n) = u;
	end
	show(n4e,[],c4n,full(U(:,N+1))); % Subroutine for plotting the solution
	end
