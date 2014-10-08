function[U, A, B, ndof] = TwoPoint(c4n,n4e,dirichlet,N,dt)
M = size(c4n,1); A = sparse(M,M); B = sparse(M,M);
x = zeros(M,1); 
	interFaceNodes = find(c4n(:,1)==1);
	FreeNodes = setdiff(1:M,unique(dirichlet));
	FreeNodes = setdiff(FreeNodes,interFaceNodes);
	ndof = length(FreeNodes);
U = sparse(M,N+1); 
ndof = length(FreeNodes);

%% Stiffnes Matrix
for j = 1:size(n4e,1)
  A(n4e(j,:),n4e(j,:)) = A(n4e(j,:),n4e(j,:))+ stima3(c4n(n4e(j,:),:));
end

%% Mass Matrix
for j = 1:size(n4e,1)
  B(n4e(j,:),n4e(j,:)) = B(n4e(j,:),n4e(j,:)) ...
     + det([1,1,1;c4n(n4e(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
end

	x = find(c4n(:,1)<1); 
	U(x,1) = sin(pi*c4n(x,1)).*sin(pi*c4n(x,2));
	x = find(c4n(:,1)>1);
	U(x,1) = -sin(2*pi*c4n(x,1)).*sin(pi*c4n(x,2));
	U(find(c4n(:,1)==1),1) = 0; 
n=2;
  b = sparse(M,1);
  for j = 1:size(n4e,1)
    b(n4e(j,:)) = b(n4e(j,:)) ...
  + 0.5*det([1,1,1; c4n(n4e(j,:),:)']) * dt*(f(sum(c4n(n4e(j,:),:))/3,n*dt)/6 + f(sum(c4n(n4e(j,:),:))/3,(n-1)*dt)/6);
  end
  %C = C + (dt^2)*A*U(:,n-2);
  b = b + (B - (dt/2)*A) * U(:,n-1);
  u = sparse(M,1);
	
  u(unique(dirichlet)) = u_d(c4n(unique(dirichlet),:),n*dt);
	u(interFaceNodes) = 0;
  b = b - (0.5*dt*A + B) * u;
  u(FreeNodes) = (0.5*dt*A(FreeNodes,FreeNodes)+ ...
      B(FreeNodes,FreeNodes))\b(FreeNodes); 
  U(:,n) = u;

for n=3:N+1
	b = sparse(M,1);
  for j = 1:size(n4e,1)
    b(n4e(j,:)) = b(n4e(j,:)) ...
  + det([1,1,1; c4n(n4e(j,:),:)']) *dt*(f(sum(c4n(n4e(j,:),:))/3,n*dt))/9;
  end
  b = b + (4*B * U(:,n-1))/3-(B*U(:,n-2))/3;
  u = sparse(M,1);
	
  u(unique(dirichlet)) = u_d(c4n(unique(dirichlet),:),n*dt);
	u(interFaceNodes) = 0;
  b = b - ((2*dt*A)/3 + B) * u;
  u(FreeNodes) = ((2/3)*dt*A(FreeNodes,FreeNodes)+ ...
      B(FreeNodes,FreeNodes))\b(FreeNodes); 
  U(:,n) = u;
end

show(n4e,[],c4n,full(U(:,N+1)));
