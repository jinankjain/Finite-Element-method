function[U, A, B, ndof] = CNFEMInterface(c4n,n4e,dirichlet,N,dt)
M = size(c4n,1); A = sparse(M,M); B = sparse(M,M); C = sparse(M,1);
FreeNodes = setdiff(1:M,unique(dirichlet)); U = sparse(M,N+1); 
ndof = length(FreeNodes);
for j = 1:size(n4e,1)
  A(n4e(j,:),n4e(j,:)) = A(n4e(j,:),n4e(j,:))+ stima3(c4n(n4e(j,:),:));
end
for j = 1:size(n4e,1)
  B(n4e(j,:),n4e(j,:)) = B(n4e(j,:),n4e(j,:)) ...
     + det([1,1,1;c4n(n4e(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
end
U(:,1) = sin(c4n(:,1)).*sin(c4n(:,2)); 
for n = 2
  b = sparse(M,1);
  for j = 1:size(n4e,1)
    b(n4e(j,:)) = b(n4e(j,:)) ...
  + 0.5*det([1,1,1; c4n(n4e(j,:),:)'])*dt*(f(sum(c4n(n4e(j,:),:))/3,n*dt)/6 + f(sum(c4n(n4e(j,:),:))/3,(n-1)*dt)/6);
  end
  b = b + (B - (dt/2)*A + (1/4)*dt^2*A) * U(:,n-1);
  u = sparse(M,1);
  u(unique(dirichlet)) = u_d(c4n(unique(dirichlet),:),n*dt);
  b = b - (0.5*dt*A + B - 0.25*dt^2*A) * u;
  u(FreeNodes) = (0.5*dt*A(FreeNodes,FreeNodes)+ ...
      B(FreeNodes,FreeNodes) - 0.25*dt^2*A(FreeNodes,FreeNodes))\b(FreeNodes); 
  U(:,n) = u;
end

for n = 3
  b = sparse(M,1);
  for j = 1:size(n4e,1)
    b(n4e(j,:)) = b(n4e(j,:)) ...
  + 0.5*det([1,1,1; c4n(n4e(j,:),:)']) * dt*(f(sum(c4n(n4e(j,:),:))/3,n*dt)/6 + f(sum(c4n(n4e(j,:),:))/3,(n-1)*dt)/6);
  end
  b = b + (B - (dt/2)*A + (3/4)*dt^2*A) * U(:,n-1) + 0.5*(dt^2)*A*U(:,1);
  u = sparse(M,1);
  u(unique(dirichlet)) = u_d(c4n(unique(dirichlet),:),n*dt);
  b = b - (0.5*dt*A + B - 0.25*dt^2*A) * u;
  u(FreeNodes) = (0.5*dt*A(FreeNodes,FreeNodes)+ ...
      B(FreeNodes,FreeNodes) - 0.25*dt^2*A(FreeNodes,FreeNodes))\b(FreeNodes); 
  U(:,n) = u;
end
C1 = 0.5*(dt^2)*A*U(:,1);
C = C1;
for n = 4:N+1
  b = sparse(M,1);
  for j = 1:size(n4e,1)
    b(n4e(j,:)) = b(n4e(j,:)) ...
  + 0.5*det([1,1,1; c4n(n4e(j,:),:)']) * dt*(f(sum(c4n(n4e(j,:),:))/3,n*dt)/6 + f(sum(c4n(n4e(j,:),:))/3,(n-1)*dt)/6);
  end
  C = C + (dt^2)*A*U(:,n-2);
  b = b + (B - (dt/2)*A + (3/4)*dt^2*A) * U(:,n-1) + C;
  u = sparse(M,1);
  u(unique(dirichlet)) = u_d(c4n(unique(dirichlet),:),n*dt);
  b = b - (0.5*dt*A + B - 0.25*dt^2*A) * u;
  u(FreeNodes) = (0.5*dt*A(FreeNodes,FreeNodes)+ ...
      B(FreeNodes,FreeNodes) - 0.25*dt^2*A(FreeNodes,FreeNodes))\b(FreeNodes); 
  U(:,n) = u;
end
show(n4e,[],c4n,full(U(:,N+1)));
