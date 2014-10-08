function[U, A, B, ndof] = FEMPARABOLIC(c4n,n4e,dirichlet,N,dt)
M = size(c4n,1); d = size(c4n,2); A = sparse(M,M); B = sparse(M,M); x = zeros(M,1); 
interFaceNodes = find(c4n(:,1)==1);
FreeNodes = setdiff(1:M,unique(dirichlet));
FreeNodes = setdiff(FreeNodes,interFaceNodes);
SecondHalf = find(c4n(:,1)>1);
FreeNodes = setdiff(FreeNodes,SecondHalf);
U = sparse(M,N+1);
ndof = length(FreeNodes);
for j = 1:size(n4e,1)
  A(n4e(j,:),n4e(j,:)) = A(n4e(j,:),n4e(j,:))...
      + stima3(c4n(n4e(j,:),:));
end
for j = 1:size(n4e,1)
  B(n4e(j,:),n4e(j,:)) = B(n4e(j,:),n4e(j,:)) ...
      + det([1,1,1;c4n(n4e(j,:),:)'])*[2,1,1;1,2,1;1,1,2]/24;
end
x = find(c4n(:,1)<1); 
U(x,1) = sin(pi*c4n(x,1)).*sin(pi*c4n(x,2));%%%zeros(M,1);%% sin(pi*c4n(:,1)).*sin(pi*c4n(:,2)); % zeros(M,1); 

U(find(c4n(:,1)==1),1) = 0;
for n = 2:N+1
  b = sparse(M,1);
  for j = 1:size(n4e,1)
    b(n4e(j,:)) = b(n4e(j,:)) + det([1,1,1; c4n(n4e(j,:),:)']) * dt*f(sum(c4n(n4e(j,:),:))/3,n*dt)/6;
  end
  b = b + B * U(:,n-1);
  u = sparse(size(c4n,1),1);
  u(unique(dirichlet)) = 0;
  u(SecondHalf) = 0;
  u(interFaceNodes) = 0;
  b = b - (dt * A + B) * u;
  u(FreeNodes) = (dt*A(FreeNodes,FreeNodes)+ ...
      B(FreeNodes,FreeNodes))\b(FreeNodes); 
  U(:,n) = u;
end
%show(n4e,[],c4n,full(U(:,N+1)));
end






