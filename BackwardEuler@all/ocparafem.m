function [maxL2e, L2H1e] = ocparafem(N,c4n,n4e,n4sDb,T,dt,U)
%[U,A,B,ndof] = FEMPARABOLIC(c4n,n4e,unique(n4sDb),N,dt);
C = zeros(N+1,1);
D = zeros(N+1,1);
D1 = zeros(N,1);
for l = 1:N+1
   [C(l,1), D(l,1)] = Errj(c4n, n4e, U(:,l), l, dt);
end
D;
D1(1:N,1) = D(1:N,1);
maxL2e = max(C);
L2H1e = sqrt(sum(dt.*D1.^2));


