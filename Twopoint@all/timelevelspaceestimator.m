function value =  timelevelspaceestimator(l,c4n,n4e,n4sDb,T,N,dt,U)
%[U,A,B,ndof] = FEMPIDE(c4n,n4e,unique(n4sDb),N,dt);
val2 = zeros(1,N+1);
 for j = 2:l
     val2(1,j) = spaceestimator(j,c4n,n4e,n4sDb,T,N,dt,U);
 end
v2 = sum(val2);
value =  v2;
