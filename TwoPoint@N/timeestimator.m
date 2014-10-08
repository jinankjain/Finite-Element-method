function Estm = timeestimator(j,c4n,n4e,n4sDb,T,N,dt,U)
%[U,A,B,ndof] = FEMPARABOLIC(c4n,n4e,unique(n4sDb),N,dt);
etaf2 = .5*sqrt(sum(Residualdiff1(n4e,c4n,U,2,T,dt)));
 if j == 2
     Estm = etaf2;
 else
     Estm = .5*dt*sqrt(sum(Residualdiff1(n4e,c4n,U,j,T,dt)));
 end





   
