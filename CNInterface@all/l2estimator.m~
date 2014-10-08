function Estm = l2estimator(j,c4n,n4e,n4sDb,T,N,dt,U)
%[U,A,B,ndof] = FEMPIDE(c4n,n4e,unique(n4sDb),N,dt);
h = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
 etaE = sqrt(sum(edgeresiduall2(c4n,n4e,n4sDb,U,j)));
 
 if j == 1
     etaf = 10*h*h; 
 else
     etaf = sqrt(sum(residuall2(n4e,c4n,U,1,T,dt)));
 end
 Estm = etaE + etaf; 


