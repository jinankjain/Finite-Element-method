function Estm = spaceestimator(j,c4n,n4e,n4sDb,T,N,dt,U)
%[U,A,B,ndof] = FEMPIDE(c4n,n4e,unique(n4sDb),N,dt);
etaE(j) = sqrt(sum(edgeresidualdiff(c4n,n4e,n4sDb,U,j,dt)));
etaf(j) = sqrt(sum(Residualdiff(n4e,c4n,U,j,T,dt)));
%%+etaf(j)
%Estm = etaE(j);
Estm = etaf(j) + etaE(j);
%Estm = etaf(j);    
