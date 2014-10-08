function Estm = space1(j,c4n,n4e,n4sDb,T,N,dt,U)
h = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
Estm = dt*(h^2*abs(log(h))*744.5);
