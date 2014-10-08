function Estm = timerec2estimator(j,c4n,n4e,n4sDb,T,N,dt,U)
h = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
Estm = (dt^2/2)*(9.1+h^2*sqrt(abs(log(h)))*744.5);
