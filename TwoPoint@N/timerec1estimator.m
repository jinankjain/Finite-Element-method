function Estm = timerec1estimator(j,c4n,n4e,n4sDb,T,N,dt,U)
h = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
Estm = (dt^2/sqrt(30))*(h*744.5+55.5);
