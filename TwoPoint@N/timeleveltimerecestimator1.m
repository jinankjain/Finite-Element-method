function value =  timeleveltimerecestimator1(c4n,n4e,n4sDb,T,N,dt,U)
%[U,A,B,ndof] = FEMPIDE(c4n,n4e,unique(n4sDb),N,dt);
val2 = zeros(1,N+1);
 for j = 2:N+1
     val2(1,j) = timerec1estimator(j,c4n,n4e,n4sDb,T,N,dt,U);
 end
value = sum(val2);
%value =  v2;
