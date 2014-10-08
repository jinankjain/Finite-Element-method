function value =  timeleveltimerecestimator2(c4n,n4e,n4sDb,T,N,dt,U)
%[U,A,B,ndof] = FEMPIDE(c4n,n4e,unique(n4sDb),N,dt);
value = zeros(1,N+1);
 for j = 2:N+1
     value(1,j) = timerec2estimator(j,c4n,n4e,n4sDb,T,N,dt,U);
 end

