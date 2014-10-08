function value =  timelevell2spaceestimate(c4n,n4e,n4sDb,T,N,dt,U)
for j = 1:N+1
  value = l2estimator(j,c4n,n4e,n4sDb,T,N,dt,U);
end
