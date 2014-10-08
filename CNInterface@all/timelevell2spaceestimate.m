function value =  timelevell2spaceestimate(l,c4n,n4e,n4sDb,T,N,dt,U)
for j = 1:l
  value = l2estimator(j,c4n,n4e,n4sDb,T,N,dt,U);
end
