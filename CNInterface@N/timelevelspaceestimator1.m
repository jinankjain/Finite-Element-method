function value =  timelevelspaceestimator1(c4n,n4e,n4sDb,T,N,dt,U)
val = zeros(1,N+1);
for j = 2:N+1
    val(1,j)  = space1(j,c4n,n4e,n4sDb,T,N,dt,U);
end
value = sum(val);




