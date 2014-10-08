function value =  timelevelH1spaceestimate(c4n,n4e,n4sDb,T,N,dt,U,k)
val = zeros(1,N+1);
for j = 1:N+1
     val(1,j) = h1estimator(j,c4n,n4e,n4sDb,T,N,dt,U,k)^2;
end
value = sum(val);
