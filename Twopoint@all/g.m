function value =  g(c4n,n4e,U,T,dt,N)
value = zeros(1,N+1);
for j = 2:N+1
     value(j) = L2Residualdiff(j,n4e,c4n,U,T,dt);
end
