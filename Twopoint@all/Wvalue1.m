function value =  Wvalue1(c4n,n4e,n4sDb,T,N,dt)
value = zeros(1,N+1);
for j = 1:N+1
     value(1,j) = WL2(j,c4n,n4e,n4sDb,T,N,dt);
end