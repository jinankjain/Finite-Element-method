function value =  Wvalue(c4n,n4e,n4sDb,T,N,dt)
value = zeros(1,N+1);
for j = 1:N+1
     value(1,j) = WH1(j,c4n,n4e,n4sDb,T,N,dt);
end
