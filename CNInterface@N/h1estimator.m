function Estm = h1estimator(j,c4n,n4e,n4sDb,T,N,dt,U,k)
	h = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
	etaE = (sqrt(h))*sqrt(sum(edgeresidualh1(c4n,n4e,n4sDb,U,j)));
	if(j==1)
		etaf = (5*pi*pi*(1+2.5*pi*pi))*h;
	else
		etaf = sqrt(sum(residualh1(n4e,c4n,U,k,j,dt)));
    end
    %Estm = etaE(j);
    Estm = etaE+etaf;
end



    

