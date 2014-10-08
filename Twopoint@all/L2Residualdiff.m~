function etaf = L2Residualdiff(k,Elem,Coord,U,T,dt)
area4e = getArea4e(Coord,Elem);
if(k == 2)
    for j = 1:size(Elem,1)
        curnodes = Elem(j,:);
        curcoords=Coord(curnodes,:);
        P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
        mk=area4e(j); mp23=1/2*(P2+P3);  mp13=1/2*(P1+P3);   mp12=1/2*(P1+P2); % midpoints of edges
        h = 2*sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
        L2f1 = mk/3*((f1(mp23,k*dt))^2+(f1(mp13,k*dt))^2+(f1(mp12,k*dt))^2);
        etaf(j,1) = (L2f1/dt^2);
    end
else
    for j = 1:size(Elem,1)
    curnodes = Elem(j,:);
    curcoords=Coord(curnodes,:);
    P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mk=area4e(j); mp23=1/2*(P2+P3);  mp13=1/2*(P1+P3);   mp12=1/2*(P1+P2); % midpoints of edges
    h = 2*sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
    L2f = (mk/3)*((f1(mp23,k*dt) - f1(mp23,(k-1)*dt))^2 ...
    + (f1(mp13,k*dt)- f1(mp13,(k-1)*dt))^2 ...
    + (f1(mp12,k*dt)- f1(mp12,(k-1)*dt))^2);
    etaf(j,1) =  (L2f/dt^2);
    end
end
etaf = sqrt(sum(etaf));
