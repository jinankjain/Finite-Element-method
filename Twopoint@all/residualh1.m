function etaf = residualh1(Elem,Coord,U,k,T,dt)
area4e = getArea4e(Coord,Elem);
for j=1:size(Elem,1),
        area = area4e(j);
end
if(k == 1)
    for j = 1:size(Elem,1)
        etaf(j,1) = 0;
    end
else
    for j=1:size(Elem,1),
    curnodes = Elem(j,:);
    curcoords=Coord(curnodes,:);
    P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mk=area4e(j); mp23=1/2*(P2+P3);  mp13=1/2*(P1+P3);   mp12=1/2*(P1+P2); % midpoints of edges
    h = 2*sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
    Q1=curnodes(:,1);
    Q2=curnodes(:,2);
    Q3=curnodes(:,3);
    U12N1 = (U(Q1,k) + U(Q2,k))/2;
    U13N1 = (U(Q1,k) + U(Q3,k))/2;
    U23N1 = (U(Q2,k) + U(Q3,k))/2;
    U12N = (U(Q1,k-1) + U(Q2,k-1))/2;
    U13N = (U(Q1,k-1) + U(Q3,k-1))/2;
    U23N = (U(Q2,k-1) + U(Q3,k-1))/2;
    bacder1 = (U12N1 - U12N)/dt;
    bacder2 = (U13N1 - U13N)/dt;
    bacder3 = (U23N1 - U23N)/dt;
    L2f = mk/3*((f(mp12,k*dt) - bacder1)^2+(f(mp13,k*dt) - bacder2)^2+(f(mp23,k*dt)-bacder3)^2);
    etaf(j,1) = h^2.*L2f;
end
end
end
