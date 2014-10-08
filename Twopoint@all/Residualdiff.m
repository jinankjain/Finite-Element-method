function etaf = Residualdiff(Elem,Coord,U,k,T,dt)
area4e = getArea4e(Coord,Elem);
if(k == 2)
    for j = 1:size(Elem,1)
        curnodes = Elem(j,:);
        curcoords=Coord(curnodes,:);
        P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
        mk=area4e(j); mp23=1/2*(P2+P3);  mp13=1/2*(P1+P3);   mp12=1/2*(P1+P2); % midpoints of edges
        h = 2*sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
        Q1=curnodes(:,1);
        Q2=curnodes(:,2);
        Q3=curnodes(:,3);
        U12N3 = (U(Q1,k) + U(Q2,k))/2;
        U13N3 = (U(Q1,k) + U(Q3,k))/2;
        U23N3 = (U(Q2,k) + U(Q3,k))/2;
        U12N2 = (U(Q1,k-1) + U(Q2,k-1))/2; 
        U13N2 = (U(Q1,k-1) + U(Q3,k-1))/2;  
        U23N2 = (U(Q2,k-1) + U(Q3,k-1))/2;  
        bacder21 = (U12N3 - U12N2)/dt;
        bacder22 = (U13N3 - U13N2)/dt;
        bacder23 = (U23N3 - U23N2)/dt;
        L2f1 = mk/3*((f(mp23,k*dt) - bacder23)^2+(f(mp13,k*dt) - bacder22)^2+(f(mp12,k*dt)-bacder21)^2);
        etaf(j,1) = h^4*L2f1/dt^2;
    end
else
    for j = 1:size(Elem,1)
    curnodes = Elem(j,:);
    curcoords=Coord(curnodes,:);
    P1=curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mk=area4e(j); mp23=1/2*(P2+P3);  mp13=1/2*(P1+P3);   mp12=1/2*(P1+P2); % midpoints of edges
    h = 2*sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
    Q1=curnodes(:,1);
    Q2=curnodes(:,2);
    Q3=curnodes(:,3);
    U12N3 = (U(Q1,k) + U(Q2,k))/2;
    U13N3 = (U(Q1,k) + U(Q3,k))/2;
    U23N3 = (U(Q2,k) + U(Q3,k))/2;
    U12N2 = (U(Q1,k-1) + U(Q2,k-1))/2;
    U13N2 = (U(Q1,k-1) + U(Q3,k-1))/2;
    U23N2 = (U(Q2,k-1) + U(Q3,k-1))/2;
    U12N1 = (U(Q1,k-2) + U(Q2,k-2))/2;
    U13N1 = (U(Q1,k-2) + U(Q3,k-2))/2;
    U23N1 = (U(Q2,k-2) + U(Q3,k-2))/2;
    bacder1 = (U12N3 - U12N2)/dt;
    bacder2 = (U13N3 - U13N2)/dt;
    bacder3 = (U23N3 - U23N2)/dt;
    bacder4 = (U12N2 - U12N1)/dt;
    bacder5 = (U13N2 - U13N1)/dt;
    bacder6 = (U23N2 - U23N1)/dt;
    L2f = (mk/3)*((f(mp23,k*dt) - f(mp23,(k-1)*dt) - bacder3 + bacder6)^2 ...
    + (f(mp13,k*dt)- f(mp13,(k-1)*dt) - bacder2 + bacder5)^2 ...
    + (f(mp12,k*dt)- f(mp12,(k-1)*dt)- bacder1 + bacder4)^2);
    etaf(j,1) =  h^4*abs(log(h))*L2f/dt^2;
    end
end