function etaf = Residualdiff1(Elem,Coord,U,k,T,dt)
area4e = getArea4e(Coord,Elem);
% for j=1:size(Elem,1),
%         area = area4e(j);
% end
if(k == 2)
    for j = 1:size(Elem,1)
        curnodes  =  Elem(j,:);
        curcoords = Coord(curnodes,:);
        P1=curcoords(1,:);
		P2=curcoords(2,:);
		P3=curcoords(3,:);
        mk=area4e(j); 
		mp23=1/2*(P2+P3);
		mp13=1/2*(P1+P3);
		mp12=1/2*(P1+P2); % midpoints of edges
        %h = 2*sqrt(det([1 1 1;Coord(Elem(1,:),:)'])/2);
        Q1=curnodes(:,1);
        Q2=curnodes(:,2);
        Q3=curnodes(:,3);
        U12N3 = (U(Q1,k) + U(Q2,k))/2;
        U13N3 = (U(Q1,k) + U(Q3,k))/2;
        U23N3 = (U(Q2,k) + U(Q3,k))/2;
        U12N2 = (U(Q1,k-1) + U(Q2,k-1))/2;   % 0
        U13N2 = (U(Q1,k-1) + U(Q3,k-1))/2;  % 0
        U23N2 = (U(Q2,k-1) + U(Q3,k-1))/2;   % 0
		%full(U23N2)
        bacder1 = (full(U12N3) - full(U12N2))/dt;
        bacder2 = (full(U13N3) - full(U13N2))/dt;
        bacder3 = (full(U23N3) - full(U23N2))/dt;
        L2f = mk/3*((f(mp23,k*dt) - bacder3)^2+(f(mp13,k*dt) - bacder2)^2+(f(mp12,k*dt)-bacder1)^2);
        etaf(j,1) = L2f;
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
    L2f = mk/3*((f(mp23,k*dt) - f(mp23,(k-1)*dt) - bacder3 + bacder6)^2 ...
    + (f(mp13,k*dt)- f(mp13,(k-1)*dt) - bacder2 + bacder5)^2 ...
    + (f(mp12,k*dt)- f(mp12,(k-1)*dt)- bacder1 + bacder4)^2);
    etaf(j,1) = L2f;
end
end
end
