function etaf = W_H1(c4n,n4e,l,dt)
WH1norm = 0; 
for j=1:size(n4e,1)
    curnodes = n4e(j,:);
    curcoords = c4n(curnodes,:);
    for k = 1:3
        Cuh2(k) = .5*(f(curcoords(k,:),l*dt)-f(curcoords(k,:),(l-1)*dt))/dt;
    end
    %Cuh2 = .5*(f(curcoords,l*dt)-f(curcoords,(l-1)*dt));
    P1 = curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mk=1/2*det([1 P1;1 P2;1 P3]);
    L1=[ones(1,3);curcoords']'\[1;0;0];  %% L1 = inv(A)*[1;0;0]; where A = [ones(1,3); curcoords']'
    L2=[ones(1,3);curcoords']'\[0;1;0];
    L3=[ones(1,3);curcoords']'\[0;0;1];
    uxh2=Cuh2(1)*L1(2)+Cuh2(2)*L2(2)+Cuh2(3)*L3(2);
    uyh2=Cuh2(1)*L1(3)+Cuh2(2)*L2(3)+Cuh2(3)*L3(3);
    WH1norm = WH1norm+mk*uxh2^2+mk*uyh2^2;
end
etaf = sqrt(WH1norm);