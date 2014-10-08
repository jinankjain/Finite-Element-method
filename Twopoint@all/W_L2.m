function etaf = W_L2(c4n,n4e,l,dt)
WL2norm = 0; 
for j=1:size(n4e,1)
    curnodes = n4e(j,:);
    curcoords = c4n(curnodes,:);
    for k = 1:3
        Cuh(k) = .5*(f(curcoords(k,:),l*dt)-f(curcoords(k,:),(l-1)*dt))/dt;
    end
    P1 = curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
    mk=1/2*det([1 P1;1 P2;1 P3]);
    Wmp12=1/2*(Cuh(1)+Cuh(2)); Wmp13=1/2*(Cuh(1)+Cuh(3));
    Wmp23=1/2*(Cuh(3)+Cuh(2));
    WL2norm = WL2norm+mk/3*(Wmp12^2+Wmp13^2+Wmp23^2);
end
etaf = sqrt(WL2norm);

% L2e=0;   H1e=0; 
% for j=1:size(n4e,1)
%     curnodes = n4e(j,:);
%     curcoords = c4n(curnodes,:);
%     Cuh = U(curnodes);
%     P1 = curcoords(1,:); P2=curcoords(2,:); P3=curcoords(3,:);
% 
%     uhmp12=1/2*(Cuh(1)+Cuh(2)); uhmp13=1/2*(Cuh(1)+Cuh(3));
%     uhmp23=1/2*(Cuh(3)+Cuh(2)); %%uhcg=1/3*(Cuh(1)+Cuh(2)+Cuh(3));
%     mk=1/2*det([1 P1;1 P2;1 P3]);
%     L2e=L2e+mk/3*((ue(mp12,l*dt)- uhmp12)^2+(ue(mp13,l*dt)-uhmp13)^2+(ue(mp23,l*dt)-uhmp23)^2);
% end
%    L2e=sqrt(L2e); H1e=sqrt(H1e);
% end