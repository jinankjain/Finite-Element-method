function eta4edge = edgeresidualdiff(c4n,n4e,n4sDb,U,j,dt)
nrElems = size(n4e,1); %nrElems := No. of elements
gradU1 = zeros(nrElems,2); 
gradU = zeros(nrElems,2);
h = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
%% Compute normals for all triangles simultaneously
normal4e = computeNormal4e(c4n,n4e);
    %% Compute the gradient
    for elem = 1 : nrElems
        grads = [1,1,1;c4n(n4e(elem,:),:)']\[0,0;eye(2)];
        gradU1(elem,:) = U(n4e(elem,:),j)' * grads;
        gradU(elem,:) = U(n4e(elem,:),j-1)' * grads;
    end
    %% Compute s4e,n4s,s4n
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];  %% Collecting all the sides
    [b,ind,back] = unique(sort(allSides,2),'rows','first'); %% Avoiding repeatation of sides
    [n4sInd, sortInd] = sort(ind); % by the way: n4s = allSides(n4sInd,:)
    n4s = allSides(n4sInd,:);
    S = size(n4s,1);
    N = max(max(n4e));
    s4n = sparse(n4s(:,1),n4s(:,2),1:S,N,N);
    % Up to here, s4n is not yet symmetric as each side has only been
    % considered once. The following makes sure that s4n is symmetric.
    s4n = s4n + s4n';
    sideNr(sortInd) = 1:length(ind); % sideNr(back): numbers for allSides
    s4e = reshape(sideNr(back),size(n4e));
%% Computing the edge residual
jump4s1 = zeros(size(n4s,1),1);
jump4s = zeros(size(n4s,1),1);
%% inner jumps
    for elem = 1 : nrElems   %%Update jump
		a = s4e(elem,:);
		b1 = n4s(a(1),:);
		b2 = n4s(a(2),:);
		b3 = n4s(a(3),:);
		if(c4n(b1(1))<=1 && c4n(b1(2))<=1 && c4n(b2(1))<=1 && c4n(b2(2))<=1 && c4n(b3(1))<=1 && c4n(b3(2))<=1)
        	jump4s(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
                              + normal4e(:,:,elem)*gradU(elem,:)';
			jump4s1(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
                              + normal4e(:,:,elem)*gradU1(elem,:)';
		else %if(c4n(b1(1))>=1 && c4n(b1(2))>=1 && c4n(b2(1))>=1 && c4n(b2(2))>=1 && c4n(b3(1))>=1 && c4n(b3(2))>=1)
			jump4s(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
                              + 0.5*normal4e(:,:,elem)*gradU(elem,:)';
        	jump4s1(s4e(elem,:)) = jump4s(s4e(elem,:)) ...
                              + 0.5*normal4e(:,:,elem)*gradU1(elem,:)';
        end
    end
 
    %% Dirichlet jumps = 0
    jump4s(diag(s4n(n4sDb(:,1),n4sDb(:,2)))) = 0;
    jump4s1(diag(s4n(n4sDb(:,1),n4sDb(:,2)))) = 0;
    length4s = sqrt(sum((c4n(n4s(:,2),:)-c4n(n4s(:,1),:)).^2, 2)); 
    jumpIntdiff4s = (jump4s1(:,1) - (1 + dt)*jump4s(:,1)).^2.*length4s.^4;
    %h = 2*sqrt(det([1 1 1;c4n(n4e(1,:),:)'])/2);
    eta4edge = h^3*abs(log(h))*(sum(jumpIntdiff4s(s4e),2));
end  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
