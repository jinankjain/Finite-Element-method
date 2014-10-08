function normal4e = computeNormal4e(c4n,n4e)
%% computeNormal4e - normals for elements.
%   computeNormal4e(c4n, n4e) computes the three outer unit normal vectors of
%                         each element of the decomposition. c4n and n4e
%                         are as specified in the documentation.
   if isempty(n4e)
        normal4e = zeros(0,2);
        return;
   end

    %% Compute normal4e.
    allSides = [n4e(:,[1 2]); n4e(:,[2 3]); n4e(:,[3 1])];
    c4start  = c4n(allSides(:,1),:);
    c4end    = c4n(allSides(:,2),:);
    lengths  = sqrt(sum((c4end-c4start).^2,2));
    tangents = (c4end - c4start)./[lengths lengths];
    normals  = [tangents(:,2), -tangents(:,1)];
    normal4e(1,:,:) = normals(1:size(n4e,1),:)';
    normal4e(2,:,:) = normals(size(n4e,1)+1:2*size(n4e,1),:)';
    normal4e(3,:,:) = normals(2*size(n4e,1)+1:3*size(n4e,1),:)';
end
