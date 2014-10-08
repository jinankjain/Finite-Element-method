function n4sMarked = markUniform(n4e)
%% markUniform - Mark given parts uniformly.
%   n4sMarked = markUniform(n4p) marks every given part for refinement,
%       regardless of estimated errors. n4p contains the nodes for the
%       parts: 2 entries per row for sides, 3 entries per row for elements.
%       The output is a list of marked sides given by their end nodes.

    nrParts = size(n4e,1);

    %% Uniform criterion: mark everything
    I = 1:nrParts;

    %% Mark sides
        allSidesMarked = [n4e(I,[1 2]);n4e(I,[2 3]);n4e(I,[3 1])];
        % Eliminate duplicates.
        [b, ind]   = unique(sort(allSidesMarked,2), 'rows');
        n4sMarked = allSidesMarked(sort(ind),:);
end

