function show(elements3,elements4,coordinates,u)

% Octave inbuilt function used to generate the mesh
trisurf(elements3,coordinates(:,1),coordinates(:,2),u','facecolor','interp')
hold on
trisurf(elements4,coordinates(:,1),coordinates(:,2),u','facecolor','interp')
hold off
view(10,40);
title('Finite Element Solution of the Problem')
