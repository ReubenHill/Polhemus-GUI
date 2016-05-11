% Head mesh transformation sequence for polhemus gui
% landmarks = Inion, Nasion, Ar, Al, Cz measured with Polhemus
% pts = landmarks of atlas (from refpts_landmarks.mat)
% mesh = adult atlas scalp mesh

[A,B] = affinemap(pts,landmarks);

mesh_trans = mesh;
mesh_trans.node = affine_trans_RJC(mesh.node,A,B);

%Then plot the transformed mesh as visual reference for further points...
trisurf(mesh_trans.face, mesh_trans.node(:,1), mesh_trans.node(:,2), mesh_trans.node(:,3),'FaceColor',[0.7 0.7 0.7]);
axis equal
hold on;
plot3(landmarks(:,1), landmarks(:,2), landmarks(:,3), 'm.', 'MarkerSize', 20);
hold off;