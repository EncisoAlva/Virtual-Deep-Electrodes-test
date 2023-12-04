%% 
% creation of an example
% this select the upper part of the pial

xport = xport_pre;

N = size(xport.Vertices,1);
idx = xport.Vertices(:,3)>0.01;

xport.Vertices = xport.Vertices(idx,:);
xport.VertConn = xport.VertConn(idx,idx);
xport.VertNormals = xport.VertNormals(idx,:);
xport.Curvature = xport.Curvature(idx,:);
xport.SulciMap = xport.SulciMap(idx,:);

idx_rev = zeros(N,1 );
counter = 1;
for i = 1:N
    if idx(i)
        idx_rev(i) = counter;
        counter = counter+1;
    end
end
for i = 1:size(xport.Faces,1)
    xport.Faces(i,1) = idx_rev( xport.Faces(i,1) );
    xport.Faces(i,2) = idx_rev( xport.Faces(i,2) );
    xport.Faces(i,3) = idx_rev( xport.Faces(i,3) );
end

idx_fac = (xport.Faces(:,1)==0);
for i = 1:size(idx_fac)
    if (xport.Faces(i,1)~=0) && (xport.Faces(i,2)~=0) && (xport.Faces(i,3)~=0)
        idx_fac(i) = true;
    else
        idx_fac(i) = false;
    end
end
xport.Faces = xport.Faces(idx_fac,:);

cortex_vert = xport.Vertices;
cortex_conn = xport.Faces;
if(dev)
    % DELETE
    figure()
    trisurf(cortex_conn, cortex_vert(:,1),cortex_vert(:,2),cortex_vert(:,3))
    scatter3(cortex_vert(:,1),cortex_vert(:,2),cortex_vert(:,3),'filled')
    title('Only cortex surface (mesh vertices)')
    xlabel('x')
    ylabel('y')
    zlabel('z')
end


example1 = xport;

save('example1.mat','example1')