nodes = surf_new.Vertices;
edges = surf_new.Faces;

[nodes2, edges2] = meshcheckrepair(nodes, edges);

[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'dupnodes');
[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'dupelemen');
[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'dup');
[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'isolated');
[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'open');
[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'deep');
[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'meshfix');
[nodes2, edges2] = meshcheckrepair(nodes2, edges2, 'intersect');

surf_new.Vertices = nodes2;
surf_new.Faces    = edges2;
surf_new.Comment  = 'cortex_concat04';