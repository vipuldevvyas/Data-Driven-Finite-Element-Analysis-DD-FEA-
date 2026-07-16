function K = assembleK(mesh, C)

ndof    = mesh.ndof;
numElem = mesh.numElements;

K = sparse(ndof, ndof);

for e = 1:numElem

    conn   = mesh.elements(e,:);      % 1x10 for Tet10
    coords = mesh.nodes(conn,:);      % 10x3

    % Get all 4 Gauss point B matrices
    [B_all, V_all, ~] = tet10_Bmatrix(coords);

    % Element stiffness = sum over 4 Gauss points
    Ke = zeros(30, 30);

    for gp = 1:4
        B  = B_all{gp};
        Vg = V_all(gp);
        Ke = Ke + Vg * (B' * C * B);
    end

    edofs = mesh.elementDOFs{e};      % 1x30 for Tet10

    K(edofs, edofs) = K(edofs, edofs) + Ke;

end

K = (K + K') / 2;

end