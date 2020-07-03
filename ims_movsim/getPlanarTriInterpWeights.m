function [weights, tri_faces, ptInTri] = getPlanarTriInterpWeights(tri_pts, intrp_pts)

% number of interpolation points
nInt = size(intrp_pts,1);
% number of points forming the triangles
nTriPts = size(tri_pts,1);
% form triangles by Delaunay triangulation
tri_faces = delaunay(tri_pts(:,1), tri_pts(:,2));
tri_faces = sortTriangles(tri_faces);
nTri = size(tri_faces,1);

% compute barycentric coords for all triangles
invT = zeros(2,2,nTri);
for nt=1:nTri
    r1 = tri_pts(tri_faces(nt,1),1:2);
    r2 = tri_pts(tri_faces(nt,2),1:2);
    r3 = tri_pts(tri_faces(nt,3),1:2);
    invT(:,:,nt) = inv([    r1(1)-r3(1)    r2(1)-r3(1);
        r1(2)-r3(2)    r2(2)-r3(2)]);
end
L = zeros(nInt,3,nTri);
for ni=1:nInt
    for nt=1:nTri
        r3 = tri_pts(tri_faces(nt,3),1:2);
        L(ni,1:2,nt) = invT(:,:,nt)*[intrp_pts(ni,1:2)-r3]';
        L(ni,3,nt) = 1-sum(L(ni,1:2,nt),2);
    end
end

% find if point is in triangle (and in which triangle)
% in case the point is outside all triangles, find which triangle it is
% closer to
ptInTri = zeros(nInt,2);
for ni=1:nInt
    for nt=1:nTri
        if abs( sum(abs(L(ni,:,nt)),2)-1 )<10e-5
            ptInTri(ni,1) = 1;
            ptInTri(ni,2) = nt;
            break
        else
            [~,ptInTri(ni,2)] = min(squeeze(sum(abs(L(ni,:,:)),2)));
        end
    end
end

% compute inteporlation weights
% for points outside the triangles, project the point to the closest
% triangle edge, and compute interpolation weights from the distance of the
% projected point to the two vertices of that edge
G = zeros(nInt,nTriPts);
for ni=1:nInt
    if ptInTri(ni,1)
        g = L(ni,:,ptInTri(ni,2));
        g(g<0)=0;
        G(ni,tri_faces(ptInTri(ni,2),:)) = g;
    else
        active_tripts = find(L(ni,:,ptInTri(ni,2))>0);
        active_pts = tri_faces(ptInTri(ni,2),active_tripts);
        vec12 = tri_pts(active_pts(2),:) - tri_pts(active_pts(1),:);
        mag12 = sqrt(sum(vec12.^2));
        vec1x = intrp_pts(ni,:) - tri_pts(active_pts(1),:);
        proj1x_12 = (vec1x*vec12')/mag12;
        G(ni,active_pts(1)) = 1-proj1x_12/mag12;
        G(ni,active_pts(2)) = proj1x_12/mag12;
    end
end
weights = G;
end

function tri_faceidx_sorted = sortTriangles(tri_faceidx)
nTri = size(tri_faceidx,1);
tri_faceidx_sorted = zeros(nTri,3);
for nt=1:nTri
    [~,sorted] = sort(tri_faceidx(nt,:));
    tri_faceidx_sorted(nt,:) = tri_faceidx(nt,sorted);
end
[~,sorted2] = sort(tri_faceidx_sorted(:,1));
tri_faceidx_sorted = tri_faceidx_sorted(sorted2,:);
end
