using LinearAlgebra
using SparseArrays

cross2(a, b) = a[1]*b[2] - a[2]*b[1]


function intersect_seg_line(B1, B2, G1, G2; eps=1e-14)
    d1 = B2-B1
    d2 = G2-G1
    denom = cross(d1, d2)
    if norm(denom) < eps
        return B2
    end
    t = norm(cross(G1-B1, d2)) / norm(denom)
    return B1 + t * d1
end


XA = [0.0 0.0 0.0;
      1.0 0.0 0.0;
      1.0 1.0 0.0;
      0.0 1.0 0.0]
TA = [1 2 3;
      1 3 4]

TB = [1 2 4;
      2 3 4]

# all the vertices and subvertices
verts = []

for i in 1:size(TA,1)
    V_A = XA[TA[i,:], :]
    # iterate of edges of A
    for eA in 1:3
        B1 = V_A[eA, :]
        B2 = V_A[mod1(eA+1,3), :]

        # iterate over triangles of B
        for j in 1:size(TB,1)
            V_B = XA[TB[j,:], :]
            # iterate over edges of B
            for eB in 1:3
                G1 = V_B[eB, :]
                G2 = V_B[mod1(eB+1,3), :]

                P = intersect_seg_line(B1, B2, G1, G2)
                push!(verts, P)
            end
        end
    end
end