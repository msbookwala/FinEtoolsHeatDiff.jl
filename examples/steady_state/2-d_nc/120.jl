using FinEtools
using FinEtools.AlgoBaseModule: matrix_blocked, vector_blocked
using FinEtoolsHeatDiff
using FinEtoolsHeatDiff.AlgoHeatDiffModule
using FinEtools.MeshExportModule.VTK: vtkexportmesh
using LinearAlgebra
using SparseArrays
include("utilities.jl")

function f(x, y)
    return 1.0 + x^2 + 2.0*y^2
end

function Q!(forceout, XYZ, tangents, feid, qpid)
    forceout[1] = -6.0
    return forceout
end

function qbar(x, y, nx, ny)
    return 2.0*x*nx + 4.0*y*ny
end
r = 2
lam_order = 0
material = MatHeatDiff([1.0 0.0; 0.0 1.0])
fi = ForceIntensity(Float64, 1, Q!)
# function run(r=0, lam_order=0)
skew = 0.5*tan(pi/6)
    n1 = 3 * 2^r
    n2 = 3 * 2^r
    n3 = 4 * 2^r

    nlam_v = 6 * 2^r
    nlam_h = 3 * 2^r

    # ============================================================
    yi_v = collect(linearspace(0.0, 1.0, nlam_v+1))
    xi_v = 0.5*ones(nlam_v+1)
    xi_v .+=  (xi_v).*skew.*yi_v.*(yi_v.<0.5) +
                 (xi_v).*skew.*(1 .-yi_v).*(yi_v.>=0.5)
    fensi_v, fesi_v = L2blockx2D(xi_v, yi_v)
    # ============================================================
    xi_h = collect(linearspace(0.0, 0.5, nlam_h+1))
    yi_h = 0.5*ones(nlam_h+1)
    xi_h .+=  (xi_h).*skew.*(1 .-yi_h).*(yi_h.>=0.5)
    fensi_h, fesi_h = L2blockx2D(xi_h, yi_h)

    # ============================================================
    fens1, fes1 = T3block(0.5, 0.5, n1, n1)
    bfes1 = meshboundary(fes1)
    edge_fes_1v = subset(bfes1, selectelem(fens1, bfes1, box=[0.5, 0.5, 0.0, 0.5], inflate=1e-8))
    edge_fes_1h = subset(bfes1, selectelem(fens1, bfes1, box=[0.0, 0.5, 0.5, 0.5], inflate=1e-8))

    fens1.xyz[:, 1] .+= fens1.xyz[:, 1].*skew.*fens1.xyz[:, 2]

    geom1 = NodalField(fens1.xyz)
    T_1 = NodalField(zeros(size(fens1.xyz, 1),1))
    dbc_nodes_1 = selectnode(fens1; box=[0.0, 0.0, 0.0, 0.0], inflate=1e-8)
    for i in dbc_nodes_1
        setebc!(T_1, [i], 1, f(fens1.xyz[i,1], fens1.xyz[i,2]))
    end
    applyebc!(T_1)
    numberdofs!(T_1)
    femm1 = FEMMHeatDiff(IntegDomain(fes1, TriRule(3)), material)
    K1 = conductivity(femm1, geom1, T_1)
    Kff1 = matrix_blocked(K1, nfreedofs(T_1), nfreedofs(T_1))[:ff]
    Kfd1 = matrix_blocked(K1, nfreedofs(T_1), nfreedofs(T_1))[:fd]
    F1 = distribloads(femm1, geom1, T_1, fi, 2)

    nbc_elem1b  = subset(bfes1,selectelem(fens1,bfes1,
                     box=[0.0, 1.0, 0.0, 0.0], inflate=1e-8))
    femm1b = FEMMBase(IntegDomain(nbc_elem1b,  GaussRule(1, 2)))
    nx = 0.0
    ny = -1.0
    fi1b = ForceIntensity(Float64, 1,
        function (forceout, XYZ, tangents, feid, qpid)
            forceout[1] = qbar(XYZ[1], XYZ[2], nx, ny)
            return forceout
        end
    )
    F1 += distribloads(femm1b, geom1, T_1, fi1b, 1) 

    nbc_elem1l = subset(bfes1,selectelem(fens1,bfes1,
                     box=[0.0, 0.0, 0.0, 1.0], inflate=1e-8))
    femm1l = FEMMBase(IntegDomain(nbc_elem1l,  GaussRule(1, 2)))
    nx = -1.0
    ny = 0.0
    fi1l = ForceIntensity(Float64, 1,
        function (forceout, XYZ, tangents, feid, qpid)
            forceout[1] = qbar(XYZ[1], XYZ[2], nx, ny)
            return forceout
        end
    )
    F1 +=  distribloads(femm1l, geom1, T_1, fi1l, 1) 

    Ff1 = vector_blocked(F1, nfreedofs(T_1))[:f] -
             Kfd1*gathersysvec(T_1, :d)

    fensu1v, fesu1v = build_union_mesh(fensi_v,fesi_v, fens1, 
                                        edge_fes_1v, 1; 
                                        lam_order=lam_order, 
                                        to_trim =true, dir = 2)
    D1v,_,_ = build_D_matrix(fensu1v, fesu1v, 
                            fensi_v, fesi_v, 
                            fens1, edge_fes_1v;lam_order=lam_order)
    
    flam_1v = -D1v[:, dbc_nodes_1]* gathersysvec(T_1, :d)
    D1v = D1v[:, setdiff(1:size(fens1.xyz, 1), dbc_nodes_1)]


    fensu1h, fesu1h = build_union_mesh(fensi_h,fesi_h, fens1, 
                                        edge_fes_1h, 1; 
                                        lam_order=lam_order, 
                                        to_trim =true, dir = 1) 
    D1h,_,_ = build_D_matrix(fensu1h, fesu1h, 
                            fensi_h, fesi_h,
                            fens1, edge_fes_1h;lam_order=lam_order)
    flam_1h = -D1h[:, dbc_nodes_1]* gathersysvec(T_1, :d)
    D1h = D1h[:, setdiff(1:size(fens1.xyz, 1), dbc_nodes_1)]



    # ============================================================
    fens2, fes2 = T3block(0.5, 0.5, n2, n2)
    fens2.xyz[:, 2] .+= 0.5
    bfes2 = meshboundary(fes2)
    edge_fes_2v = subset(bfes2, selectelem(fens2, bfes2, box=[0.5, 0.5, 0.5, 1.0], inflate=1e-8))
    edge_fes_2h = subset(bfes2, selectelem(fens2, bfes2, box=[0.0, 0.5, 0.5, 0.5], inflate=1e-8))
    fens2.xyz[:, 1] .+= fens2.xyz[:, 1].*skew.*(1 .-fens2.xyz[:, 2])

    geom2 = NodalField(fens2.xyz)
    T_2 = NodalField(zeros(size(fens2.xyz, 1),1))
    dbc_nodes_2 = selectnode(fens2; box=[0.0, 0.0, 1.0, 1.0], inflate=1e-8)
    for i in dbc_nodes_2
        setebc!(T_2, [i], 1, f(fens2.xyz[i,1], fens2.xyz[i,2]))
    end
    applyebc!(T_2)
    numberdofs!(T_2)
    femm2 = FEMMHeatDiff(IntegDomain(fes2, TriRule(3)), material)
    K2 = conductivity(femm2, geom2, T_2)
    Kff2 = matrix_blocked(K2, nfreedofs(T_2),nfreedofs(T_2))[:ff]
    Kfd2 = matrix_blocked(K2, nfreedofs(T_2), nfreedofs(T_2))[:fd]
    F2 = distribloads(femm2, geom2, T_2, fi, 2)

    nbc_elem2t  = subset(bfes2,selectelem(fens2,bfes2,
                    box=[0.0, 1.0, 1.0, 1.0], inflate=1e-8))
    femm2t = FEMMBase(IntegDomain(nbc_elem2t,  GaussRule(1, 2)))
    nx = 0.0
    ny = 1.0
    fi2t = ForceIntensity(Float64, 1,
        function (forceout, XYZ, tangents, feid, qpid)
            forceout[1] = qbar(XYZ[1], XYZ[2], nx, ny)
            return forceout
        end
    )
    F2 += distribloads(femm2t, geom2, T_2, fi2t, 1)


    nbc_elem2l = subset(bfes2,selectelem(fens2,bfes2,
                     box=[0.0, 0.0, 0.0, 1.0], inflate=1e-8))
    femm2l = FEMMBase(IntegDomain(nbc_elem2l,  GaussRule(1, 2)))
    nx = -1.0
    ny = 0.0
    fi2l = ForceIntensity(Float64, 1,
        function (forceout, XYZ, tangents, feid, qpid)
            forceout[1] = qbar(XYZ[1], XYZ[2], nx, ny)
            return forceout
        end
    )
    F2 +=  distribloads(femm2l, geom2, T_2, fi2l, 1) 


    Ff2 = vector_blocked(F2, nfreedofs(T_2))[:f] - Kfd2*gathersysvec(T_2, :d)

    fensu2h, fesu2h = build_union_mesh(fensi_h,fesi_h, fens2, 
                                        edge_fes_2h, 1; 
                                        lam_order=lam_order, 
                                        to_trim =true, dir = 1)
    D2h,_,_ = build_D_matrix(fensu2h, fesu2h,
                            fensi_h, fesi_h,
                            fens2, edge_fes_2h;lam_order=lam_order)
    flam_2h = D2h[:, dbc_nodes_2]* gathersysvec(T_2, :d)
    D2h = -D2h[:, setdiff(1:size(fens2.xyz, 1), dbc_nodes_2)]
    
    fensu2v_all = unique(trunc.(vcat(fensi_v.xyz, fens2.xyz[unique(connasarray(edge_fes_2v)),:]), digits=12), dims=1)
    # sort and then take only where y>=0.5
    fensu2v_all = fensu2v_all[sortperm(fensu2v_all[:, 2]), :]
    fensu2v_nodes = fensu2v_all[fensu2v_all[:, 2] .>= 0.5, :]
    fensu2v, fesu2v = L2blockx2D(fensu2v_nodes[:, 1], fensu2v_nodes[:, 2])
    D2v,_,_ = build_D_matrix(fensu2v, fesu2v,
                            fensi_v, fesi_v,
                            fens2, edge_fes_2v;lam_order=lam_order)
    flam_2v = -D2v[:, dbc_nodes_2]* gathersysvec(T_2, :d)
    D2v = D2v[:, setdiff(1:size(fens2.xyz, 1), dbc_nodes_2)]


    # ============================================================
    fens3, fes3 = T3block(0.5, 1.0, ceil(Int,n3/2), n3)
    fens3.xyz[:, 1] .+= 0.5
    bfes3 = meshboundary(fes3)

    edge_fes_3v = subset(bfes3, selectelem(fens3, bfes3, box=[0.5, 0.5, 0.0, 1.0], inflate=1e-8))
    fens3.xyz[:,1] .+= (1 .-fens3.xyz[:,1]).*skew.*fens3.xyz[:,2].*(fens3.xyz[:,2].<0.5) + 
                    (1 .-fens3.xyz[:,1]).*skew.*(1 .-fens3.xyz[:,2]).*(fens3.xyz[:,2].>=0.5)


    geom3 = NodalField(fens3.xyz)
    T_3 = NodalField(zeros(size(fens3.xyz, 1),1))
    dbc_nodes_3 = selectnode(fens3; box=[1.0, 1.0, 0.0, 1.0], inflate=1e-8)
    for i in dbc_nodes_3
        setebc!(T_3, [i], 1, f(fens3.xyz[i,1], fens3.xyz[i,2]))
    end
    applyebc!(T_3)
    numberdofs!(T_3)
    femm3 = FEMMHeatDiff(IntegDomain(fes3, TriRule(3)), material)
    K3 = conductivity(femm3, geom3, T_3)
    Kff3 = matrix_blocked(K3, nfreedofs(T_3), nfreedofs(T_3))[:ff]
    Kfd3 = matrix_blocked(K3, nfreedofs(T_3), nfreedofs(T_3))[:fd]
    F3 = distribloads(femm3, geom3, T_3, fi,2)

    nbc_elem3t  = subset(bfes3,selectelem(fens3,bfes3,
                     box=[0.0, 1.0, 1.0, 1.0], inflate=1e-8))
    femm3t = FEMMBase(IntegDomain(nbc_elem3t,  GaussRule(1, 2)))
    nx = 0.0
     ny = 1.0
    fi3t = ForceIntensity(Float64, 1,
        function (forceout, XYZ, tangents, feid, qpid)
            forceout[1] = qbar(XYZ[1], XYZ[2], nx, ny)
            return forceout
        end
    )
    F3 += distribloads(femm3t, geom3, T_3, fi3t, 1)

    nbc_elem3b  = subset(bfes3,selectelem(fens3,bfes3,
                     box=[0.0, 1.0, 0.0, 0.0], inflate=1e-8))
    femm3b = FEMMBase(IntegDomain(nbc_elem3b,  GaussRule(1, 2)))
    nx = 0.0
    ny = -1.0
    fi3b = ForceIntensity(Float64, 1,
        function (forceout, XYZ, tangents, feid, qpid)
            forceout[1] = qbar(XYZ[1], XYZ[2], nx, ny)
            return forceout
        end
    )
    F3 += distribloads(femm3b, geom3, T_3, fi3b, 1)
    Ff3 = vector_blocked(F3, nfreedofs(T_3))[:f] - Kfd3*gathersysvec(T_3, :d)    
    
    fensu3v, fesu3v = build_union_mesh(fensi_v,fesi_v, fens3, 
                                        edge_fes_3v, 1; 
                                        lam_order=lam_order, 
                                        to_trim =false, dir = 2)
    D3v,_,_ = build_D_matrix(fensu3v, fesu3v, 
                            fensi_v, fesi_v, 
                            fens3, edge_fes_3v;lam_order=lam_order)
    flam_3v = D3v[:, dbc_nodes_3]* gathersysvec(T_3, :d)
    D3v = -D3v[:, setdiff(1:size(fens3.xyz, 1), dbc_nodes_3)]

    K = [Kff1 spzeros(length(Ff1), length(Ff2)+length(Ff3));
         spzeros(length(Ff2),    length(Ff1)) Kff2 spzeros(length(Ff2), length(Ff3));
         spzeros(length(Ff3), length(Ff1)+length(Ff2)) Kff3]
    F = vcat(Ff1, Ff2, Ff3)
    Flam = vcat(flam_1v+flam_2v+flam_3v, flam_1h+flam_2h)
    D = [D1v D2v D3v;
         D1h D2h spzeros(size(D1h, 1), size(D3v, 2))]
    A = [K D'; D spzeros(size(D, 1), size(D, 1))]
    rhs = vcat(F, Flam)
    sol = A\rhs
    scattersysvec!(T_1, sol[1:length(Ff1)])
    scattersysvec!(T_2, sol[length(Ff1)+1 : length(Ff1)+length(Ff2)])
    scattersysvec!(T_3, sol[length(Ff1)+length(Ff2)+1 : length(Ff1)+length(Ff2)+length(Ff3)])

    if !isdir("120")
        mkdir("120")
    end

    vtkexportmesh(
        "120/1.vtk",
        fens1, fes1;
        scalars=[
            ("Temperature", T_1.values),
            # ("Err", err1.values)
        ]
    )

    vtkexportmesh(
        "120/2.vtk",
        fens2, fes2;
        scalars=[
            ("Temperature", T_2.values),
            # ("Err", err2.values)
        ]
    )

    vtkexportmesh(
        "120/3.vtk",
        fens3, fes3;
        scalars=[
            ("Temperature", T_3.values),
            # ("Err", err3.values)
        ]   
    )
    vtkexportmesh(
        "120/v.vtk",
        fensi_v, fesi_v;
        scalars=[
            # ("lambda_kink", uk.values)
        ]
    )

    vtkexportmesh(
        "120/h.vtk",
        fensi_h, fesi_h;
        scalars=[
            # ("lambda_horizontal", uh.values)
        ]
    )

    vtkexportmesh(
        "120/3v.vtk",
        fensu3v, fesu3v;
        scalars=[
            # ("lambda_vertical", u3v.values)
        ]
    )

    vtkexportmesh(
        "120/1v.vtk",
        fensu1v, fesu1v;
        scalars=[
            # ("lambda_vertical", u1v.values)
        ]
    )
    vtkexportmesh(
        "120/1h.vtk",
        fensu1h, fesu1h;
        scalars=[
            # ("lambda_horizontal", u1h.values)
        ]
    )
    vtkexportmesh(
        "120/2h.vtk",
        fensu2h, fesu2h;
        scalars=[
            # ("lambda_horizontal", u2h.values)
        ]
    )
    vtkexportmesh(
        "120/2v.vtk",
        fensu2v, fesu2v;
        scalars=[
            # ("lambda_vertical", u2v.values)
        ]
    )
# end
# run()