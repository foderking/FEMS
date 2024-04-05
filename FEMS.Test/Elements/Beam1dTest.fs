module Test.Elements.Beam1d

open Xunit
open FEMS.Core
open MathNet.Numerics.LinearAlgebra
open Test.Base

[<Fact>]
let ``beam with supports at both ends and downward uniformly distributed load``() =
    (*
                  q
        ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓

        1    3    5    7    9 <- angular displacement  (θ)
        |----|----|----|----|
        0    2    4    6    8 <- vertical displacement (w)
        ▲                   ▲ <- fixed supports (w(0) = w(L) = 0) 
    *)
    let n_elems = 4;
    let EI = 20000.0
    let q = -5.0
    let L = 2.0
    let w_0, w_l = 0, 8
    let θ_0, θ_l = 1, 9
    let w_mid = 4

    let elem = Elements.Beam1d(L, n_elems, EI, beamType=Elements.UniformlyDistributed(q))
    let n_nodes = elem.totalNodes

    let global_k = DenseMatrix.create n_nodes n_nodes 0.0
    let local_k  = (elem :> Element).localStiffNess
    let global_f = DenseVector.create n_nodes 0.0
    let local_f  = (elem :> Element).localForce

    for i in 0..n_elems-1 do
        Operations.assemble(local_k, global_k, local_f, global_f,[2*i..2*i+3])

    Operations.enforceBC(global_k,global_f, [[w_0;0];[w_l;0]],[])

    let X = Solve.staticEquilibrum(global_k,global_f)
    let w_mid_exact = (5.0*q*L*L*L*L)/(384.0*EI)
    let θ_0_exact = (q*L*L*L)/(24.0*EI)

    Assert.True(looseEqual(X[θ_l], -X[θ_0]))
    Assert.True(looseEqual(w_mid_exact, X[w_mid]))
    Assert.True(looseEqual(θ_0_exact, X[θ_0]))

[<Fact>]
let ``cantilever fixed at left end and downward uniformly distributed load``() =
    (*
                  q
        ↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓↓
      |-
      |-1    3    5    7    9 <- angular displacement  (θ)
      |-|----|----|----|----|
      |-0    2    4    6    8 <- vertical displacement (w)
      |-
 fixed end     
    *)
    let n_elems = 4;
    let EI = 20000.0
    let q = -5.0
    let L = 2.0
    let w_0, w_l = 0, 8
    let θ_0, θ_l = 1, 9
    let w_mid = 4

    let elem = Elements.Beam1d(L, n_elems, EI, beamType=Elements.UniformlyDistributed(q))
    let n_nodes = elem.totalNodes

    let global_k = DenseMatrix.create n_nodes n_nodes 0.0
    let local_k  = (elem :> Element).localStiffNess
    let global_f = DenseVector.create n_nodes 0.0
    let local_f  = (elem :> Element).localForce

    for i in 0..n_elems-1 do
        Operations.assemble(local_k, global_k, local_f, global_f,[2*i..2*i+3])

    Operations.enforceBC(global_k,global_f, [[w_0;0];[θ_0;0]],[])

    let X = Solve.staticEquilibrum(global_k,global_f)
    let w_l_exact = (q*L*L*L*L)/(8.0*EI)
    let θ_l_exact = (q*L*L*L)/(6.0*EI)

    Assert.True(looseEqual(w_l_exact, X[w_l]))
    Assert.True(looseEqual(θ_l_exact, X[θ_l]))
