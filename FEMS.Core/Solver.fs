module FEMS.Core.Solve

open MathNet.Numerics.LinearAlgebra

type MatrixType = 
    | Regular = 0
    | Banded  = 1

let staticEquilibrum(k: MatrixD, f: VectorD): VectorD =
    k.Solve(f)

let eigenvalueProblem(k: MatrixD, m: MatrixD, w: double): unit =
    ()

