module Solver

open Xunit
open MathNet.Numerics.LinearAlgebra

[<Fact>]
let ``My test`` () =
    let A: Matrix<double> = DenseMatrix.ofRowList [[2;-3;1]; [1;-1;-2]; [3;1;-1]]
    let b: Vector<double> = DenseVector.ofList [7;-2;0]
    let expected: Vector<double> = DenseVector.ofList [1;-1;2]
    let actual = FEMS.Core.Solve.staticEquilibrum(A, b)
    
    Assert.True(expected.Subtract(actual).Maximum() < 0.000001)


