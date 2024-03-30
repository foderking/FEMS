module Assemble

open Xunit
open MathNet.Numerics.LinearAlgebra

[<Fact>]
let ``basic test assembling global matrix`` () =
    let globalM = DenseMatrix.create 8 8 0.0
    FEMS.Core.Assemble.assemble( 
        DenseMatrix.ofRowList [[1;2;3;4];[5;6;7;8];[9;10;11;12];[13;14;15;16]], globalM, [3;0;1;4]
    )
    Assert.Equal(globalM, matrix [
        [6.0;7;0;5;8;0;0;0];[10;11;0;9;12;0;0;0];[0;0;0;0;0;0;0;0];[2;3;0;1;4;0;0;0];
        [14;15;0;13;16;0;0;0];[0;0;0;0;0;0;0;0];[0;0;0;0;0;0;0;0];[0;0;0;0;0;0;0;0]
    ])
