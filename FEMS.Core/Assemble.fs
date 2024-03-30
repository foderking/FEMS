module FEMS.Core.Assemble

open MathNet.Numerics.LinearAlgebra

let assemble(localMatrix: MatrixD, globalMatrix: MatrixD, guide: List<int>): unit =
    let n = guide.Length

    for local_i in 0..n-1 do
        for local_j in 0..n-1 do
            let global_i = guide[local_i]
            let global_j = guide[local_j]
            globalMatrix[global_i,global_j] <- globalMatrix[global_i,global_j] + localMatrix[local_i,local_j]
// TODO: use element by element techniques to save storage