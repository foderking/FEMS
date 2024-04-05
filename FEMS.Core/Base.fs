namespace FEMS.Core

open MathNet.Numerics.LinearAlgebra

type MatrixD = Matrix<double>
type VectorD = Vector<double>

type Element =
    abstract localForce    : VectorD
    abstract localStiffNess: MatrixD

module Constants =
    let zerosVec4: VectorD = DenseVector.zero 4
    let zerosVec2: VectorD = DenseVector.zero 2

module Operations =
    // TODO: use element by element techniques to save storage
    let assembleMatrix(localMatrix: MatrixD, globalMatrix: MatrixD, guide: List<int>): unit =
        let n = guide.Length

        for local_i in 0..n-1 do
            for local_j in 0..n-1 do
                let global_i = guide[local_i]
                let global_j = guide[local_j]
                globalMatrix[global_i,global_j] <- globalMatrix[global_i,global_j] + localMatrix[local_i,local_j]

    let assembleVector(localVector: VectorD, globalVector: VectorD, guide: List<int>): unit =
        let n = guide.Length
        for local_i in 0..n-1 do
            let global_i = guide[local_i]
            globalVector[global_i] <- globalVector[global_i] + localVector[local_i]

    let enforceBC(globalMatrix: MatrixD, globalVector: VectorD, matrixConstr: List<List<int>>, vecConstr: List<List<int>>): unit =
        for [bcIndex; bcValue] in vecConstr do 
            globalVector[bcIndex] <- bcValue

        for [bcIndex; bcValue] in matrixConstr do 
            for i in 0..globalMatrix.ColumnCount-1 do
                globalMatrix[bcIndex, i] <- 0
            globalMatrix[bcIndex,bcIndex] <- 1
            globalVector[bcIndex] <- bcValue
