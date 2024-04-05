namespace FEMS.Core

open MathNet.Numerics.LinearAlgebra

type MatrixD = Matrix<double>
type VectorD = Vector<double>

type Element =
    abstract LocalForce    : VectorD
    abstract LocalStiffNess: MatrixD
