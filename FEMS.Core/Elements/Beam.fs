module FEMS.Core.Elements

open MathNet.Numerics.LinearAlgebra

type BeamForceType =
    | LinearlyDistributed  of double
    | UniformlyDistributed of double
    | BoundaryForce  of double*double*double*double

type Beam1d(L: double, n: double, EI: double, beamType: BeamForceType) =
    let l = L/n;
    interface Element with
        member this.LocalForce     = 
            match beamType with
            | LinearlyDistributed(q) ->
                (q*l/60.0) * DenseVector.ofList [9.0;2.0*l;21.0;-3.0*l]

            | UniformlyDistributed(q) ->
                (q*l/12.0) * DenseVector.ofList [6.0;l;6.0;-l]

            | BoundaryForce(F1,M1,F2,M2) ->
                DenseVector.ofList [F1;M1;F2;M2]

        member this.LocalStiffNess = (EI/l**3) * DenseMatrix.ofRowList [
            [12.0 ; 6.0*l   ;-12   ; 6.0*l   ];
            [6.0*l; 4.0*l**2;-6.0*l; 2.0*l**2];
            [-12.0;-6.0*l   ; 12.0 ;-6.0*l   ];
            [6.0*l; 2.0*l**2;-6.0*l; 4.0*l**2]
        ]
