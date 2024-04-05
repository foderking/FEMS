module FEMS.Core.Elements

open MathNet.Numerics.LinearAlgebra

/// Specifies type of force acting on beam element
type BeamForceType =
    ///  Force is liner around the length of beam
    | LinearlyDistributed  of double
    ///  Force is constant around the length of beam. Zero at left end, maximum at right end
    | UniformlyDistributed of double
    // (F1,M1,F2,M2) f and m represents forces at moments at left end (1) and right end (2)
    | BoundaryForce  of double*double*double*double

type Beam1d(L: double, n: double, EI: double, beamType: BeamForceType) =
    let l = L/n;

    interface Element with
        member this.localForce     = 
            match beamType with
            | LinearlyDistributed(q) ->
                (q*l/60.0) * DenseVector.ofList [9.0;2.0*l;21.0;-3.0*l]

            | UniformlyDistributed(q) ->
                (q*l/12.0) * DenseVector.ofList [6.0;l;6.0;-l]

            | BoundaryForce(F1,M1,F2,M2) ->
                DenseVector.ofList [F1;M1;F2;M2]

        member this.localStiffNess = (EI/l**3) * DenseMatrix.ofRowList [
            [12.0 ; 6.0*l   ;-12   ; 6.0*l   ];
            [6.0*l; 4.0*l**2;-6.0*l; 2.0*l**2];
            [-12.0;-6.0*l   ; 12.0 ;-6.0*l   ];
            [6.0*l; 2.0*l**2;-6.0*l; 4.0*l**2]
        ]

    member this.totalNodes = 4 + 2*(int n-1)

    member this.shape(x: double): VectorD =
        DenseVector.ofList [
            1.0 - 3.0*(x/l)**2 + 2.0*(x/l)**3;
            x - 2.0*x**2/l + x**3/l**2;
            3.0*(x/l)**2 - 2.0*(x/l)**3;
            -x*x/l + x*x*x/l/l;
        ]
