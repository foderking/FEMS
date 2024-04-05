module Test.Base

let looseEqual(a: double, b: double): bool =
    let p = (a-b)/a
    p < 1e-14

