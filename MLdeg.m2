restart;
needsPackage("PHCpack");
--given monomial list, compute ML degree
--Uses Birchs theorem
MLdeg = (exps,coeffs) -> (
    --exps is A-matrix, coeffs is list of coeffs

    n := #coeffs;
    U := (for i in 1..n list ( random (1,10^10)));
    N := sum U;
    R := QQ[vars {1..(numgens target exps)}];
    A := substitute(exps,R);
    X := for i in 0..(n-1) list (
        N * coeffs_i * R_(entries exps_i)
    );
    system := entries (A*transpose((matrix {X}) - substitute(matrix {U},R)))_0;
    --return (system,degree saturate(ideal system,ideal product gens R));
    return degree saturate(ideal system,ideal product gens R));

    --pointList := solveSystem system;
    --toReturn := 0;
    --for p in pointList do (
        --if norm(2,p)>10e-8 then toReturn = toReturn + 1
        --else print "hihihi";
    --);
    --return (toReturn,pointList,system);
)

m = matrix {{1,1,1,1},{0,1,2,3}};
c = {1,1,1,1};
m = matrix {{1,1,1,1,1,1,1,1,1},{1,1,1,0,0,0,0,0,0},{0,0,0,1,1,1,0,0,0},{0,0,0,0,0,0,1,1,1},{1,0,0,1,0,0,1,0,0},{0,1,0,0,1,0,0,1,0},{0,0,1,0,0,1,0,0,1}};
c = {1,2,3,1,2,4,5,6,7};
m = matrix {{1,1,1,1,1,1},{0,1,2,0,1,0},{0,0,0,1,1,2}};
c = {1,2,1,1,3,1};
print m
asdf = MLdeg(m,c)
print asdf
