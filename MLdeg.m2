--------------------------------------------------------------------------------
----Function to calculate the maximum likelihood equations of a toric model-----
----using the equation given by Birch`s theorem. -------------------------------
-------* exps is the A-matrix of exponents.-------------------------------------
-------* coeffs is a list of scaling coefficients for the monomial map.---------
-------* U is the data vector.--------------------------------------------------
-------* Setting coeffRing allows you to choose the coefficient ring of the-----
-------- returned system.-------------------------------------------------------
--------------------------------------------------------------------------------
MLeqs = {coeffRing=>QQ} >> o -> (exps,coeffs,U) -> (
    n := #coeffs;
    R := (o.coeffRing)[vars {53..(52+(numgens target exps))}];
    A := substitute(exps,R); --Convert A matrix to correct ring
    N := sum U;  --Sample size
    U = transpose matrix {for i in U list i_R};
    -- X is the tuple of the (scaled) monomial map
    X := for i in 0..(n-1) list coeffs_i * R_(entries exps_i);
    X = transpose matrix {X};
    system := ideal (A*(N * X - U));
    return first entries gens system;
)

--------------------------------------------------------------------------------
----Constructs the monomial map from the given data.----------------------------
----Returns a function that applies the monomial map to a Point or a matrix.----
--------* exps is the A-matrix of exponents.------------------------------------
--------* coeffs is a list of scaling coefficients for the monomial map.--------
--------------------------------------------------------------------------------
needsPackage "NAGtypes";
MLmap = (exps,coeffs) -> (
    n := #coeffs;
    R := QQ[vars {53..(52+(numgens target exps))}];
    X := for i in 0..(n-1) list coeffs_i * R_(entries exps_i);
    return (pt -> first entries evaluate(matrix {X},pt));
)

isPosReal := pt -> (
    if not isRealPoint(pt) then return false;
    return all(coordinates pt,a->(realPart a)>0);
);

--------------------------------------------------------------------------------
----Using MLeqs, computes the maximum likelihood degree of a (scaled) toric-----
----model.----------------------------------------------------------------------
-------* exps is the A-matrix of exponents.-------------------------------------
-------* coeffs is a list of scaling coefficients for the monomial map.---------
--------------------------------------------------------------------------------
MLdeg = (exps,coeffs) -> (
    U := (for i in 1..#coeffs list (random(1,10^5)));
    system := MLeqs(exps,coeffs,U);
    return degree saturate(ideal system,product gens ring first system);
);

--------------------------------------------------------------------------------
----Uses PHCpack to perform a homotopy from-------------------------------------
------* the ML equations based on A-matrix exps, scaling cWin, and a random-----
------- data vector, to---------------------------------------------------------
------* the ML equations with the same A-matrix, scaling cStat, and with--------
 ------ data vector U-----------------------------------------------------------
--------------------------------------------------------------------------------
try needsPackage "PHCpack";
performPHCHomotopy = (exps,cStat,cWin,U) -> (
    startSyst := MLeqs(exps,cWin,U,coeffRing => CC);
    --return startSyst;
    sols := solveSystem(startSyst,Verbose=>false);
    --return sols;
    sols = select(sols,isPosReal);
    targetSyst := MLeqs(exps,cStat,U,coeffRing => CC);
    toReturn = (timing trackPaths(targetSyst,startSyst,sols,gamma=>1));
    --print (MLmap(exps,cStat))(toReturn#1#0);
    return toReturn;
);

--------------------------------------------------------------------------------
----Same as performPHCHomotopy, just with Bertini. In this case we add in the---
----homotopy parameter manually.------------------------------------------------
--------------------------------------------------------------------------------
try needsPackage "Bertini";
performBertiniHomotopy = (exps,cStat,cWin,U) -> (
    startSyst := MLeqs(exps,cWin,U,coeffRing => CC);
    --sols := select(bertiniZeroDimSolve(startSyst,Verbose=>false),isPosReal);
    sols := select(bertiniZeroDimSolve(startSyst,Verbose=>true),isPosReal);
    R = CC[t]**(ring first startSyst);
    targetSyst := MLeqs(exps,cStat,U,coeffRing => CC);
    targetSyst = (targetSyst / (a->(1-t)*substitute(a,R)));
    startSyst = (startSyst / (a->t*substitute(a,R)));
    toSolve := (for i in 0..(#startSyst-1) list (targetSyst_i+startSyst_i));
    return timing bertiniTrackHomotopy(t,toSolve,sols,Verbose=>false);
);

----Rational normal curve example-----------------------------------------------
--m = matrix{{1,1,1,1},{0,1,2,3}};
--cStat = {1,1,1,1};
--cWin = {1,3,3,1};

--print MLdeg(m,cWin)
--print MLdeg(m,cStat)
--print performPHCHomotopy(m,cStat,cWin,{2,3,5,7})
--print performBertiniHomotopy(m,cStat,cWin,{2,3,5,7})

