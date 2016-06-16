restart
installPackage"Bertini"
--%%--numU is the number of u's there will be. This is the same as the number of monomials.
--%%%%--The code is set to work for A matrices with 2 rows. 
--%%%%--If you have more than 2 rows, then the ourMonomials line and xList need to be adjusted. 
numU=4
uList=for i to numU-1 list "u"|i
cList=for i to numU-1 list "c"|i
R=QQ[s,t,mu]**QQ[uList]**QQ[cList]
AMatrix=matrix{{1,1,1,1},{0,1,2,3}}
uList=for i to numU-1 list value("u"|i)
cList=for i to numU-1 list value("c"|i)
xList={s,t}
ourMonomials=for i in entries transpose AMatrix list (xList_0)^(i_0)*(xList_1)^(i_1)

parameterization=for i to #cList-1 list cList_i*ourMonomials_i

leftSide=(sum uList) *sub(AMatrix,R) *transpose matrix {parameterization}
rightSide=AMatrix *transpose matrix {uList}

likelihoodEquations=flatten entries (leftSide-rightSide)
cWin={1,3,3,1}
cSub=for i to #cList-1 list cList_i=>cWin_i
uFix=for i to #uList-1 list random RR
uSub=for i to #uList-1 list uList_i=>uFix_i

bertiniZeroDimSolve(likelihoodEquations,
    AffVariableGroup=>xList,
    B'Constants=>uSub|cSub
    )
moveB'File(storeBM2Files,"nonsingular_solutions","start")
writeParameterFile(storeBM2Files,cWin,NameParameterFile=>"start_parameters")

cStat=for i to #cList-1 list random RR
writeParameterFile(storeBM2Files,cStat)
makeB'InputFile(storeBM2Files,
    B'Polynomials=>likelihoodEquations,
    B'Configs=>{{ParameterHomotopy,2},{MPType,2}},
    B'Constants=>uSub,    
    AffVariableGroup=>xList,
    ParameterGroup=>{c0,c1,c2,c3}    
    )
runBertini(storeBM2Files)
finalPoint=importSolutionsFile(storeBM2Files)
