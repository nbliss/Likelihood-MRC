import numpy as np
#------------------------------------------------------------------------------#
def iterativeProportionalScaling(A,u,eps,verbose=True):
    """
    Python / NumPy implementation of iterative proportional scaling, which is
    Algorithm 2.1.9 in Lectures on Algebraic Statistics,
    https://math.berkeley.edu/~bernd/owl.pdf
    NOTE: assumes column sums of A are all the same.
    """
    global intermediateArray
    dataType = "float128"
    #dataType = "float64"
    A = np.array(A)
    d,k = A.shape
    u = np.array(u,dtype=dataType)
    u.shape = (k,1)
    a = np.sum(A[:,0])

    # initialize v
    v = np.full((k,1),np.linalg.norm(u,ord=1)/k)

    intermediateArray = np.empty((d,k),dtype=dataType) #we'll reuse it for speed
    def applyPhi(x):
        """
        Apply each column of A to x as exponents, then take column products.
        """
        intermediateArray[:,:] = x
        applied = np.power(intermediateArray,A)
        applied = np.prod(applied,axis=0)
        applied.shape = (k,1)
        return applied

    Au = np.dot(A,u)
    for i in xrange(90000):
        #print '\r'+str(i),
        Av = np.dot(A,v)
        error = np.linalg.norm(Av-Au,ord=1)
        #print error,
        if error<=eps:
            if verbose:print "\nSucceeded in",i,"iterations, with residual",error
            return v
        v = v * np.power(applyPhi(np.divide(Au,Av)),1.0/a)
    else:
        if verbose:print "\nGave up after",i,"iterations, with residual",error
        return v
    return v

#------------------------------------------------------------------------------#
def ratNormScroll(degList,dataType):
    # return the matrix for the RNS w/ degrees in degList
    blocks = []
    numBlocks = len(degList)
    maxVal = max(degList[:-1]+[degList[-1]-1]) + 1
    for i in xrange(len(degList)):
        thisBlock = np.zeros((numBlocks,degList[i]+1),dtype=dataType)
        # ith row = ones
        for j in xrange(degList[i]+1):
            thisBlock[i][j]=1
            thisBlock[numBlocks-1][j] = j
        blocks.append(thisBlock)
    toReturn = np.bmat(blocks)
    lastRow = maxVal - np.sum(toReturn,axis=0)
    return np.bmat([[toReturn],[lastRow]])


#------------------------------------------------------------------------------#
testing = False
if testing:
    np.set_printoptions(linewidth=128,precision=9)
    degList = [20 for i in xrange(15)] # ni's are 20, d-1 = 15
    u = np.random.randint(1,20,size=sum(degList) + len(degList))
    A = ratNormScroll(degList,"int")
    eps = 10e-12
    from timeit import timeit
    toTime = lambda:iterativeProportionalScaling(A,u,eps,verbose=False)
    numTrials = 7
    averageTime = timeit(toTime,number=numTrials)/numTrials
    print "Average time from",numTrials,"trials:",averageTime
#------------------------------------------------------------------------------#
