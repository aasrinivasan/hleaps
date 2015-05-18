# Author: Arun Srinivasan
# Email: ara.arun@gmail.com


# Test for: basic inputs
# Expected: hleaps should run without error and produce the 3 best models for each number of terms
#          adj r squared is the criteria
# Observed: hleaps correctly finds the best models and the adj r squared match the outputs from LM

set.seed(21564)
x1 = rnorm(100)
x2 = rnorm(100)
x3 = rnorm(100)
x4 = rnorm(100)
y = rnorm(100)
fr = cbind(x1,x2,x3,x4,y)
outStandard = hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE)
outStandardNoAlt = hleaps(y~(x1+x2+x3+x4)^2, altOut = FALSE)
outLM3 = lm(y~x3) # MATCH
outLM34 = lm(y~x3+x4) # MATCH
outLM23 = lm(y~x2+x3+x2:x3) # MATCH

# Test for: basic inputs factors
# Expected: hleaps should run without error and produce the 3 best models for each number of terms
#          adj r squared is the criteria
# Observed: hleaps correctly finds the best models and the adj r squared match the outputs from LM

f1 = c(rep("1",10), rep("2",30), rep("3",40))
f2 = c(rep("1",40), rep("2",0), rep("3",40))
f3 = c(rep("1",25), rep("2",30), rep("3",25))
y2 = c(rep("1", 30), rep("2",10), rep("3", 40))

outFactors = hleaps(y2 ~ (f1+f2+f3)^2, altOut=TRUE)
outFactorsLM1 = lm(y2~f1) # MATCH
outFactorsLM13 = lm(y2~f1+f3) # MATCH

# Test for: basic inputs combination of factors and non-factors
# Expected: hleaps should run wihtout error and produce the requested best models for each number
#           of terms with adj r squared as the criteria.
# Observed: hleaps correctly finds the best models and the adj r squared match the lm outputs.

f4 = c(rep("1",25), rep("2",25), rep("3",50))

outMixed = hleaps(y ~ (x1+x2+x3+f4)^2, altOut = TRUE)
outMixedLM = lm(y ~ x3+f4) # MATCH
outMixedLM34Int = lm(y ~ x3+f4+x3:f4) # MATCH

# Test for: colinearity check
# Expected: hleaps should run and there should be no models where there are both x1 and x5 or the
#          interaction
# Observed: hleaps correctly finds the best models and the colinearity check ensures x1 and x5 are
#          not in the same model

x5 = x1
outColin = hleaps(y~(x1+x2+x3+x5)^2, altOut=TRUE)

# Test for: NA in data
# Expected: In the first case, hleaps should run normally, as NA.Action defautls to 
#           na.omit. In case 2, na.fail should cause and error
# Observed: hleaps correctly runs in case one, then hleaps errors when na.action is set
#          to na.fail

library(testthat)
x6 = x4
x6[1] = NA

outNAOmit = hleaps(y ~ (x1 + x2 + x3 + x6)^2, altOut=TRUE)
lm(y~x6) # MATCH
expect_error(hleaps(y ~ (x1 + x2 + x3 + x6)^2, altOut=TRUE, na.action="na.fail"))

# Test for: fewer observations than possible interactions
# Expected: hleaps should not error out and produce results identical to those 
#          from lm
# Observed: hleaps correctly produces results identical to lm

x7 = rnorm(10)
x8 = rnorm(10)
x9 = rnorm(10)
x10 = rnorm(10)
y3 = rnorm(10)

outDiff = hleaps(y3~(x7+x8+x9+x10)^3, altOut=TRUE)
outDiffLM = lm(y3~x7)

# Test for: no intercept
# Expected: hleaps should drop the intercept and calculate criteria
# Observed: hleaps drops the intercept and the adj r squared calculations
#          match what is returned lm

outNoInt = hleaps(y~0+(x1+x2+x3+x4)^2, altOut=TRUE)
outNoIntLM = lm(y~0+x3)

# Test for: incorrect timeOut input. Case 1: negative time out
# Expected outcome: hleaps should show and error message that time
# Test for: Incorrect timeOut variable. Case 2: input is incorrect type (non-numeric)
# Expected outcome: Should trigger an error message stating that time must be a numeric
#                  input.
# Observed outcome: Both displayed the correct error

outTimeOutChar = hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, timeOut ="v")
outTimeOutNeg = hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, timeOut = -5)


# Test for: incorrect mem input. Case 1: negative mem
# Expected outcome: hleaps should show and error message that time and return.
# Test for: Incorrect mem input Case 2: input is incorrect type (non-numeric)
# Expected outcome: Should trigger an error message stating that mem must be a numeric
#                  input and return.
# Test for: incorrect mem input. Case 3: less than 1 mem
# Expected outcome: hleaps should show the error message and return
# Observed outcome: All displayed the correct error and all cases ran

outMemChar = hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, mem ="v")
outMemNeg = hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, mem = -5)

# Test for: improper method selections. The user inputs an unsupported method: "SS" and TRUE
# Expected outcome: error message that the method is not supported and default to adjr2.

outCrit <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, 
                  altOut=TRUE,  method = "SS")
outCritTF <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, 
                    altOut=TRUE,  method = TRUE)

# Test for: max < 0
# Expected outcome: prompts user that maxSize must be greater than zero, returns to default value of NULL
# Test for: min < 0
# Expected outcome: prompts user that minSize must be at least 1 and then runs with default values of NULL.
# Observed outcome: The error message correctly is shown in both cases.

outNegMax <-hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE,  maxSize = -1)
outNegMin <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE, minSize = -2)

# Test for: minSize > maxSize
# Expected outcome: prompts user that maxSize cannot be less than minSize. Returns to default values.
# Observed outcome: The error message correctly is shown and the defaults are run

outMaxMin <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE,  maxSize = 1, minSize = 3)


# Test for: minSize > the possible size
# Expected outcome: Error message should state that minSize > than number of possible terms. Then
#                   should run with default to 1.

outLargeMin = hleaps(y~(x1 + x2 + x3 + x4)^2, altOut=TRUE, minSize = 10000)


# Test for: nBest = 0
# Expected outcome: Error message should state that nBest should be at least 1, and default to
#                   3.
# Test for: nBest is negative
# Expected outcome: Error message should state that nBest should be at least 1, and default to
#                   3.
# Test for: nBest is an incorrect type - non-numeric
# Expected outcome: Error message should state that nBest should be numeric and should default to 3.
# Observed outcome: All cases produce correct error message.

outBadnBest <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, nBest = "2c")
out0nBest <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, nBest = 0)
outnBestNeg <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, nBest = -1)

# Test for: minSize is an incorrect type - non-numeric
# Expected outcome: Error message should state that minSize should be numeric and should default to 1.
# Test for: maxSize is an incorrect type - non-numeric
# Expected outcome: Error message should state that maxSize should be numeric and should default to num.terms.
# Observed outcome: Correct error message was shown and correctly defaulted.

outBadminSize <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut=TRUE, minSize = "3c")
outBadmaxSize <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut=TRUE, maxSize = "4cdd")

# Test for: altOut is an incorrect type - non-logical
# Expected outcome: Error message should state that altout should be numeric and should default to FALSE.
# Observed outcome: Error message correctly stated that altOut must be logical and defaulted to FALSE

outBadAltOut <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut="r")

# Test for: n < p
# Expected outcome: The algorithm should run with the correct dropping of terms taken care of by code
# Observed outcome: algorithm correctly runs with no error

x12 = rnorm(10)
x13 = rnorm(10)
x14 = rnorm(10)
x15 = rnorm(10)
y3 = rnorm(10)

outLess <- hleaps(y3~(x12+x13+x14+x15)^3, altOut = TRUE)

# Test for: lm output - subset
# Expected outcome: The test with subset 1:10 should be the same as the test with only rows
#                   1:10 in the dataframe
# Observed outcome: the subset command works correctly. The LM outputs match the outputs from
#                   the hleaps outcome.

fr2 = fr[1:10,]
outSub = hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE, data=as.data.frame(fr), subset = 1:10)
outSubDF = hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE, data=as.data.frame(fr2))

outLMSub = lm(y~x2, data=as.data.frame(fr), subset=1:10) # MATCH
expect_equal(outSub$modelInfo, outSubDF$modelInfo) # PASS

# Test for: lm output - subset
# Expected outcome: The test with subset 1:10 should be the same as the test with only rows
#                   1:10 in the dataframe
# Observed outcome: the subset command works correctly. The LM outputs match the outputs from
#                   the hleaps outcome.

x18 = runif(n=100, 1, 2)
outWeights = hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE, data=as.data.frame(fr), weights = x18)

outLMSub = lm(y~x2, data=as.data.frame(fr), subset=1:10) # MATCH
expect_equal(outSub$modelInfo, outSubDF$modelInfo) # PASS

# Test for: lm output - weights
# Expected outcome: The test with weights should match the lm command with weights for that model
# Observed outcome: The weights command does not seem to be working as the outputs are remaining the
#                   same with weights added and removed. The weights seem to be making it to hpmodsetup
#                   correctly.

x16 = runif(n=100, 1, 2)
outWeights = hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE, weights = x18)
outLMWeight = lm(y~x4, data=as.data.frame(fr), weights = x18) # MATCH

# Test for: lm output - offset
# Expected outcome: The test with offset should match the lm command with offset
# Observed outcome: The offset commmand does not work correctly. The offset seems to be correctly
#                   be passed into the hpmodsetup; however, the offset does not seem to 
#                   be done

x18 = runif(n=100, 1, 2)
outOff = hleaps(y~(x1+x2+x3+x4)^2, altOut = TRUE, offset = x18)

outLMSub = lm(y~x4, data=as.data.frame(fr), offset = x18) # FAIL