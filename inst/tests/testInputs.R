context("hleaps inputs")

test_that("test inputs", {
  set.seed(21435)
  y = rnorm(100)
  x1 = rnorm(100)
  x2 = rnorm(100)
  x3 = rnorm(100)
  x4 = rnorm(100)
  
  # Test for: standard correct inputs
  # Expected outcome: should not not error
  outG <- hleaps(y~(x1 + x2 + x3 + x4)^2,  altOut = TRUE,
                 method = "adjr2", nbest = 4, minSize = 2, maxSize = 5)
  
  
  # Test for: correct processing of factors in both greedy and non-greedy
  # Expected outcome: no errors in greedy and non-greedy
  f1 = c(rep("1",30), rep("2",25), rep("3",45))
  f2 = c(rep("1",12), rep("2",25), rep("3",63))
  f3 = c(rep("1",44), rep("2",11), rep("3",45))
  y = f1
  
  outFactors <- hleaps(y ~ f1+f2+f3,  altOut=TRUE)
  
  
  # Test for: incorrect timeOut input. Case 1: negative time out
  # Expected outcome: hleaps should show and error message that time
  # Test for: Incorrect timeOut variable. Case 2: input is incorrect type (non-numeric)
  # Expected outcome: Should trigger an error message stating that time must be a numeric
  #                  input. Should not run searches.
  
  
  outTimeLG <- hleaps(y~(x1 + x2 + x3 + x4)^2,  altOut = TRUE,
                      timeOut ="v")
  outTimeNG <- hleaps(y~(x1 + x2 + x3 + x4)^2,altOut = TRUE,
                      timeOut = -5)
  outTimeG <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE,
                     timeOut = Inf)
  
  expect_identical(outTimeLG$modelInfo, outTimeG$modelInfo)
  expect_identical(outTimeNG$modelInfo, outTimeG$modelInfo)
  
  # Test for: improper method selections. The user inputs an unsupported method: "SS" and TRUE
  # Expected outcome: error message that the method is not supported and default to adjR2.
  
  outCrit <- hleaps(y ~ (x1 + x2 + x3 + x4)^2,altOut=TRUE,  method = "SS")
  outCritTF <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE,  method = TRUE)
  outCritAdjR2 <- hleaps(y ~ (x1 + x2 + x3 + x4)^2,altOut=TRUE,  method = "adjR2")
  
  expect_identical(outCrit$modelInfo, outCritAdjR2$modelInfo)
  expect_identical(outCritTF$modelInfo, outCritAdjR2$modelInfo)
  
  # Test for: max < 0
  # Expected outcome: prompts user that maxSize must be greater than zero, returns to default value of NULL
  # Test for: min < 0
  # Expected outcome: prompts user that minSize must be at least 1 and then runs with default values of NULL.
  
  outNegMax <-hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE,  maxSize = -1)
  outNegMin <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE, minSize = -2)
  outCorrect <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE)
  
  expect_identical(outNegMax$modelInfo, outCorrect$modelInfo)
  expect_identical(outNegMin$modelInfo, outCorrect$modelInfo)
  
  # Test for: minSize > maxSize
  # Expected outcome: prompts user that maxSize cannot be less than minSize. Returns to default values.
  
  outMaxMin <- hleaps(y ~ (x1 + x2 + x3 + x4)^2,  altOut=TRUE,  maxSize = 1, minSize = 3)
  
  outMaxMinCorrect <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE,  minSize = 3)
  
  expect_identical(outMaxMin$modelInfo, outMaxMinCorrect$modelInfo)
  
  # Test for: minSize > that possible size
  # Expected outcome: Error message should state that minSize > than number of possible terms. Then
  #                   should run with default to 1.
  
  outLargeMin = hleaps(y~(x1 + x2 + x3 + x4)^2,altOut=TRUE, minSize = 10000)
  outLargeMinCorrect = hleaps(y~(x1 + x2 + x3 + x4)^2, altOut=TRUE,minSize = 1)
  
  expect_identical(outLargeMin$modelInfo, outLargeMinCorrect$modelInfo)
  
  # Test for: nBest = 0
  # Expected outcome: Error message should state that nBest should be at least 1, and default to
  #                   3.
  # Test for: nBest is negative
  # Expected outcome: Error message should state that nBest should be at least 1, and default to
  #                   3.
  # Test for: nBest is an incorrect type - non-numeric
  # Expected outcome: Error message should state that nBest should be numeric and should default to 3.
  
  outBadnBest <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, maxSize = 4, minSize = 3, nBest = "2c")
  out0nBest <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, nBest = 0)
  outnBestNeg <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, nBest = -1)
  outnBestCorrect <- hleaps(y~(x1 + x2 + x3 + x4)^2, , altOut = TRUE,nBest = 3)
  
  expect_identical(out0nBest$modelInfo, outnBestCorrect$modelInfo)
  expect_identical(outBadnBest$modelInfo, outnBestCorrect$modelInfo)
  expect_identical(outnBestNeg$modelInfo, outnBestCorrect$modelInfo)
  
  # Test for: minSize is an incorrect type - non-numeric
  # Expected outcome: Error message should state that minSize should be numeric and should default to 1.
  # Test for: maxSize is an incorrect type - non-numeric
  # Expected outcome: Error message should state that maxSize should be numeric and should default to num.terms.
  
  outBadminSize <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut =TRUE, minSize = "3c", nBest = 3)
  outBadmaxSize <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, maxSize = "4cdd", nBest = 3)
  outBadMaxMinCorrect <- hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, nBest = 3)
  
  expect_identical(outBadminSize$modelInfo, outBadMaxMinCorrect$modelInfo)
  expect_identical(outBadmaxSize$modelInfo, outBadMaxMinCorrect$modelInfo)
  
  
  
  # Test for: NA exists in the data when passed in.
  # Expected outcome: no error should occur when NA is in the data, autos to na.action = na.omit
  
  y = rnorm(900)
  x1 = rnorm(900)
  x1[2] = NA
  x2 = rnorm(900)
  x3 = gl(3,10,900)
  x4 = gl(10,3,900)
  
  outNAG <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE)
  
  # Test for: na.action inputs
  # Expected outcome: na.action set to na.fail should error, na.pass should error, na.exclude should continue
  #                   as usual.
  
  expect_error(hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, na.action = "na.pass"))
  expect_error(hleaps(y~(x1 + x2 + x3 + x4)^2, altOut = TRUE, na.action = "na.fail"))
  outNAExclude <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE, na.action = "na.exclude")
  outNAOmit <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE, na.action = "na.omit")
  # Test for: entire term is NA.
  # Expected outcome: hleaps will not run due to the NA term f1. 
  
  f1 = NA
  f2 = c(rep("1",12), rep("2",25), rep("3",63))
  f3 = c(rep("1",44), rep("2",11), rep("3",45))
  y = f1
  
  expect_error(hleaps(y ~ f1+f2+f3, altOut=TRUE))
  
  # Test for: column in a factor variable is NA
  # Expected outcome: default na.action is na.omit; therefore, the rows in the factor should 
  #                   be omitted. Search should proceed as usual after omit.
  
  f1 = c(rep("1",30), rep("2",25), rep("3",45))
  f2 = c(rep(NA,12), rep("2",25), rep("3",63))
  f3 = c(rep("1",44), rep("2",11), rep("3",45))
  y = f1
  
  outNAcolG <- hleaps(y ~ f1+f2+f3, altOut=TRUE)
  
  # Test for: data-frame as an input
  # Expected: the greedy search with data frame should be identical with the greedy search with
  #           data terms as inputs. Likewise for non-greedy.
  
  y = rnorm(900)
  x1 = rnorm(900)
  x2 = rnorm(900)
  x3 = gl(3,10,900)
  x4 = gl(10,3,900)
  fr <- data.frame(y, x1, x2, x3, x4)
  
  outDFG <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE, data = fr)
  outNoDFG <- hleaps(y ~ (x1 + x2 + x3 + x4)^2, altOut=TRUE)
  
  expect_identical(outDFG$modelInfo, outNoDFG$modelInfo)
  
  # Test for: n < p
  # Expected outcome: The algorithm should run with the correct dropping of terms taken care of by code
  
  x12 = rnorm(10)
  x13 = rnorm(10)
  x14 = rnorm(10)
  x15 = rnorm(10)
  y3 = rnorm(10)
  
  outLess <- hleaps(y3~(x12+x13+x14+x15)^3, altOut = TRUE)
  
})