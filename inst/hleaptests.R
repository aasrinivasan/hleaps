# Initialize Variables
set.seed(21435)
y = rnorm(100)
x1 = rnorm(100)
x2 = rnorm(100)
x3 = rnorm(100)
x4 = rnorm(100)
x5 = x1

x11 = rnorm(100)
x11[2] = NA

f1 = c(rep("1",10), rep("2",30), rep("3",40))
f2 = c(rep("1",40), rep("2",0), rep("3",40))
f3 = c(rep("1",25), rep("2",30), rep("3",25))
y2 = c(rep("1", 30), rep("2",10), rep("3", 40))
f4 = c(rep("1",25), rep("2",25), rep("3",50))
f11 = NA
f21 = c(rep(NA,12), rep("2",5), rep("3",63))

fr <- data.frame(y, x1, x2, x3, x4)
fr2 = fr[1:10,]
x12 = rnorm(10)
x13 = rnorm(10)
x14 = rnorm(10)
x15 = rnorm(10)
y3 = rnorm(10)

x16 = runif(n=100, 1, 2)

automatedHleapsTests(y, x1, x2, x3, x4, x11, f1, f2, f3, y2, f11, f21
               , x12, x13, x14, x15, y3, fr2, x16, x5, f4)

remove(y, x1, x2, x3, x4, x11, f1, f2, f3, y2, f11, f21
        , x12, x13, x14, x15, y3, fr2, x16, x5, f4)