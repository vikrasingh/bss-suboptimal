load("C:/Users/vikra/Downloads/Leukemia.RData") # load the data from the file
write.csv(Leukemia$x), file = "xMatrix.txt", col.names= FALSE)
write.csv(Leukemia$y, file = "yArray.txt", col.names= FALSE)