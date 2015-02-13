
doInstall <- TRUE  # Change to FALSE if you don't want packages installed.
toInstall <- c("ggplot2", "reshape2","OutbreakTools","ape","grid","gridExtra")
if(doInstall){install.packages(toInstall, repos = "http://cran.us.r-project.org")}
lapply(toInstall, library, character.only = TRUE)


# Generate a random matrix
# This can be any type of numeric matrix or characters
nRow <- 9
nCol <- 3
 myData.1 <- matrix(rnorm(nRow * nCol), ncol = nCol)
rownames(myData.1) <- letters[1:nRow]
colnames(myData.1) <- LETTERS[1:nCol]
# a discrete dataset of characters
nRow <- 9
nCol <- 2
char1<-c("A","B","C","C","A","C","B","A","B","B","B","C","B","A","C","B","A","B")
myData.2 <- matrix(char1, ncol = nCol)
rownames(myData.2) <- letters[1:nRow]
colnames(myData.2) <- LETTERS[1:nCol]
# a DNAbin dataset
cat("((((a:0.05,b:0.1):0.02,((c:0.03,d:0.04):0.05,e:0.05):0.3):0.02,f:0.06):0.01,(g:0.3,(h:0.02,i:0.01):0.01):0.02):0.01;",
   file = "letter.tre", sep = "\n")
tree.letters <- read.tree("letter.tre")

x <- structure(c("a", "j", "c", "e", "AATTCAATGCTCGGGAAGCAAGGAAAGCTGGGGACCAACTTCTCTTGGAGACATGAGCTTAGTGCAGTTAGATCGGAAGAGCA", "AATTCCTAAAACACCAATCAAGTTGGTGTTGCTAATTTCAACACCAACTTGTTGATCTTCACGTTCACAACCGTCTTCACGTT", "AATTCACCACCACCACTAGCATACCATCCACCTCCATCACCACCACCGGTTAAGATCGGAAGAGCACACTCTGAACTCCAGTC", "AATTCTATTGGTCATCACAATGGTGGTCCGTGGCTCACGTGCGTTCCTTGTGCAGGTCAACAGGTCAAGTTAAGATCGGAAGA"), .Dim = c(4L, 2L))
y <- t(sapply(strsplit(x[,2],""), tolower))
rownames(y) <- x[,1]
dna_y<-as.DNAbin(y) 

#list of dataframes and DNAbin
ff <- list(myData.1, myData.2,dna_y)
#different widths of the plots as porportions
size<-c(0.2,0.2,0.15,0.45)
titles<-c("Serology","Age Group","Polymerase")
palettes<-list(c("#D55E00", "#CC79A7"),c("#56B4E9", "#F0E442", "#D55E00"),NULL)
plot_ggphydata(tree.letters,ff,size,titles,palettes)
palettes<-list(c("#D55E00", "#CC79A7"),c("#56B4E9", "#F0E442", "#D55E00"),c("red","blue","orange","green"))
plot_ggphydata(tree.letters,ff,size,titles,palettes)

titles<-c("Serology","Age Group")
palettes<-list(c("#D55E00", "#CC79A7"),c("#56B4E9", "#F0E442", "#D55E00"))
ff <- list(myData.1, myData.2)
size<-c(0.2,0.2,0.15)
plot_ggphydata(tree.letters,ff,size,titles,palettes)

# randomly create age groups
# x <- as.data.frame(sample( LETTERS[1:3], 81, replace=TRUE, prob=c(0.1, 0.2, 0.7) ))
# zebov_tree <- read.tree("gires_zebov_2014.bayes.tree")
# rownames(x)<-zebov_tree$tip.label
# x<-list(x)

#size<-c(0.7,0.3)
#plot_ggphydata(zebov_tree,x,size)

