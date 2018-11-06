#tools for working with self made montecarlo distribution
#we know that he p value for something single tailed is equal to 1 minus the cumulative distribution of that variable
#We can calculate this by using a handy dandy R tool ecdf to calculate the empirical cumulative density function
library(randomNames)
library(stats)

sample_numb <- 10000
df <- data.frame(dice_value = rnorm(sample_numb,sd=1),
                 name = randomNames(n=sample_numb)
                 )

cdf <- ecdf(df$dice_value)
p_value <- 1 - cdf(3)
p_value



#for creating the name co occurrence matrix later, this is for the more complicated example
number_of_gene_types <- 430
number_of_plasmids <- 8000
totaldata <- 60000

plasmid_names <- randomNames(n=number_of_plasmids)
gene_names <- randomNames(n= number_of_gene_types)

df <- data.frame(gene= sample(gene_names, totaldata, replace = TRUE),
                 plasmid = sample(plasmid_names, totaldata, replace = TRUE),
                 location = sample(1e6, totaldata),
                 circumference = 1e6,
                 all = "total_count")
library(reshape2)
dat2 <- melt(df)
w <- dcast(dat2, gene~plasmid)


plasmids_containing<- as.data.frame(w[, -1])
plasmids_containing[is.na(plasmids_containing)] <- 0
plasmids_containing <- apply(plasmids_containing, 2,  function(plasmids_containing) as.numeric(plasmids_containing > 0))  #recode as 0/1
occur_for_each <- rowSums(plasmids_containing)
occur_for_each <- as.data.frame(occur_for_each)
rownames(occur_for_each)<- w$gene



w <- dcast(dat2, gene~plasmid)
x <- as.matrix(w[,-1])
x[is.na(x)] <- 0
x <- apply(x, 2,  function(x) as.numeric(x > 0))  #recode as 0/1
v <- x %*% t(x)                                   #the magic matrix 
dimnames(v) <- list(w[, 1], w[,1])                #name the dimensions

#occur_for_each <- diag(v)

View(diag(v))

matrix1<- do.call("cbind",(rep(occur_for_each,length(rownames(occur_for_each)))))
matrix2<- do.call("rbind",(rep(occur_for_each,length(rownames(occur_for_each)))))
divisor <- matrix1 + matrix2
vx2 <- v*2
dicemat <- vx2/divisor
library(purrr)
#now get the list of all of the values in upper or lower
#triangle so we can construct our p values from our list of dice values
length(as_vector(dicemat))

dice_list <- as_vector(dicemat[upper.tri(dicemat)])
p_val_function <- ecdf(dice_list)
p_list <- as.data.frame(1-p_val_function(dice_list))
colnames(p_list)<- "values"
significant_vals <- p_list[p_list$values <= 0.05,]


