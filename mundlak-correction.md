# Mundlak correction

Input data: a dataframe with the following columns (assuming 3 covariates and 5 chr-PRS. in our analyses, cov's were 21, and PRSs 22)

`FID IID PHENO cov1 cov2 cov3 PRS1 PRS2 PRS3 PRS4 PRS5`

Then, the following code adjusts the PRS and adds the panel mean as a covariate. 

```r

## mundlak correction
library(dplyr)

args <- commandArgs(trailingOnly=TRUE)

#locations and variable names, all dependent on how you scale-up.
datafile <- args[1]
foldindir <- args[2]
mergefile <- args[3]
wdir <- args[4]
pheno <- args[5]
stats_out <- args[6]


#format data described above
data <- read.table(datafile, header=TRUE, stringsAsFactors=FALSE)


#creating a vector to include the membership of each individual to a specific fold(cluster vector)


#make a vector that has the folds as numbers for the individuals
indicator <- data.frame("IID"= numeric(), "fold"=numeric())

for (fold in seq(1,10)) {
	file <- paste0(foldindir, "/test_",fold,".txt")
	df <- read.table(file, header=FALSE, stringsAsFactors=FALSE)

	indic_new <- data.frame("IID"= df, "fold"=rep(fold,nrow(df)))

	indicator <- rbind(indicator, indic_new)

}
colnames(indicator) <- c("IID","fold")

#add the fold indication vector to the data.
data <- inner_join(data, indicator, by="IID" )

print("matched with fold nr's")

# panel average = average per fold per prs, then average those averages. 
## so first those perfold per chr averages
chrmeans <- data.frame("fold"=numeric(),"chr"=numeric(),"mean"=numeric())

for (i in seq(1,10)){
	for (chr in seq(1,22)){
		current_col  <- which(colnames(data)==paste0(pheno,"_PRS",chr))

		data_col <- filter(data, fold==i)[current_col]

		mean <- sum(data_col)/nrow(data_col)

		current_add <- data.frame("fold"=i,"chr"=chr,"mean"=mean)

		chrmeans <- rbind(chrmeans,current_add)
		
	}

}

print("calculated per fold means")
#now calculate panel means
panelmeans <- data.frame("fold"=numeric(),"mean"=numeric())

for (i in seq(1,10)){
	data_col <- filter(chrmeans, fold==i)

	mean <- sum(data_col$mean)/nrow(data_col)

	current_add <- data.frame("fold"=i,"panelmean"=mean)
	panelmeans <- rbind(panelmeans,current_add)

}

print("calculated panel means")
## add the panelmeans to the data variable to create fold&panelmean vector

data <- inner_join(data,panelmeans,by="fold")

#create output dataframe
lastcol <- which(colnames(data)=="cov21")
newdata <- data[c("1":lastcol)]


## remove those means from the PRSs, add those to the output frame
for (chr in seq(1,22)){
	col <- data[[paste0(pheno,"_PRS",chr)]] - data["panelmean"]

	colnames(col) <- paste0(paste0(pheno,"_PRS",chr))

	newdata[paste0(pheno,"_PRS",chr)] <- col

}

newdata$panelmean <- data$panelmean

print("did mundlak adjustment")

write.table(newdata, mergefile, quote=FALSE, col.name=TRUE, row.name=FALSE, sep=" ")

print("did mundlak adjustment")

```