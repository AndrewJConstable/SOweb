###Blue whales

rm(list=ls())

# parameter file #

param <- data.frame(min= c(200000,100,0.0044, 1682, 0.7,0.25,0.25, 0.7, 0.1, 0.1, 8264),max= c(300000,150,0.1905,4130, 0.9 ,0.75,0.75, 0.9, 0.3, 0.5, 416667), base= c(250000,120,0.0391, 2917, 0.8, 0.5, 0.5, 0.8, 0.2, 0.3, 140599))
row.names(param) <- c('N','d','k','c','fw','pw','ew','fk', 'pk', 'ek', 'u')


labels <- rownames (param)

# Function to compute primary production based on a given parameter set (11 parameters)

computePP <- function(param_set){
  return ((((param_set[1]* param_set[2]* param_set[3]*(0.23* param_set[4]) * #amount of Fe consumed as prey in kg
               (( param_set[5]*   #proportion of Fe released in whale faeces
                  param_set[6]*   #proportion of Fe in whale faeces remaining in the photic zone
                  param_set[7])+  #bioavailability of  whale faecal iron
                 (param_set[8]*   #proportion of Fe released in krill faecal pellet
                  param_set[9]*   #proportion of Fe in krill faecal pellet remaining in the photic zone
                  param_set[10])))/   #bioavailability of  krill faecal iron
                  55.845)*   #to get mol of Fe
                  12.011* param_set[11])/  #to get mol of C based on C:Fe ratio of phytoplankton
                  2e+13) 
}

##Base value
basevalue <- computePP(param$base)
maxvalue <- computePP(param$max)
minvalue <- computePP(param$min)

# GENERATE THE TABLE WITH ALL PARAMETER SETS TO BE USED IN THE SENSI ANALYSIS

table <- data.frame(base=param$base)
names <- 'Base'
for(j in seq(1,2)){
  for (i in seq(1,11)){    
    temp <- param$base
    # REPLACE i-th PARAMETER BY ITS MIN (j=1) or MAX (j=) VALUE IN THE VECTOR OF MODEL PARAMETER
    temp[i] <- param[i,j]
    table <- data.frame(table, temp)
    # NAME THE COLUMN ACCORDING TO PARAMETER NAME + MIN / MAX
    names(table)[ncol(table)] <-  paste(rownames(param)[i],names(param)[j],sep='_')
  }
}
table

apply(table,2, computePP)
min <- apply(table,2, computePP)[2:12]
max <- apply(table,2, computePP)[13:23]
min
max

##how to plot?
library(ggplot2) 
d <- data.frame(labels=factor(labels,labels[order(abs(max-min))]),
                min=min,
                max=max,
                base= basevalue)

italic.text <- element_text(face = "italic", color = "black", size = 8)

(b <- ggplot(d, aes(labels, base, ymin = min, ymax=max)) + 
  geom_linerange(mapping=aes(ymin=min,ymax=base),size=3,color="grey70") + 
  geom_linerange(mapping=aes(ymin=base,ymax=max),size=3,color="grey40") + 
  coord_flip() + theme_bw() + theme(legend.position = "none") + xlab("") +
  ylab(expression(paste(Estimate~of~primary~production~(g~C~m^-2~yr^-1)))) + theme(axis.text.y = italic.text))


(c <- b + geom_hline(yintercept = basevalue, show_guide = FALSE,
                    linetype = "dashed", color = "black", size=0.5) + ggtitle('Blue whales')) 



basevalue
minvalue
maxvalue

