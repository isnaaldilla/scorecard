# IV
#
# This is an function to generate Information Value
# Information Value (IV) is predictive strength of each characteristic / variable

IV <- function(data,gbi,wt,vardisc,varcont){

  data$good[data$gbi==0] <- data$wt
  data$good[data$gbi==1] <- 0
  data$bad[data$gbi==1] <- data$wt
  data$bad[data$gbi==0] <- 0

  gtot <- sum(data$good,na.rm=TRUE)
  btot <- sum(data$bad,na.rm=TRUE)
  gbtot <- gtot + btot

  #discrete
  IV <- array(1:length(vardisc))
  char_analysis_disc <- NULL
  for(i in 1:length(vardisc))
  {
    a = data[,paste(vardisc[i])]
    bad <- aggregate(list(bad=data$bad), by=list(val=a), FUN=sum)
    good <- aggregate(list(good=data$good), by=list(val=a), FUN=sum)
    goodbad <- merge(good,bad,by="val")

    b <- ((goodbad$good/gtot)-(goodbad$bad/btot))*log((goodbad$good/gtot)/(goodbad$bad/btot))
    b[is.infinite(b)]<-NA
    IV[i] <- sum(b,na.rm=TRUE)
    IV[i] <- round(IV[i], digits = 4)
    goodbad$var <- vardisc[i]
    goodbad <- rename(goodbad, c(val="mylo"))
    goodbad$myhi <- goodbad$mylo
    char_analysis_disc <- rbind(char_analysis_disc, goodbad)
  }
  char_analysis_disc <- char_analysis_disc[order(char_analysis_disc$var,char_analysis_disc$mylo), ]
  ddply(char_analysis_disc, .(var), transform, .no = seq_along(var))
  char_analysis_disc$no <- with(char_analysis_disc, ave(var, var, FUN = seq_along))
  char_analysis_disc <- char_analysis_disc[c("no","var", "mylo", "myhi", "good", "bad")]
  IV_disc <- data.frame(vardisc,IV)

  rm(bad)
  rm(good)
  rm(goodbad)
  rm(a)
  rm(b)
  rm(IV)

  #continous
  IV <- array(1:length(varcont))
  char_analysis_cont <- NULL
  for(i in 1:length(varcont))
  {
    a = data[,paste(varcont[i])]
    x <- aggregate(list(tot=data$wt), by=list(var=a), FUN=sum)
    y <- aggregate(list(good=data$good), by=list(var=a), FUN=sum)
    z <- aggregate(list(bad=data$bad), by=list(var=a), FUN=sum)
    summary <- merge(x, y, by="var")
    summary <- merge(summary,z,by = "var")

    IV_cont_grp <- NULL
    tot <- 0
    mylo <- 0
    myhi <- 0
    good <- 0
    bad <- 0
    prior <- 0.05*gbtot
    j <- 1
    while (j <= nrow(summary))
    {
      row <- summary[j,]
      if (j == 1)
      {
        mylo <- row$var
        tot <- row$tot
        myhi <- row$var
        good <- row$good
        bad <- row$bad
        if (tot >= prior)
        {
          IV_cont_grp <- rbind(IV_cont_grp, data.frame(mylo,myhi,good, bad, tot))
          summary <- summary[j+1:nrow(summary),]
          summary <- summary[complete.cases(summary), ]
          j <- 1
        }
        else
        {
          j <- j+1
        }
      }
      else if ( j > 1 & j != nrow(summary) & tot < prior)
      {
        mylo <- mylo
        tot <- tot + row$tot
        myhi <- row$var
        good <- good + row$good
        bad <- bad + row$bad
        if (tot < prior)
        {
          j <- j+1
        }
        else
        {
          IV_cont_grp <- rbind(IV_cont_grp, data.frame(mylo,myhi,good, bad, tot))
          summary <- summary[j+1:nrow(summary),]
          summary <- summary[complete.cases(summary), ]
          j <- 1
        }
      }
      else if (j == nrow(summary))
      {
        if (tot < prior)
        {
          mylo <- mylo
          tot <- tot + row$tot
          myhi <- row$var
          good <- good + row$good
          bad <- bad + row$bad
        }
        else
        {
          mylo <- row$var
          tot <- row$tot
          myhi <- row$var
          good <- row$good
          bad <- row$bad
        }
        IV_cont_grp <- rbind(IV_cont_grp, data.frame(mylo,myhi,good, bad, tot))
        j <- j+1
      }
      else
      {
        IV_cont_grp <- rbind(IV_cont_grp, data.frame(mylo,myhi,good, bad, tot))
        summary <- summary[j+1:nrow(summary),]
        summary <- summary[complete.cases(summary), ]
        j <- 1
      }
    }
    IV_cont_grp$var <- varcont[i]
    char_analysis_cont <- rbind(char_analysis_cont, IV_cont_grp)
    b <- ((IV_cont_grp$good/gtot)-(IV_cont_grp$bad/btot))*log((IV_cont_grp$good/gtot)/(IV_cont_grp$bad/btot))
    b[is.infinite(b)]<-NA
    IV[i] <- sum(b,na.rm=TRUE)
    IV[i] <- round(IV[i], digits = 4)
  }
  char_analysis_cont <- char_analysis_cont[order(char_analysis_cont$var,char_analysis_cont$mylo), ]
  ddply(char_analysis_cont, .(var), transform, .no = seq_along(var))
  char_analysis_cont$no <- with(char_analysis_cont, ave(var, var, FUN = seq_along))
  char_analysis_cont <- char_analysis_cont[c("no","var", "mylo", "myhi", "good", "bad")]
  IV_cont <- data.frame(varcont,IV)

  # Join IV disc+cont
  IV_disc <- rename(IV_disc, c(vardisc="var"))
  IV_cont <- rename(IV_cont, c(varcont="var"))
  IV_all <- rbind(IV_disc, IV_cont)
  IV_all <- IV_all[order(IV_all$IV, decreasing = TRUE), ]

  # Join char analysis disc+cont
  char_analysis_all <- rbind(char_analysis_cont, char_analysis_disc)

  char_analysis_all$gp <- char_analysis_all$good / gtot
  char_analysis_all$bp <- char_analysis_all$bad / btot
  char_analysis_all$tp <- (char_analysis_all$good+char_analysis_all$bad) / gbtot
  char_analysis_all$rowrate <- char_analysis_all$bad / (char_analysis_all$good+char_analysis_all$bad)

  char_analysis_all$gp <- round(char_analysis_all$gp*100, digits = 2)
  char_analysis_all$bp <- round(char_analysis_all$bp*100, digits = 2)
  char_analysis_all$tp <- round(char_analysis_all$tp*100, digits = 2)
  char_analysis_all$rowrate <- round(char_analysis_all$rowrate*100, digits = 2)

  # Sort variable in char analysis by IV(decreasing)
  char_analysis_all <- merge(char_analysis_all, IV_all, by="var")
  char_analysis_all <- char_analysis_all[order(-char_analysis_all$IV,as.numeric(char_analysis_all$no)), ]
  char_analysis_all <- char_analysis_all[c("no", "var", "mylo", "myhi", "gp", "bp","tp","rowrate")]

  rm(char_analysis_cont,char_analysis_disc,IV_cont,IV_cont_grp,IV_disc,row,summary,x,y,z)
  rm(a,b,bad,btot,gbtot,good,gtot,i,IV,j,myhi,mylo,prior,tot)

  return(char_analysis_all)

}



