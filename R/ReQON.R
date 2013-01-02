ReQON <-
function(in_bam, out_bam, region, max_train = -1, SNP = "", 
		RefSeq = "", nerr = 2, nraf = 0.05, plotname = "", temp_files = 0){

  message("\nStart time: ", date(), "\n\n") 
  diagnostics <- list()
	
  # obtain training region from BAM unpacker
   .jpackage("ReQON")  
    message("Reading in bam file")
    s1 <- .jarray(c(paste("INPUT=", in_bam, sep=""), "OUTPUT=train_temp.txt", paste("REGION=", region, sep="")))       
    bup <- .jnew("org/renci/sequencing/util/BAMUnpacker")
    .jcall(bup, "V", "main",s1)
    dat <- read.table("train_temp.txt",nrows=max_train)  
    unlink("train_temp.txt")

  # data set-up
    if(min(dat[,5]) > 32){ 	  
      dat[,5] <- dat[,5] - 33
      dat[,6] <- dat[,6] - 33
    }
	startpos <- min(dat[,3]) 
	endpos <- max(dat[,3])
	readlen <- endpos - startpos + 1 
	if (startpos == 0) {
		dat[,3] <- dat[,3] + 1
	}
	dat[,3] <- abs(dat[,3]-(readlen + 1)*(dat[,7] == '-'))
	message("\nNumber of bases in training set: ", dim(dat)[1], "\n")
	
	
  # find unique positions and remove positions with coverage < 3 if RefSeq unspecified
	if(RefSeq == ""){
      r <- rle(dat[,2])
	  w <- which(r$lengths<=2)
	  s <- r$values[w]
      m <- merge(dat[,2],cbind(s,s=1),by=c(1),all.x=T)
	  w <- which(m[,2]==1)
	  dat <- dat[-w,]
	  r <- c()
	  w <- c()
	  s <- c()
	  m <- c()
	}  
    unipos <- unique(dat[,1:2])
  
  # if given, read in SNP file, remove SNPs from dat, readpos
    if(SNP != ""){
      message("Removing known SNPs from training set \n")
	  if(length(grep(".Rdata",SNP)) > 0){
	    load(SNP)
	  } else{
        snp <- read.table(SNP)
	  }
      snp <- unique(snp[,1:2])
	  m <- merge(unipos,cbind(snp,snp=1),by=c(1,2),all.x=T)	  
	    w <- which(!is.na(m[,3])==1)
        unipos <- unipos[-w,]
	  m <- merge(dat[,1:2],cbind(snp,snp=1),by=c(1,2),all.x=T)	  
	    w <- which(!is.na(m[,3])==1)
        dat <- dat[-w,]
	  m <- c()
	  w <- c()
	  snp <- c()	# save memory
    }

  # calculate sequencing errors 
    message("Determining Sequencing Errors \n")
    Y <- matrix(0,length(dat[,2]),1) 	
    if(RefSeq == ""){			# RefSeq not specified
	  m <- match(unipos[,2],dat[,2])
	  m <- c(m,dim(dat)[1]+1)
	  useless <- matrix(0,dim(dat)[1],1)
      for(i in 1:nrow(unipos)){ 
        tempseq <- c() 
        tempseq <- dat[m[i]:(m[i+1]-1),4] 
        flagA <- (tempseq=='A')   
        flagC <- (tempseq=='C')   
        flagG <- (tempseq=='G')
        flagT <- (tempseq=='T') 
        useless[m[i]:(m[i+1]-1)] <- (1 - (flagA + flagC + flagG + flagT))
        k <- max(sum(flagA),sum(flagC),sum(flagG),sum(flagT)) 
        temp <- matrix(1,length(tempseq),1)
        if(sum(flagA) == k){ 
          temp[flagA] <- 0
	    } 
        if(sum(flagC) == k){ 
          temp[flagC] <- 0
	    } 
        if(sum(flagG) == k){ 
          temp[flagG] <- 0
	    }
        if(sum(flagT) == k){ 
          temp[flagT] <- 0
	    }
        Y[m[i]:(m[i+1]-1)] <- temp 
      }
    } else if(RefSeq != ""){		# RefSeq file given
	  if(length(grep(".Rdata",RefSeq)) > 0){
	    load(RefSeq)
	  } else{
        ref <- read.table(RefSeq)
	  }
	  ref <- ref[(ref[,1]==dat[1,1]),]
	  useless <- matrix(0,dim(dat)[1],1)	
	  dat2ref <- findInterval(dat[,2],ref[,2])
	  if (max(dat2ref) == max(ref[,2])){		# remove rows with no reference
	    warning("\t Removing bases with no matching location in RefSeq.")
	    m <- merge(dat[,1:2],cbind(ref[,1:2],ref=1),by=c(1,2),all.x=T)		
	    w <- which(is.na(m[,3]) == 1)
  		dat <- dat[-w,]
		warning("\t Removed", length(w), "bases without given RefSeq.\n") 	
	    dat2ref <- findInterval(dat[,2],ref[,2])
		useless <- matrix(0,dim(dat)[1],1)	
		m <- c()
		w <- c()
	  }	
	  flag <- (dat[,4] == 'A')|(dat[,4] == 'C')|(dat[,4] == 'G')|(dat[,4] == 'T')
	  useless <- 1 - flag
	  tempseq <- (ref[dat2ref,3] == 'A') + 2*(ref[dat2ref,3] == 'C') + 3*(ref[dat2ref,3] == 'G') + 4*(ref[dat2ref,3] == 'T')
	  dattemp <- (dat[,4] == 'A') + 2*(dat[,4] == 'C') + 3*(dat[,4] == 'G') + 4*(dat[,4] == 'T')
	  Y <- (tempseq != dattemp)
	  dat2ref <- c()
	  flag <- c()
	  tempseq <- c()
	  dattemp <- c()
	  ref <- c()	# save memory
    } 

  # remove ambiguous positions (> max(nerr,coverage*nraf) errors)
    coverage <- rle(dat[,2])
	thresh <- pmax(nerr, coverage$lengths*nraf)
	w <- which(Y == 1)
	r <- rle(dat[w,2])
    m <- merge(coverage$values,cbind(r$values,matrix(1,length(r$values),1)),by=c(1),all.x=T)	  
       w <- which(!is.na(m[,2])==1)	
	w <- which(r$lengths > thresh[w])
	s <- r$values[w]
    m <- merge(dat[,2],cbind(s,s=1),by=c(1),all.x=T)
	useless <- (useless | !is.na(m[,2]))
	
  # calculate error rate
#    if(sum(useless)>0){
      w <- which(useless == 1)
      Y <- Y[-w]
      dat <- dat[-w,]
#	}
	w <- c()
	useless <- c()
    err <- sum(Y) / length(Y)
    message("Training set error rate = ", sprintf("%.5f",err), "\n") 
	message("Number of bases in filtered training set: ", length(Y), "\n")	

    flag <- which(Y == 1) 
    rp <- dat[flag,3] 
	flag <- c()
    h <- hist(rp,breaks = c(0:readlen),plot=FALSE)
    diagnostics$ReadPosErrors <- h$counts


  # find flagged positions
    FlagPos <- c() 
	pos <- c(1:readlen) 
    for(i in 1:readlen){
      if(diagnostics$ReadPosErrors[i] > 1.5*(sum(Y)/readlen)){ 
        FlagPos <- c(FlagPos,pos[i])
      }
    }


  # Find base called 
    baseI <- cbind(1*(dat[,4]=='A'), 1*(dat[,4]=='C'), 1*(dat[,4]=='G')) 


  # set-up for logistic regression
    flag <- matrix(0,length(dat[,3]),length(FlagPos)) 
    for(i in 1:length(FlagPos)){
      flag[,i] <- (dat[,3] == FlagPos[i]) 
    }
    X <- cbind(dat[,5], 1*(dat[,5] == 0), dat[,6], baseI, dat[,3], flag) 
	baseI <- c()		# save memory
    flag <- c()
	dat <- c()

  
  # fit logistic regression
    message("Calculating Training Parameters \n")
    m <- ceiling(nrow(X)/10000000)
    out.c <- c()
    for(j in 1:m){
      mark <- c(round((nrow(X)/m)*(j-1)+1),min(round((nrow(X)/m)*j),nrow(X))) 
      out <- glm.fit(cbind(matrix(1,mark[2]-mark[1]+1,1),X[mark[1]:mark[2],]),Y[mark[1]:mark[2]],family = binomial())
      out.c <- cbind(out.c,out$coefficient)
    }
    out <- c()
    coeff <- c()
    for(j in 1:nrow(out.c)){
      coeff[j] <- median(out.c[j,])
    }      
	coeff[which(is.na(coeff))] <- 0		# Replace NA's with 0's


  # save coefficients and flagged positions
    diagnostics$FlagPos <- FlagPos
	diagnostics$coeff <- coeff  
    write.table(FlagPos,file="FlagPos.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE)
    write.table(coeff,file="Coeff.txt",sep="\t",quote=FALSE,row.names = FALSE,col.names = FALSE)

  # read in recalibrated data
    inner <- cbind(X=1,X) %*% coeff
    recalqc <- 1/(1+exp(-1*inner))
    recalqc <- floor(-10*log10(recalqc))
	inner <- c()
		
  # calculate quality score frequencies
    maxqc <- max(cbind(recalqc,X[,1]))
    h <- hist(X[,1],breaks=c(-1:maxqc),plot=FALSE)
	counts <- sum(h$counts)
    diagnostics$QualFreqBefore <- h$counts/counts
    h <- hist(recalqc,breaks=c(-1:maxqc),plot=FALSE)
    diagnostics$QualFreqAfter <- h$counts/counts
	
  # calculate empirical error rates
    flag <- which(Y==1)
    h <- hist(X[flag,1],breaks=c(-1:maxqc),plot=FALSE)
    Nbefore_err <- h$counts
    diagnostics$ErrRatesBefore <- -10*log10(Nbefore_err / (diagnostics$QualFreqBefore*counts))
	diagnostics$ErrRatesBefore[diagnostics$ErrRatesBefore > maxqc] <- maxqc
    h <- hist(recalqc[flag],breaks=c(-1:maxqc),plot=FALSE)
    Nafter_err <- h$counts
    diagnostics$ErrRatesAfter <- -10*log10(Nafter_err / (diagnostics$QualFreqAfter*counts))
	diagnostics$ErrRatesAfter[diagnostics$ErrRatesAfter > maxqc] <- maxqc
	h <- c()
	flag <- c()		
	
  # calculate FWSE
    diagnostics$FWSE <- c()
    temp1 <- (diagnostics$ErrRatesBefore - c(0:maxqc))^2
	temp2 <- temp1 * diagnostics$QualFreqBefore
	flag2 <- which(is.na(temp2))
	if (length(flag2)>0){
	  diagnostics$FWSE[1] <- sum(temp2[-flag2])
	} else {
	  diagnostics$FWSE[1] <- sum(temp2)
	}
    temp1 <- (diagnostics$ErrRatesAfter - c(0:maxqc))^2
	temp2 <- temp1 * diagnostics$QualFreqAfter
	flag2 <- which(is.na(temp2))
	if (length(flag2)>0){
	  diagnostics$FWSE[2] <- sum(temp2[-flag2])
	} else {
	  diagnostics$FWSE[2] <- sum(temp2)
	}	
		
  # make recalibration plots
    if(plotname != ""){
      message("Making recalibration plots \n")

      # set up plot
        pdf(plotname)
        par(mfrow=c(2,2))

      # plot (1,1) - errors by read position
        ReadPosErrorPlot(diagnostics$ReadPosErrors, startpos = startpos)
 
      # plot (1,2) - distribution of scores
		QualFreqPlot(diagnostics$QualFreqBefore, diagnostics$QualFreqAfter)

      # plot (2,1) - Before error rates
        f1 <- FWSEplot(diagnostics$ErrRatesBefore, diagnostics$QualFreqBefore, col = "blue", lim = c(0, maxqc),
		   main_title = "Reported vs. Empirical Quality: Before")

      # plot (2,2) - After error rates
        f2 <- FWSEplot(diagnostics$ErrRatesAfter, diagnostics$QualFreqAfter, col = "red", lim = c(0, maxqc),
		   main_title = "Reported vs. Empirical Quality: After")

      dev.off()
	  
    }
	X <- c()		# save memory

  
  # run recalibrator
    message("Recalibrating bam File")
    s2 <- .jarray(c(paste("INPUT=", in_bam, sep=""), paste("BAM_OUTPUT=", out_bam, sep=""), "COEFF=Coeff.txt", "FLAGGED=FlagPos.txt"))       
    bre <- .jnew("org/renci/sequencing/util/BAMRecalibrator")
    .jcall(bre, "V", "main",s2)
    message("Recalibration complete \n")

	
  # delete temporary files
    if(temp_files == 0){
      message("Deleting temporary files \n")
      unlink("FlagPos.txt")
      unlink("Coeff.txt")
    }


  message("\nEnd time: ", date(), "\n") 
  
  diagnostics 

}

