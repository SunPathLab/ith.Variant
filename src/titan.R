#!/srv/gsfs0/projects/curtis/ruping/tools/R/bin/Rscript

## this is for WGS or WXS titan run

inputpar <- commandArgs(TRUE)
if (length(inputpar) < 16) stop("Wrong number of input parameters: 'path sampleName alleleCount tumorWig normalWig gcWig mapWig plp plpe normalc normalcm symmetric transtate tranclone exons(if WXS)'")


path <- inputpar[1]
src <- inputpar[2]
sampleName <- inputpar[3]
gender <- inputpar[4]
alleleCount <- inputpar[5]
tumorWig <- inputpar[6]
normalWig <- inputpar[7]
gcWig <- inputpar[8]
mapWig <- inputpar[9]
plp <- inputpar[10]
plpe <- inputpar[11]
normalc <- inputpar[12]
normalcm <- inputpar[13]
symmetric <- inputpar[14]
transtate <- inputpar[15]
tranclone <- inputpar[16]
exons <- inputpar[17]

library(TitanCNA)
library(HMMcopy)
library(doMC)
library(S4Vectors)

setwd(path)

#run titan
runTitan <- function(sampleName, gender, snpFile, tumWig, normWig, gc, map, plp, plpe, normalc, normalcm, symmetric, transtate, tranclone, exons="SRP") {

    #prepare data
    snpData <- loadAlleleCounts(snpFile, symmetric=symmetric)
    
    if (exons != "SRP") {
      cnData <- correctReadDepth(tumWig, normWig, gc, map, targetedSequence = exons)
    } else {
      cnData <- correctReadDepth(tumWig, normWig, gc, map)
    }
    logR <- getPositionOverlap(snpData$chr, snpData$posn, cnData)
    snpData$logR <- log(2^logR) #transform the log ratio to natural logs
    snpData <- filterData(snpData, c(1:22,"X"), minDepth = 10, maxDepth = 500, positionList = NULL)
    #prepare data

    registerDoMC(cores = 2)
    titancnaresults <- vector("list",2)
    
    for (j in 1:2) {
        numClusters <- j
        params <- loadDefaultParameters(copyNumber = 8, numberClonalClusters = numClusters)
        K <- length(params$genotypeParams$alphaKHyper)
        params$genotypeParams$alphaKHyper <- rep(1000, K)
        params$normalParams$n_0 <- normalc
        params$ploidyParams$phi_0 <- plp

        convergeParams <- runEMclonalCN(snpData, params = params,
                                        maxiter = 20, maxiterUpdate = 1500,
                                        useOutlierState = FALSE, txnExpLen = transtate,
                                        txnZstrength = tranclone,
                                        normalEstimateMethod = normalcm,
                                        estimateS = TRUE, estimatePloidy = plpe)
        
        optimalPath <- viterbiClonalCN(snpData, convergeParams)
        if (length(unique(optimalPath)) == 1) next
        results <- outputTitanResults(snpData, convergeParams, optimalPath,
                                      filename = NULL, posteriorProbs = FALSE,
                                      subcloneProfiles = TRUE)

        convergeParams <- results$convergeParams   #reset convergeParams
        results <- results$corrResults             #corrected result table
        norm <- tail(convergeParams$n, 1)
        ploidy <- tail(convergeParams$phi, 1)
        results$AllelicRatio = 1-as.numeric(results$AllelicRatio)    # reverse the allelic ratio
        cellularity = 1 - convergeParams$s[, ncol(convergeParams$s)] # estimated cellular prevalence
        
        titancnaresults[[j]] <- list(S_DbwIndex=computeSDbwIndex(results)$S_DbwIndex,results=results,
                                     convergeParams=convergeParams)

        #generate segmentation files, original/uncorrect one
        segmenttmp = titancna2seg(results, convergeParams)
        write.table(segmenttmp, file=paste(sampleName,"_nclones",numClusters,".TitanCNA.segments.txt",sep=""),
                    quote = F, row.names = F, sep = "\t")


        # under construction for correcting intCN for chrX
        outigv = paste(sampleName,"_nclones",numClusters,".TitanCNA.segments.igv.txt",sep="")
        segs <- outputTitanSegments(results, sampleName, convergeParams, filename = NULL, igvfilename = outigv)
        corrIntCN.results <- correctIntegerCN(results, segs, 1 - norm, ploidy, maxCNtoCorrect.autosomes = 8,
                                                                                            maxCNtoCorrect.X = NULL, minPurityToCorrect = 0.2, gender = gender, chrs = c(1:22, 'X'))
        results.cn <- corrIntCN.results$cn
        segs <- corrIntCN.results$segs
        outfile = paste(sampleName,"_nclones",numClusters,".TitanCNA.corrIntCN.txt",sep="")
        outseg = paste(sampleName,"_nclones",numClusters,".TitanCNA.corrIntCNseg.txt",sep="")
        outparam = paste(sampleName,"_nclones",numClusters,".TitanCNA.params.txt",sep="")
        message("Writing titan seg results to ", outfile, ", ", outseg, ", ", outparam)
        write.table(results.cn, file = outfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
        write.table(segs, file = outseg, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
        outputModelParameters(convergeParams, results.cn, outparam)
        # under construction
        
        #if (j == 1) {
        #    rawTable = titancna2seg(results, convergeParams, raw=TRUE)
        #    write.table(rawTable, file=paste(sampleName,".TitanCNA.rawTable.txt",sep=""),
        #                quote = F, row.names = F, sep = "\t")
        #}
        
        #make plots
        if (exons == "SRP") {  #WGS
            
            allstate <- paste(results$Chr,results$TITANstate,results$ClonalCluster)
            changepoints <- c(1,which(allstate[-1] != allstate[-length(allstate)])+1)
            segments <- results[changepoints,c("Chr",rep("Position",2),"AllelicRatio","LogRatio")]
            names(segments)[2:3] <- c("Position1","Position2")
            segments$Position2 <- results$Position[c(changepoints[-1]-1,length(allstate))]
            segments[[2]] <- as.numeric(segments[[2]])
            segments[[3]] <- as.numeric(segments[[3]])
            segments[[4]] <- as.numeric(segments[[4]])
            segments[[5]] <- as.numeric(segments[[5]])
            segments$NumMarker <- diff(c(changepoints,length(allstate)+1))
            for (k in 1:nrow(segments)) {
                af <- as.numeric(results$AllelicRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
                segments$AllelicRatio[k] <- mean(pmax(af,1-af))
                lr <- as.numeric(results$LogRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
                segments$LogRatio[k] <- mean(lr)
            }
            xlim1 <- quantile(rep(segments$AllelicRatio,segments$NumMarker),c(0.0001,0.9999))
            ylim1 <- quantile(rep(segments$LogRatio,segments$NumMarker),c(0.0001,0.9999))

            #plotting for each chromosome
            for (chro in c(1:22,"X")) {
                pdf(paste(sampleName,"_nclones",numClusters,"_chr", chro, ".TitanCNA.pdf",sep=""),
                    width=12, height=6)
                if (is.null(titancnaresults[[j]])) next
                SD <- round(titancnaresults[[j]]$S_DbwIndex,3)
                nclones <- nrow(convergeParams$s)
                ploidy <- round(tail(convergeParams$phi, 1),2)
                ploidy2 <- ploidy * (1 - norm) + 2 * norm
                meandepth <- round(mean(as.numeric(results$Depth)),2)
                npoints <- nrow(results)
                s <- round(convergeParams$s[1,ncol(convergeParams$s)],2)

                
                layout(matrix(c(1,2,3,3),nrow=2),widths=c(2,1))
                par(pty="m")
                par(mar=c(4,4,2,1))
                plotCNlogRByChr(results.cn, chr = chro, normal=norm, ploidy = ploidy, ylim = c(-2, 2), cex=0.25,
                                main=paste(sampleName, " nc=", numClusters, sep=""),
                                xlab=paste("normC=", round(norm,3), " pl=", ploidy,
                                           " cellularity=", round(cellularity,3),
                                           " SD=",SD," s=",s," nc=",nclones," np=",npoints,
                                           " md=",meandepth,sep=""),
                                cex.lab=0.8)
                par(mar=c(4,4,2,1))
                plotAllelicRatio(results.cn, chr = chro, ylim = c(0, 1),
                                 cex = 0.25, xlab = paste("Chromosomes", chro, sep=" "),
                                 main = "", cex.lab=0.8)

                #plot bubble like
                par(mar=c(7,4,6,1))
                smoothScatter(rep(segments$AllelicRatio,segments$NumMarker),
                              colramp=colorRampPalette(terrain.colors(32)),
                              rep(segments$LogRatio,segments$NumMarker),
                              xlim=xlim1,ylim=ylim1,
                              main = sampleName, xlab="Allelic ratio",ylab="Log ratio")
                abline(h=log2(2/ploidy2), lty=3)
                dev.off()
            }
            
            #plotting for all chromosome
            pdf(paste(sampleName,"_nclones",numClusters,".TitanCNA.pdf",sep=""),width=12, height=6)
            layout(matrix(c(1,2,3,3),nrow=2),widths=c(2,1))
            par(pty="m")
            par(mar=c(4,4,2,1))
            plotCNlogRByChr(results.cn, normal=norm, ploidy = ploidy, ylim = c(-2, 2), cex=0.25, chr=c(1:22,"X"),
                            main=paste(sampleName, " nc=", numClusters, sep=""),
                            xlab=paste("normC=", round(norm,3), " pl=", ploidy, " cellularity=",
                                       round(cellularity,3), " SD=",SD," s=",s," nc=",nclones,
                                       " np=",npoints," md=",meandepth,sep=""), cex.lab=0.8)
            par(mar=c(4,4,2,1))
            plotAllelicRatio(results.cn, ylim = c(0, 1), chr=c(1:22,"X"),
                             cex = 0.25, xlab = "Chromosomes", main = "", cex.lab=0.8)

            #plot bubble like
            par(mar=c(7,4,6,1))
            smoothScatter(rep(segments$AllelicRatio,segments$NumMarker),
                          colramp=colorRampPalette(terrain.colors(32)),
                          rep(segments$LogRatio,segments$NumMarker),
                          xlim=xlim1,ylim=ylim1,
                          main = sampleName, xlab="Allelic ratio",ylab="Log ratio")
            abline(h=log2(2/ploidy2), lty=3)
            dev.off()
            
      } else if (exons != "SRP") { #WES
          
          if (is.null(titancnaresults[[j]])) next
          SD <- round(titancnaresults[[j]]$S_DbwIndex,3)
          nclones <- nrow(convergeParams$s)
          ploidy <- round(tail(convergeParams$phi, 1),2)
          ploidy2 <- ploidy * (1 - norm) + 2 * norm
          meandepth <- round(mean(as.numeric(results$Depth)),2)
          npoints <- nrow(results)
          s <- round(convergeParams$s[1,ncol(convergeParams$s)],2)

          pdf(paste(sampleName,"_nclones",numClusters,".TitanCNA.pdf",sep=""),width=12, height=6)
          layout(matrix(c(1,2,3,3),nrow=2),widths=c(2,1))
          par(pty="m")
          par(mar=c(4,4,2,1))
          plotCNlogRByChr(results.cn, ploidy = ploidy, normal=norm, ylim = c(-2, 2), cex=0.25, chr=c(1:22,"X"),
                          main=paste(sampleName, " nc=", numClusters, sep=""),
                          xlab=paste("normC=", round(norm,3), " pl=", ploidy,
                                     " cellularity=", round(cellularity,3),
                                     " SD=",SD," s=",s," nc=",nclones," np=",npoints,"
                                     md=",meandepth,sep=""),
                          cex.lab=0.8)
          plotAllelicRatio(results.cn, ylim = c(0, 1), chr=c(1:22,"X"),
                           cex = 0.25, xlab = "Chromosomes", main = "", cex.lab=0.8)

          #plot bubble like
          allstate <- paste(results$Chr,results$TITANstate,results$ClonalCluster)
          changepoints <- c(1,which(allstate[-1] != allstate[-length(allstate)])+1)
          segments <- results[changepoints,c("Chr",rep("Position",2),"AllelicRatio","LogRatio")]
          names(segments)[2:3] <- c("Position1","Position2")
          segments$Position2 <- results$Position[c(changepoints[-1]-1,length(allstate))]
          segments[[2]] <- as.numeric(segments[[2]])
          segments[[3]] <- as.numeric(segments[[3]])
          segments[[4]] <- as.numeric(segments[[4]])
          segments[[5]] <- as.numeric(segments[[5]])
          segments$NumMarker <- diff(c(changepoints,length(allstate)+1))
          for (k in 1:nrow(segments)) {
              af <- as.numeric(results$AllelicRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
              segments$AllelicRatio[k] <- mean(pmax(af,1-af))
              lr <- as.numeric(results$LogRatio[changepoints[k]:(changepoints[k]+segments$NumMarker[k]-1)])
              segments$LogRatio[k] <- mean(lr)
          }
          xlim1 <- quantile(rep(segments$AllelicRatio,segments$NumMarker),c(0.0001,0.9999))
          ylim1 <- quantile(rep(segments$LogRatio,segments$NumMarker),c(0.0001,0.9999))
          par(mar=c(7,4,6,1))
          smoothScatter(rep(segments$AllelicRatio,segments$NumMarker),colramp=colorRampPalette(terrain.colors(32)),
                rep(segments$LogRatio,segments$NumMarker),
                xlim=xlim1,ylim=ylim1,
                main = sampleName, xlab="Allelic ratio",ylab="Log ratio")
          abline(h=log2(2/ploidy2), lty=3)
          dev.off()
      }
        
    } #how many clones
    save(titancnaresults,file=paste(sampleName,".TitanCNA.RData",sep=""))
}



titancna2seg <- function(titanresult,titanparams,raw=FALSE) {

  major_cn_code <- c(0,1,2,1,3,2,4,3,2,5,4,3,6,5,4,3,7,6,5,4,8,7,6,5,4)
  minor_cn_code <- c(0,0,0,1,0,1,0,1,2,0,1,2,0,1,2,3,0,1,2,3,0,1,2,3,4)

  ploidy <- round(tail(titanparams$phi,1),3)
  n <- round(tail(titanparams$n,1),3)

  titanresult$Position <- as.integer(titanresult$Position)
  titanresult$LogRatio <- as.numeric(titanresult$LogRatio)
  titanresult$AllelicRatio <- as.numeric(titanresult$AllelicRatio)
  titanresult$CopyNumber <- as.numeric(titanresult$CopyNumber)
  titanresult$CellularPrevalence <- as.numeric(titanresult$CellularPrevalence)
  titanresult$ClonalCluster[is.na(titanresult$ClonalCluster)] <- 0

  cp2 <- c(which(titanresult$TITANstate[-1] != titanresult$TITANstate[-nrow(titanresult)] | 
                   titanresult$Chr[-1] != titanresult$Chr[-nrow(titanresult)] | 
                     titanresult$ClonalCluster[-1] !=  titanresult$ClonalCluster[-nrow(titanresult)]),
           nrow(titanresult))
  cp1 <- c(1,cp2[-length(cp2)]+1)

  cnv <- data.frame(chrom=titanresult$Chr[cp1],
                    loc.start=titanresult$Position[cp1],
                    loc.end=titanresult$Position[cp2],
                    num.mark=cp2-cp1+1,
                    seg.mean=titanresult$LogRatio[cp1],
                    copynumber=titanresult$CopyNumber[cp1],
                    minor_cn=minor_cn_code[as.integer(titanresult$TITANstate[cp1])+1],
                    major_cn=major_cn_code[as.integer(titanresult$TITANstate[cp1])+1],
                    allelicratio=titanresult$AllelicRatio[cp1],
                    LOHcall=titanresult$TITANcall[cp1],
                    cellularprevalence=titanresult$CellularPrevalence[cp1],
                    ploidy=ploidy,
                    normalproportion=n)
  for (j in 1:length(cp1)) {
    cnv$seg.mean[j] <- mean(titanresult$LogRatio[cp1[j]:cp2[j]])
    cnv$allelicratio[j] <- mean(0.5+abs(0.5-titanresult$AllelicRatio[cp1[j]:cp2[j]]))
    if (j < length(cp1)) {
      if (titanresult$Chr[cp2[j]] == titanresult$Chr[cp1[j+1]]) {
        cnv$loc.end[j] <- round((cnv$loc.end[j]+cnv$loc.start[j+1])/2)
      }
    }
    if (j > 1) {
      if (titanresult$Chr[cp1[j]] == titanresult$Chr[cp2[j-1]]) {
        cnv$loc.start[j] <- cnv$loc.end[j-1]+1
      }
    }
  }
 
  cnv$logcopynumberratio <- log2(((cnv$copynumber - 2) * cnv$cellularprevalence + 2)/2)
  cnv$logcopynumberratio[is.na(cnv$logcopynumberratio)] <- 0

  for (i in c(5,12)) {
    cnv[[i]] <- round(cnv[[i]],3)
  }
  if (raw == FALSE){
      return(cnv)
  } else {
      rawRes = data.frame(Chr=titanresult$Chr,
                          Position=titanresult$Position, LogRatio=titanresult$LogRatio,
          AllelicRatio=titanresult$AllelicRatio)
      return(rawRes)
  }

}




#process input par
if (plpe == "FALSE") {
    plpe = FALSE
} else {
    plpe = TRUE
}

if (symmetric == "FALSE") {
    symmetric = FALSE
} else {
    symmetric = TRUE
}


plp = as.numeric(plp)
message(plp)
message(plpe)
normalc = as.numeric(normalc)
message(normalc)
message(normalcm)
message(symmetric)
message(sampleName)
message(alleleCount)
message(tumorWig)
message(normalWig)
message(gcWig)
message(mapWig)
transtate = as.numeric(transtate)
message(transtate)
tranclone = as.numeric(tranclone)
message(tranclone)

if (exons != "SRP") {   #WES
    targetRegion = read.delim(exons, header=F, as.is=T)
    targetRegion = data.frame(targetRegion[,1:3])
    colnames(targetRegion) = c("chr", "start", "end")
    message(paste(colnames(targetRegion), collapse="\t"))
    message(paste(dim(targetRegion), collapse="\t"))
    targetRegion[,1] = gsub("chr","",targetRegion[,1])     # replace the chr prefix in the target bed file
    runTitan(sampleName,gender,alleleCount,tumorWig,normalWig,gcWig,mapWig,plp,plpe,normalc,normalcm,symmetric,transtate,tranclone,targetRegion)
} else if (exons == "SRP") {  #WGS 
    runTitan(sampleName,gender,alleleCount,tumorWig,normalWig,gcWig,mapWig,plp,plpe,normalc,normalcm,symmetric,transtate,tranclone)
}
