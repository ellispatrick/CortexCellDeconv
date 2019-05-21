library(EBImage)
library(mixtools)
library(parallel)
library(MASS)

## Process file names
tifFiles = list.files('.','-1040.tiff$',recursive = TRUE,full.names = TRUE)
tifFiles = tifFiles[-grep('white',tifFiles,ignore.case = TRUE)]
fileSplit = split(tifFiles,unlist(lapply(strsplit(tifFiles,'\\/'),function(x)x[4])))



set.seed(1234)
resultsOligo2 = mclapply(fileSplit,function(x){
# Read in images
  nuc = channel(readImage(x[1]),'gray')
  cel = channel(readImage(x[2]),'gray')
  
# Threshold nucleus  
  nmask = nuc>median(nuc)+ 3*mad(nuc)
  bwN = bwlabel(nmask)

# Perform background subtraction on oligo marker 
  disc = makeBrush(201, "disc")
  disc = disc / sum(disc)
  cel_bg = filter2( cel, disc )
  BG = cel/cel_bg

# Threshold oligodendrocyte marker
  cut = median(BG)+3*mad(BG)
  cmask = BG>cut
  bwC = bwlabel(cmask)
  
# Calculate size of cells and intensity of Olig2
  sC = split(BG,bwC)[-1]
  LsC = unlist(lapply(sC,length))
  MsC = log10(unlist(lapply(sC,quantile,0.95)))
  names(MsC) = names(LsC)
  MsC = resid(rlm((MsC)~log10(LsC)))

# Use mixture model to select cells that have strongest Olig2 expression
  m =   mixtools::normalmixEM((MsC[LsC>50&LsC<500]))      
  cut = (max(min(m$x[m$posterior[,2]>0.5]),min(m$x[m$posterior[,1]>0.5])))
  mask = cmask
  mask[!bwC%in%names(which(MsC[LsC>50&LsC<500]>cut))] = 0
  
# Count how many nuclei are oligodendrocyte positive  
  s2 = split(mask,bwN)[-1]
  s2 = s2[unlist(lapply(s2,length)>=200)]
  mean(lapply(s2,sum)>50)
  
},mc.cores = 50,mc.silent = TRUE)


# Create map between images and participants
fileSplit = unlist(lapply(fileSplit,function(x)x[1]))
fileMap = sapply(fileSplit,strsplit,'\\/')
fileMap = unlist(lapply(fileMap,function(x)x[2]))
results = do.call('c',resultsOligo2)
names(results) = fileMap[names(resultsOligo2)]

# Group cell proportion estimates by participant
propOligo2 = split(results,names(results))
propOligo2  = lapply(propOligo2 ,function(x)as.numeric(x[nchar(x)<30]))

# Average proportion for participants
s = unlist(lapply(propOligo2,mad,na.rm = TRUE))
m = unlist(lapply(propOligo2,median,na.rm = TRUE))
propOligo2 = unlist(lapply(propOligo2,mean,0.5,na.rm = TRUE))

# Filter poor quality images with high variance
propOligo2 = propOligo2[s<0.15&m<.22]

