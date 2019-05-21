library(EBImage)
library(parallel)
library(mixtools)

## Process file names
tifFiles = list.files('.','-1040.tiff$',recursive = TRUE,full.names = TRUE)
tifFiles = tifFiles[-grep('white',tifFiles,ignore.case = TRUE)]
fileSplit = split(tifFiles,unlist(lapply(strsplit(tifFiles,'\\/'),function(x)x[5])))


set.seed(1)
resultsIBA1 = mclapply(fileSplit,function(x){
  
# Read in images
  nuc = channel(readImage(x[1]),'gray')
  cel = channel(readImage(x[2]),'gray')

# Threshold nucleus  
  nmask = nuc>median(nuc)+ 3*mad(nuc)
  bwN = bwlabel(nmask)
  
# Perform background subtraction on microglia marker 
  disc = makeBrush(201, "disc")
  disc = disc / sum(disc)
  cel_bg = filter2( cel, disc )
  BG = cel/cel_bg

# Threshold and segment microglia  
  cut = median(BG)+3*mad(BG)
  cmask = BG>cut
  bwC = bwlabel(cmask)
  
# Calculate size of cells and intensity of IBA1
  sC = split(BG,bwC)[-1]
  LsC = unlist(lapply(sC,length))
  MsC = unlist(lapply(sC,quantile,0.95))
  names(MsC) = names(LsC)
  
# Use mixture model to select cells that have strongest IBA1 expression
  m =   mixtools::normalmixEM(log(MsC[LsC>20]))      
  cut = exp(max(min(m$x[m$posterior[,2]>0.5]),min(m$x[m$posterior[,1]>0.5])))
  mask = cmask
  mask[!bwC%in%names(which(MsC[LsC>100]>cut))] = 0
  
# Count how many nuclei are microglia positive  
  s2 = split(mask,bwN)[-1]
  s2 = s2[unlist(lapply(s2,length)>=200)]
  mean(lapply(s2,sum)>200)
  
},mc.cores = 50,mc.silent = TRUE)


# Create map between images and participants
fileSplit = unlist(lapply(fileSplitIBA1,function(x)x[1]))
fileMap = sapply(fileSplit,strsplit,'\\/')
fileMap = unlist(lapply(fileMap,function(x)x[3]))

# Group cell proportion estimates by participant
results = do.call('c',resultsIBA1)
names(results) = fileMap[names(results)]
propIBA1 = split(results,names(results))
propIBA1  = lapply(propIBA1 ,function(x)as.numeric(x[nchar(x)<30]))


# Average proportion of participants
propIBA1 = unlist(lapply(propIBA1,median,na.rm = TRUE))


