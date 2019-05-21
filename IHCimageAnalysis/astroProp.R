library(EBImage)
library(parallel)

## Process file names
tifFiles = list.files('.','-1040.tiff$',recursive = TRUE,full.names = TRUE)
tifFiles = tifFiles[-grep('white',tifFiles,ignore.case = TRUE)]
fileSplit = split(tifFiles,unlist(lapply(strsplit(tifFiles,'\\/'),function(x)x[5])))

set.seed(123)
resultsGFAP = mclapply(fileSplit,function(x){
# Read in images
  nuc = channel(readImage(x[1]),'gray')
  cel = channel(readImage(x[2]),'gray')
# Threshold nucleus  
  nuc = (nuc-min(nuc))/max(nuc-min(nuc))
  cut = median(nuc) + 3*mad(nuc)
  nmask = nuc> cut
  bwN = bwlabel(nmask)
  
  
# Perform background subtraction on astrocyte marker 
  BG = cel
  disc = makeBrush(201, "disc")
  disc = disc / sum(disc)
  cel_bg1 = filter2( BG, disc )
  BG = BG/cel_bg1
  BG = BG-min(BG)
  BG = BG/max(BG)
# Perform further blurring to impove overlap with nucleus  
  BG = gblur(BG,8)
  BG = BG-min(BG)
  BG = BG/max(BG)
  
# Threshold  astrocytes, create mask, ignore large objects
  cmask = BG > median(BG) + 3*mad(BG)
  bwC = bwlabel(cmask)
  tab = table(bwC)
  cmask[bwC%in%names(which(tab>3000))] = 0
  
# Count how many nuclei are astrocyte positive  
  s2 = split(cmask,bwN)[-1]
  s2 = s2[unlist(lapply(s2,length))>=200]
  mean(unlist(lapply(s2,sum))>200)
  
  
},mc.cores = 50,mc.silent = TRUE)


# Create map between images and participants
fileSplit = unlist(lapply(fileSplit,function(x)x[1]))
fileMap = sapply(fileSplit,strsplit,'\\/')
fileMap = unlist(lapply(fileMap,function(x)x[3]))

# Group cell proportion estimates by participant
resultsGFAP = resultsGFAP[!unlist(lapply(resultsGFAP,is.null))]
results = do.call('c',resultsGFAP)
names(results) = fileMap[names(resultsGFAP)]

# Average proportion of participants
propGFAP = split(results,names(results))
s = unlist(lapply(propGFAP,mad,na.rm = TRUE))
m = unlist(lapply(propGFAP,mean,na.rm = TRUE))
propGFAP = unlist(lapply(propGFAP,mean,.5,na.rm = TRUE))

# Filter poor quality images with high variance
propGFAP =propGFAP[s<0.08]
