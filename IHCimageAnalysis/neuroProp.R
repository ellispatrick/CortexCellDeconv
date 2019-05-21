library(EBImage)
library(parallel)

## Process file names
tifFiles = list.files('.','-1040.tiff$',recursive = TRUE,full.names = TRUE)
tifFiles = tifFiles[-grep('white',tifFiles,ignore.case = TRUE)]
tifFiles = tifFiles[-grep('ImageJ',tifFiles)]

fileSplit = split(tifFiles,unlist(lapply(strsplit(tifFiles,'\\/'),function(x)x[5])))
fileSplit = fileSplit[grep('files',names(fileSplit))]
fileSplit = lapply(fileSplit,function(x)c(x[grep('_dapi',x)][1],x[grep('_NeuN',x)][1]))
fileSplit = fileSplit[unlist(lapply(fileSplit,function(x)sum(!is.na(x))))==2]


set.seed(1234)
resultsNeuN = mclapply(fileSplit,function(x){
# Read in images
  nuc = channel(readImage(x[1]),'gray')
  cel = channel(readImage(x[2]),'gray')

# Threshold nucleus  
  nmask = nuc>median(nuc)+ 3*mad(nuc)
  bwN = bwlabel(nmask)
  
# Perform background subtraction on neuron marker 
  disc = makeBrush(201, "disc")
  disc = disc / sum(disc)
  cel_bg = filter2( cel, disc )
  BG = cel/cel_bg
  
# Threshold neuron marker
  cut = median(BG)+3*mad(BG)
  mask = BG>cut

# Count proportion of neuron positive nuclei
  s2 = split(mask,bwN)[-1]
  s2 = s2[unlist(lapply(s2,length)>=200)]
  mean(lapply(s2,sum)>200)
  
},mc.cores = 50)



# Create map between images and participants
fileSplit = unlist(lapply(fileSplit,function(x)x[1]))
fileMap = sapply(fileSplit,strsplit,'\\/')
fileMap = unlist(lapply(fileMap,function(x)x[3]))
resultsNeuN = resultsNeuN[!unlist(lapply(resultsNeuN,is.null))]

# Group cell proportion estimates by participant
results = do.call('c',resultsNeuN)
names(results) = fileMap[names(resultsNeuN)]
propNeuN = split(results,names(results))


# Average proportion for participants
s = unlist(lapply(propNeuN,mad,na.rm = TRUE))
m = unlist(lapply(propNeuN,median,na.rm = TRUE))
propNeuN = unlist(lapply(propNeuN,median,na.rm = TRUE))

# Filter poor quality images with high variance
propNeuN = propNeuN[m>0.3&s<0.08]




