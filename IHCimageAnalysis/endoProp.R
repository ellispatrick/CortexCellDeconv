
# Function for maximum correlation thresholding (MCT)
MCT = function(inputimage){
  i = inputimage*256; 
  bd = 0:256; 
  i = as.vector(i)
  hi = hist(i, breaks = bd,plot = FALSE); 
  bd = bd[-256]
  m = sum(hi$counts*bd)/sum(hi$counts); 
  
  nm = sum(hi$counts); 
  
  nbd = bd - m; 
  
  selfcc =  (t(nbd*nbd)%*%t(t(hi$counts))); 
  
  
  fbe = matrix(0,1,256) 
  for(jj in 2:256){
    beuh = matrix(0,1,256);
    beuh[1,jj:256] = hi$counts[jj:256]; 
    benh = matrix(0:255,1,256); 
    benh = benh - m; 
    
    tvu = matrix(0,1,256); tvu[1,jj:256] = hi$counts[jj:256]; 
    tvl = matrix(0,1,256); tvl[1,1:(jj-1)] = hi$counts[1:(jj-1)];
    njl = sum(tvl); 
    nju = sum(tvu);
    lv = (0-(nju/nm));
    uv = (1-(nju/nm));
    
    fbe[1,jj-1] = (beuh)%*%t(benh)/sqrt(selfcc*(((nm-nju)*(nju))/nm));
    fbe[fbe==Inf|is.na(fbe)]=0;
  }
  
  
  bev = bep = which.max(fbe);
  tv = bep/256; 
  tv
}



library(EBImage)
library(mixtools)
library(parallel)
library(MASS)

## Process file names
tifFiles = list.files('.','-1040.tiff$',recursive = TRUE,full.names = TRUE)
tifFiles = tifFiles[-grep('white',tifFiles,ignore.case = TRUE)]
fileSplit = split(tifFiles,unlist(lapply(strsplit(tifFiles,'\\/'),function(x)x[4])))


set.seed(123)
resultsPECAM = mclapply(fileSplit,function(x){
# Read in images
  nuc = channel(readImage(x[1]),'gray')
  cel = channel(readImage(x[2]),'gray')

# Background correct  
  disc = makeBrush(3, "disc")
  disc = disc / sum(disc)
  nuc3 = nuc/filter2( nuc, disc )
  nuc2 = nuc/nuc3

# Threshold nuclei  
  nmask = nuc>MCT(nuc2/max(nuc2))
  nmask = fillHull(nmask)
  bwN = bwlabel(nmask)

# Background correct endotheliel marker 
  disc = makeBrush(201, "disc")
  disc = disc / sum(disc)
  cel_bg = filter2( cel, disc )
  BG = cel/cel_bg

# Threshold endotheliel cells  
  cut = median(BG)+3*mad(BG)
  cmask = BG>cut
  bwC = bwlabel(cmask)
  
# Calculate size of cells and intensity of PECAM
  sC = split(BG,bwC)[-1]
  LsC = unlist(lapply(sC,length))
  MsC = unlist(lapply(sC,max))
  names(MsC) = names(LsC)
  MsC = resid(rlm(log10(MsC)~log10(LsC)))

# Calculate threshold as a starting point for mixtools  
  mu = c(mean(MsC[LsN>100]),mean(MsC[LsN>100])+2*sd(MsC[LsN>100]))

# Use mixture model to select cells that have strongest PECAM expression
  m =   mixtools::normalmixEM(MsC[LsC>50],mu)  
  cut = (max(min(m$x[m$posterior[,2]>0.5]),min(m$x[m$posterior[,1]>0.5])))
  mask = cmask
  mask[!bwC%in%names(which(MsC[LsC>100]>cut))] = 0
  
# Count how many nuclei are endothelial positive  
  s2 = split(mask,bwN)[-1]
  s2 = s2[unlist(lapply(s2,length)>=200)]
  mean(lapply(s2,sum)>50)
  
},mc.cores = 50,mc.silent = TRUE)


# Create map between images and participants
fileSplit = unlist(lapply(fileSplit,function(x)x[1]))
fileMap = sapply(fileSplit,strsplit,'\\/')
fileMap = unlist(lapply(fileMap,function(x)x[2]))
results = do.call('c',resultsPECAM)
names(results) = fileMap[names(resultsPECAM)]


# Group cell proportion estimates by participant
propPECAM = split(results,names(results))
propPECAM  = lapply(propPECAM ,function(x)as.numeric(x[nchar(x)<30]))

# Average proportion for participants
s = unlist(lapply(propPECAM,mad,na.rm = TRUE))
m = (unlist(lapply(propPECAM,median,na.rm = TRUE)))
propPECAM = unlist(lapply(propPECAM,mean,.5,na.rm = TRUE))

# Filter poor quality images with high variance
propPECAM =propPECAM[propPECAM<0.4&s<0.15]

