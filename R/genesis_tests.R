#== Args 
args<-commandArgs(TRUE)
#===mandatory parameters
snpinfo.file <- args[1]
genotype.files <- args[2]
output.file <- args[3]
#==optional parameters


# added these to JSON
BUFFER <- as.numeric(args[4]) 
gene.file <- args[5] 
agg.file <- args[6] 
snp.filter <- args[7] 
gene.filter <- args[8]
top.maf <- as.numeric(args[9]) 
test.stat <-  args[10] # Score, Wald, Firth
test.type  <-  args[11] # Burden, Single, SKAT

min.mac <- as.integer(args[12])
weights <- args[13]
#conditional <- args[20]
user_cores <-  as.numeric(args[14])
nullmod <-  args[15]

# loads null model and annotated data frame
load(nullmod)
ls()
# get samples included in null model
sample_ids <- nullmod$scanID


USE_AGG = F
USE_SNPINFO = F

# GLOBAL VARIABLES
collapsing.tests <- c("SKAT",  "Burden")
test.type.vals <- c("Single","SKAT", "Burden")
test.stat.vals <- c("Score", "Wald", "Firth")



getMAC <- function(gds.file){
  geno.dat <- seqGetData(gds.file,"genotype")
  apply(geno.dat,3,function(x){min(sum(x == 0,na.rm=T),sum(x == 1,na.rm=T))})
}

cat('output.file',output.file,'\n')
cat('buffer',BUFFER,'\n')
cat('genefile',gene.file,'\n')
cat('varaggfile',agg.file,'\n')
cat('snp.filter',snp.filter,'\n')
cat('top.maf',top.maf,'\n')
cat('test.stat',test.stat,'\n')
cat('test.type',test.type,'\n')

if (!(test.stat %in% test.stat.vals)){
     msg = paste("The requested test statistic:", test.stat, "is not available (Use Firth, Score, Wald!")
     stop(msg)
}


if(!test.type %in% test.type.vals){
    stop("Test type must be one of ",paste(test.type.vals,sep=','))
}


weights = eval(parse(text=weights))
cat('Weights',weights,class(weights),'\n')

suppressMessages(library(SeqArray))
suppressMessages(library(SeqVarTools))
suppressMessages(library(GWASTools))
suppressMessages(library(gap))
suppressMessages(library(Matrix))
suppressMessages(library(plyr))
suppressMessages(library(gdsfmt))
suppressMessages(library(bdsmatrix))
suppressMessages(library(parallel))
suppressMessages(library(GENESIS))
suppressMessages(library(data.table))
suppressMessages(library(doMC))

num_cores <- detectCores(logical=TRUE)


registerDoMC(cores=min(c(user_cores,num_cores)))
cat('Running Analysis with ',min(c(as.numeric(user_cores),num_cores)),' cores of ',num_cores,'\n')

## Setup
source("/home/dnanexus/pipelineFunctions.R")


## snp info
if(snpinfo.file != 'NO_SNPINFO_FILE'){
    USE_SNPINFO = T
    cat('Reading snpinfo....')
    snpinfo <- fread(snpinfo.file,data.table=F)
    if(!all(c("CHR","POS") %in% names(snpinfo))){
    
        msg = paste("SNPinfo file must contain column names 'CHR' and 'POS'")
        stop(msg)
    }
    cat('done\n')

    cat('Input SNPINFO N=',nrow(snpinfo),'\n')
    snpinfo = eval(parse(text= paste0('subset(snpinfo,',snp.filter,')')))
    cat('Output SNPINFO N=',nrow(snpinfo),'\n')
    cat('SNPinfo format\n')
    head(snpinfo,2)
    snpinfo = as.data.frame(snpinfo)[,c('CHR','POS')]
}



# For GDS files
f <- seqOpen(genotype.files)
pos = seqGetData(f, "position")
variant.ids = seqGetData(f, "variant.id")





if(agg.file != 'NO_VAR_AGG_FILE'){
    USE_AGG = T
    system.time({ agg = fread(agg.file, stringsAsFactors=F, sep=',', header=T,data.table=F) })
    if(any(grepl('chromosome',names(agg)))){
        agg$chr = sub('chr','',agg$chromosome)
    }else if(any(grepl('chr',names(agg)))){
        agg$chr = sub('chr','',agg$chr)
    }
    agg$var_long.id = paste(agg$chr,agg$pos,agg$ref,agg$alt,sep=':')

    
    pos = seqGetData(f, "position")
    
    variant.ids = seqGetData(f, "variant.id")
    chr = seqGetData(f, "chromosome")
    allele.ref = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][1])
    allele.alt = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][2])
    var_long.ids = paste(chr,pos,allele.ref,allele.alt,sep=':')
    genes <- unique(agg$group_id)
    cat('N AggUnits=',length(genes),'done\n')
}else{


    ##  Gene list
    ##
    ## Aggregation file
    ##  must contain 'start','stop' and 'CHR'
    ##  start and stop must be numeric
    ##
    cat('reading gene file...')
    if(gene.file == "NO_GENE_REGION_FILE" & test.type != 'Single'){
        stop('For aggregate tests you must provide a aggregation file.')
  
    }else if(gene.file == "NO_GENE_REGION_FILE" & test.type == 'Single'){
  
        batchsize = 500000
        if(gene.file == "NO_GENE_REGION_FILE"){
            nbatch = ceiling(max(pos)/batchsize)
            kg = data.frame('name'=paste0('batch',1:nbatch),'start'=batchsize*((1:nbatch)-1),'stop'=batchsize*1:nbatch,stringsAsFactors=F)
        }
  
        genes <- kg$name
    
#}else if(variantLists == 1){
#  system.time({ kg = fread(gene.file, stringsAsFactors=F, sep=',', header=T) })
#  kg
    
    }else{
        system.time({ kg = fread(gene.file, stringsAsFactors=F, sep=',', header=T) })
        #if(names(kg) %in% c('CHR','CHROM','Chr','Chrome','chrom','end')){
            names(kg)[names(kg) == 'CHR'] = 'chr'
            names(kg)[names(kg) == 'CHROM'] = 'chr'
            names(kg)[names(kg) == 'Chr'] = 'chr'
            names(kg)[names(kg) == 'chrom'] = 'chr'
            names(kg)[names(kg) == 'end'] = 'stop'
            #warning("Will enaming gene file columns from 'CHR','CHROM','Chr','Chrome','chrom','end' to match 'chr','start','stop'\n")
        #}
        if(!(sum(names(kg) %in% c('start','stop','chr')) == 3 & is.numeric(kg$start) & is.numeric(kg$stop))){
            stop("Aggreagation file must contain columns 'start','stop' and 'chr'.  Columns start and stop must be numeric")
        }

        if(!( is.numeric(kg$start) & is.numeric(kg$stop))){
            stop("Columns start and stop in aggregation file must be numeric")
        }

        if( ! all(grepl('chr',kg$chr))){
            stop("chr column in aggregation file must be formated as 'chr#' (e.g. chr1) ")
        }

        cat('Aggregate file filter',gene.filter,'\n')
        kg = eval(parse(text= paste0('subset(kg,',gene.filter,')')))
        genes <- kg$name
    }
    cat('NGENEs=',NROW(kg),'done\n')
}


seqClose(f)



mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
sm_obj <- 
foreach (current.gene=genes, 
         .combine=function(...){rbindlist(list(...),fill=T)},
         .multicombine = TRUE,
         .inorder=FALSE,  
         .options.multicore=mcoptions) %dopar% {

##############
             ## Apply variant filters to get a list of variants in each gene
###############
             
             if(USE_AGG){
                 geneSNPinfo = agg[agg$group_id == current.gene,]
                 snp_idx = variant.ids[var_long.ids %in% geneSNPinfo$var_long.id]
                 cat(current.gene,"\tnSNPS: ",length(snp_idx),"\n")
             }else{
                 gidx <- which(kg$name == current.gene)
                 if (length(gidx) != 1) {
                     msg=paste("The assumptions that the list of genes is unique is violated\n",
                               "The gene:", cgene, "repeats", length(gidx), "times")
                     stop(msg)
                 }
                 if(USE_SNPINFO){ 

                     geneSNPinfo = subset(snpinfo, (POS > (kg[gidx,]$start - BUFFER) & POS <= (kg[gidx,]$stop + BUFFER)))
                     snp_idx = variant.ids[which(pos %in%  geneSNPinfo$POS)]
                     cat(current.gene,'\t',' NSNPS:',length(snp_idx),': range(',(kg[gidx,]$start - BUFFER),'-',(kg[gidx,]$stop + BUFFER),')\n')
                 }else{
                     snp_idx = variant.ids[pos > (kg[gidx,]$start - BUFFER) & pos <= (kg[gidx,]$stop + BUFFER)]
                     cat(current.gene,"\tNSNPS: ",length(snp_idx),"\n")
                 }
             }
  

    
    ## extract genotypes
    f <- seqOpen(genotype.files)
    
    seqSetFilter(f,variant.id = snp_idx,sample.id = sample_ids,verbose=TRUE)
    if(length(seqGetData(f,'variant.id')) > 0){
        ## filter to maf and remove monomorphic

        
        if (test.type %in% collapsing.tests){
            #maf <- seqAlleleFreq(f)
            #maf <- ifelse(maf < 0.5, maf, 1-maf)

            ref.freq <- seqAlleleFreq(f)
            rmaf <- pmin(ref.freq, 1-ref.freq)	

            maf.filt <- rmaf<= top.maf & rmaf>0
            message(paste("Running on", sum(maf.filt), "variants with MAF <=", top.maf))
            seqSetFilter(f, variant.sel=maf.filt, action="intersect", verbose=TRUE)
           
        }else{

            filterByMAF(f, sample_ids,  max(min.mac,1), NA)
        }
      
        num.polymorphic.snps <- seqSummary(f, "genotype", check="none", verbose=FALSE)[["seldim"]][3]
    
        if(num.polymorphic.snps > 0){
      
            genotype.data <- SeqVarData(f, sampleData=annotated.frame)
            ## Collapse test
            if (test.type %in% collapsing.tests) {
                xlist = list()
                xlist[[1]] = data.frame('variant.id'=seqGetData(f,"variant.id"),
                                        'allele.index'=rep(1,length(seqGetData(f,"variant.id")))
                                        )
                system.time({
                    collapse.results <- assocTestSeq( genotype.data, 
                                         nullmod, 
                                         xlist, 
                                         weight.beta = weights,
                                         test=test.type, 
                                         burden.test=test.stat,verbose=F)
                })
                generes <- cbind(data.frame(gene=current.gene, num.variants= nrow(xlist[[1]])),cmac=NA, collapse.results$results)
                #generes <- cbind(data.frame(gene=current.gene, num.variants= nrow(xlist[[1]])),cmac=sum(mac[filtered.alleles & rmaf > 0]), collapse.results$results)
                if(!USE_AGG){
                    generes$start=kg[gidx,]$start - BUFFER
                    generes$stop=kg[gidx,]$stop + BUFFER
                    generes$chr=kg[gidx,]$chr
                }
                #generes$cmaf=sum(rmaf[filtered.alleles & rmaf > 0])
            
            } else {  # Single variant test
 
                system.time({generes <- assocTestMM(genotype.data, 
                                              nullmod, snp.block.size = 2000,
                                              test = test.stat,verbose=FALSE)})
                generes$pos=seqGetData(f, "position")
                generes$snpID = paste0(generes$chr,':',generes$pos)
                generes$ref = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][1])
                generes$alt = sapply(seqGetData(f, "allele"),FUN=function(x)strsplit(x,',')[[1]][2])      
          
            }
        }
    }
    seqClose(f)
             
  if(!exists("generes")){
    generes= data.frame();
  }
  
  generes
}
cat("\nFinished Loop\n")
cat("Total output lines: ",nrow(sm_obj),"\n")

  
out <- gzfile(output.file, "w")
write.csv(sm_obj, out, row.names=F)
close(out)
