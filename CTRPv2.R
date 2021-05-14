library(downloader)
library(PharmacoGx)
library(devtools)
library(Biobase)
library(data.table)
library(reshape2)
library(CoreGx)
library(SummarizedExperiment)
library(biocompute)

options(stringsAsFactors=FALSE)

args = commandArgs(trailingOnly=TRUE)
ORCESTRA_ID <- args

badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]|[(]|[)]"

cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))

rownames(cell_all) <- cell_all$unique.cellid

curationCell <- cell_all[which(!is.na(cell_all[ , "CTRPv2.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "CTRPv2.cellid")]

rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "CTRPv2.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "CTRPv2.drugid","smiles")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]

standardize <- args[grep("filtered", args)]

standardizeRawDataConcRange <- function(sens.info, sens.raw){
  unq.drugs <- unique(sens.info$drugid)
  
  conc.m <- data.table(melt(sens.raw[,,1], as.is=TRUE))
  conc.m[,drugid := sens.info$drugid[match(Var1, rownames(sens.info))]]
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  conc.ranges[,Var1 := NULL]
  conc.ranges <- conc.ranges[,unique(.SD), drugid]	
  # conc.ranges[,N := .N, drugid]
  conc.ranges.disj <- conc.ranges[, {sq <- sort(unique(c(l,r))); 
  l = sq[seq(1,length(sq)-1)];
  r = sq[seq(2,length(sq))];
  .(l=l,r=r)}, drugid]
  ## Function below returns all consecutive ranges of ints between 1 and N
  returnConsInts <- function(N) {
    stopifnot(N>0)
    unlist(sapply(seq(1,N), function(ii) return(sapply(seq(ii, N), function(jj) return(seq(ii,jj))))), recursive=FALSE)
  }
  rangeNoHoles <- function(indicies, lr.tbl){
    if(length(indicies) == 1) return(TRUE)
    sq <- seq(indicies[1], indicies[length(indicies)]-1)
    all(lr.tbl[["l"]][sq+1] <= lr.tbl[["r"]][sq])
  }
  per.drug.range.indicies <- sapply(conc.ranges.disj[,.N,drugid][,N], returnConsInts)
  
  names(per.drug.range.indicies) <- conc.ranges.disj[,unique(drugid)] ## checked this: conc.ranges.disj[,.N,drugid][,drugid] == conc.ranges.disj[,unique(drugid)]
  
  
  # Check if there are any holes in the chosen range combination
  per.drug.range.indicies <- sapply(names(per.drug.range.indicies), function(drug){
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    per.drug.range.indicies[[drug]][sapply(per.drug.range.indicies[[drug]], rangeNoHoles, lr.tbl = lr.tbl)]
    
  })
  per.drug.range.indicies.2 <- sapply(names(per.drug.range.indicies), function(drug){
    
    lr.tbl <- conc.ranges.disj[drugid == drug]
    res <- t(sapply(per.drug.range.indicies[[drug]], function(x) return(c(lr.tbl[x[1],l], lr.tbl[x[length(x)],r]))))
    colnames(res) <- c("l", "r")
    res <- data.frame(res)
    res <- cbind(drugid = drug, res)
  }, simplify=FALSE)
  per.drug.range.indicies.dt <- rbindlist(per.drug.range.indicies.2)
  
  conc.ranges <- conc.m[,.(l = min(value, na.rm=T), r = max(value, na.rm=T)), c("drugid", "Var1")]
  setkey(conc.m, Var1)
  conc.m <- na.omit(conc.m)
  setkey(conc.m, drugid, Var1, value)
  setkey(conc.ranges, drugid, l, r)
  # tic()
  ## NOTE:: Data.table used for maximum speed. Probably possible to do this more intelligently by 
  ## NOTE:: being aware of which conditions overlap, but its fast enough right now as it is.
  chosen.drug.ranges <- lapply(unq.drugs, function(drug){
    num.points.in.range <- apply(per.drug.range.indicies.dt[drugid==drug, .(l,r)], 1, function(rng){
      conc.m[drugid==drug][conc.ranges[drugid==drug][l<=rng["l"]][r>=rng["r"],Var1], on="Var1"][value >= rng["l"]][value <= rng["r"],.N]
      # conc.m[drugid==drug][, Var1]
    })
    max.ranges <- per.drug.range.indicies.dt[drugid==drug][which(num.points.in.range==max(num.points.in.range))]
    max.ranges[which.max(log10(r) - log10(l)), ]
  })
  # toc()
  names(chosen.drug.ranges) <- sapply(chosen.drug.ranges, `[[`, "drugid")
  removed.experiments <- unlist(lapply(unq.drugs, function(drug){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    exp.out.range <- conc.ranges[drugid==drug][l>rng["l"] | r<rng["r"],Var1]
    return(exp.out.range)
  }))
  
  sens.raw[removed.experiments,,] <- NA_real_
  conc.ranges.kept <- conc.ranges[!Var1 %in% removed.experiments]
  
  for(drug in unq.drugs){
    rng <- unlist(chosen.drug.ranges[[drug]][,.(l,r)])
    myx <- conc.ranges.kept[drugid==drug,Var1]
    doses <- sens.raw[myx, ,"Dose"]
    which.remove <- (doses < rng["l"] | doses > rng["r"])
    sens.raw[myx, ,"Dose"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    sens.raw[myx, ,"Viability"][which(which.remove,arr.ind=TRUE)] <- NA_real_
    
    ## Annotate sens info with chosen range
    sens.info[sens.info$drugid==drug,"chosen.min.range"] <- rng["l"]
    sens.info[sens.info$drugid==drug,"chosen.max.range"] <- rng["r"]
  }
  sens.info$rm.by.conc.range <- FALSE
  sens.info[removed.experiments,"rm.by.conc.range"] <- TRUE
  
  return(list("sens.info" = sens.info, sens.raw = sens.raw))
}


#filter noisy curves from PSet (modified function to take into account standardized conc range)
filterNoisyCurves2 <- function(pSet, epsilon=25 , positive.cutoff.percent=.80, mean.viablity=200, nthread=1) {
  acceptable <- mclapply(rownames(sensitivityInfo(pSet)), function(xp) {
    #for(xp in rownames(sensitivityInfo(pSet))){
    drug.responses <- as.data.frame(apply(pSet@sensitivity$raw[xp , ,], 2, as.numeric), stringsAsFactors=FALSE)
    if (!all(is.na(drug.responses))){
      
      
      drug.responses <- drug.responses[complete.cases(drug.responses), ]
      doses.no <- nrow(drug.responses)
      drug.responses[,"delta"] <- .computeDelta(drug.responses$Viability)
      
      delta.sum <- sum(drug.responses$delta, na.rm = TRUE)
      
      max.cum.sum <- .computeCumSumDelta(drug.responses$Viability)
      
      if ((table(drug.responses$delta < epsilon)["TRUE"] >= (doses.no * positive.cutoff.percent)) &
          (delta.sum < epsilon) &
          (max.cum.sum < (2 * epsilon)) &
          (mean(drug.responses$Viability) < mean.viablity)) {
        return (xp)
      }
    }
    
  }, mc.cores=nthread)
  acceptable <- unlist(acceptable)
  noisy <- setdiff(rownames(sensitivityInfo(pSet)), acceptable)
  return(list("noisy"=noisy, "ok"=acceptable))
}

.computeDelta <- function(xx ,trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc)
  {
    return(c(pmin(100, xx[2:length(xx)]) - pmin(100, xx[1:length(xx)-1]), 0))
  }else{
    return(c(xx[2:length(xx)] - xx[1:length(xx)-1]), 0)
  }
}

#' @importFrom utils combn
.computeCumSumDelta <- function(xx, trunc = TRUE) {
  xx <- as.numeric(xx)
  if(trunc) {
    xx <- pmin(xx, 100)
  }
  tt <- t(combn(1:length(xx), 2 , simplify = TRUE))
  tt <- tt[which(((tt[,2] - tt[,1]) >= 2) == TRUE),]
  if (is.null(nrow(tt))){
    tt <- matrix(tt, ncol = 2)
  }
  cum.sum <- unlist(lapply(1:nrow(tt), function(x){xx[tt[x,2]]-xx[tt[x,1]]}))
  return(max(cum.sum))
}   

matchToIDTable <- function(ids,tbl, column, returnColumn="unique.cellid") {
  sapply(ids, function(x) {
                          myx <- grep(paste0("((///)|^)",Hmisc::escapeRegex(x),"((///)|$)"), tbl[,column])
                          if(length(myx) > 1){
                            stop("Something went wrong in curating ids, we have multiple matches")
                          }
        if(length(myx) == 0){return(NA_character_)}
                          return(tbl[myx, returnColumn])
                        })
}
  
  
  load("/pfs/CTRPv2RawSensitivity/drug_post.RData")
  load("/pfs/ctrpv2ProfilesAssemble/profiles.RData")
  
  recomputed <- res
  
    
    
# ctrp.sensitivityInfo <- read.delim("/pfs/ctrpv2raw/v20.meta.per_experiment.txt")

# repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])


# for(exp in repExps){

#   myx <- ctrp.sensitivityInfo$experiment_id == exp
#   duplicates <- duplicated(ctrp.sensitivityInfo$experiment_id) & myx
#   first <- myx & !duplicates
  
#   ctrp.sensitivityInfo[first,] <- apply(ctrp.sensitivityInfo[myx,], 2, function(x) paste(unique(x), collapse="//"))

#   ctrp.sensitivityInfo <- ctrp.sensitivityInfo[!duplicates,]
# }

# ctrp.sensitivityInfo[,"cellid"] <- ctrp.cells$cellid[match(ctrp.sensitivityInfo$master_ccl_id, ctrp.cells$master_ccl_id)]
                                        
                                        
  #######################################################
###### Published AUC and EC50 values 
#######################################################

ctrp.sensitivityPub <- read.delim("/pfs/ctrpv2raw/v20.data.curves_post_qc.txt")


ctrp.sensitivityPub[,"cellid"] <- sensitivityInfo[match(ctrp.sensitivityPub$experiment_id, sensitivityInfo$experiment_id), "cellid"] 
ctrp.sensitivityPub[,"drugid"] <- ctrp.drugs$drugid[match(ctrp.sensitivityPub$master_cpd_id, ctrp.drugs$master_cpd_id)] 
ctrp.sensitivityPub[,"culture_media"] <- sensitivityInfo[match(ctrp.sensitivityPub$experiment_id, sensitivityInfo$experiment_id), "culture_media"] 

experimentIdsPub <- paste(ctrp.sensitivityPub[,"cellid"], ctrp.sensitivityPub[,"drugid"],ctrp.sensitivityPub[,"culture_media"], ctrp.sensitivityPub[,"experiment_id"], sep="_")
ctrp.sensitivityPub$experimentIds <- experimentIdsPub

rownames(ctrp.sensitivityPub) <- ctrp.sensitivityPub$experimentIds 

sensitivityProfiles <- data.frame("aac_published" = ctrp.sensitivityPub$area_under_curve, "ec50_published"=ctrp.sensitivityPub$apparent_ec50_umol)

rownames(sensitivityProfiles) <- ctrp.sensitivityPub$experimentIds 

stopifnot(setequal(rownames(sensitivityProfiles), rownames(raw.sensitivity)))

stopifnot(setequal(rownames(sensitivityProfiles), rownames(sensitivityInfo)))


raw.sensitivity <- raw.sensitivity[rownames(sensitivityInfo),,]

sensitivityProfiles <- sensitivityProfiles[rownames(sensitivityInfo),]

#sensitivityProfiles[,c("auc_recomputed", "ic50_recomputed", "Hill-Slope", "E_inf", "EC50")] <- cbind(sensitivityProfiles, recomputed[,"auc_recomputed"][rownames(sensitivityProfiles)], recomputed[,"ic50_recomputed"][rownames(sensitivityProfiles)], recomputed[,"HS"][rownames(sensitivityProfiles)], recomputed[,"E_inf"][rownames(sensitivityProfiles)],recomputed[,"EC50"][rownames(sensitivityProfiles)] )

sensitivityProfiles <- cbind(sensitivityProfiles, recomputed[rownames(sensitivityProfiles),])
sensitivityProfiles[,"aac_recomputed"] <- as.numeric(sensitivityProfiles[,"aac_recomputed"])
sensitivityProfiles[,"ic50_recomputed"] <- as.numeric(sensitivityProfiles[,"ic50_recomputed"])
sensitivityProfiles[,"HS"] <- as.numeric(sensitivityProfiles[,"HS"])
sensitivityProfiles[,"E_inf"] <- as.numeric(sensitivityProfiles[,"E_inf"])
sensitivityProfiles[,"EC50"] <- as.numeric(sensitivityProfiles[,"EC50"])



## ctrp.cells are not all in sensitivity info. filtering only to relevant ones

ctrp.cells <- ctrp.cells[ctrp.cells$cellid %in% sensitivityInfo$cellid, ]


#MATCH SENSITIVITY INFO CELL/DRUG ID TO UNIQUE-ID (Both published and recomputed)

## Cell lines first
## Doing this in two steps because regex is not as efficient as exact string match

mapTable <- cbind(unique(sensitivityInfo$cellid), matchToIDTable(ids=unique(sensitivityInfo$cellid), tbl=curationCell, column = "CTRPv2.cellid", returnColumn="unique.cellid"))

reps <- mapTable[match(sensitivityInfo$cellid, mapTable[,1]),2]

stopifnot(!anyNA(reps))
sensitivityInfo$cellid <- reps


reps <- matchToIDTable(ids=unique(ctrp.cells$cellid), tbl=curationCell, column = "CTRPv2.cellid", returnColumn="unique.cellid")
stopifnot(!anyNA(reps))

ctrp.cells$cellid <- reps



## Now drugs

mapTable <- cbind(unique(sensitivityInfo$drugid), matchToIDTable(ids=unique(sensitivityInfo$drugid), tbl=curationDrug, column = "CTRPv2.drugid", returnColumn="unique.drugid"))

reps <- mapTable[match(sensitivityInfo$drugid, mapTable[,1]),2]

stopifnot(!anyNA(reps))
sensitivityInfo$drugid <- reps


reps <- matchToIDTable(ids=unique(ctrp.drugs$drugid), tbl=curationDrug, column = "CTRPv2.drugid", returnColumn="unique.drugid")
stopifnot(!anyNA(reps))

ctrp.drugs$drugid <- reps

# #tipifarnib-P1 & P2 have same unique id. Combined them and separated by /// into one row (row names need unique ids)
# sensitivityInfo$drugid[sensitivityInfo$drugid == "Tipifarnib-P1"] <- "Tipifarnib"
# sensitivityInfo$drugid[sensitivityInfo$drugid == "Tipifarnib-P2"] <- "Tipifarnib"
# drug_match <- ctrp.drugs$drugid
# ctrp.drugs[218,] <- paste0(ctrp.drugs[218,], sep="///" , ctrp.drugs[219,])
# ctrp.drugs <- ctrp.drugs[-219,]
# drug_match[218] <- "Tipifarnib"
# drug_match <- drug_match[-219]
# ctrp.drugs$drugid <- drug_match

collapseRows4 <- function(x, rows, skip.cols = c("unique.cellid", "PharmacoDB.id", "unique.tissueid")){
    skip.cols <-  colnames(x) %in% skip.cols
    xNew <- lapply(x[rows, !skip.cols], function(x) {
       xx <- na.omit(x)
       if (length(xx) == 0) {
         xx <- NA
       }
       if (length(unique(xx)) > 1) {
         xx <- paste(xx, collapse="///")
       } else {xx <- xx[1]}
       return(as.vector(xx))
       })
     xNew <- as.data.frame(xNew, as.is = TRUE)
     x[rows[1], !skip.cols] <- xNew
     x <- x[-rows[-1], ]
     return(x)
}

# #tipifarnib-P1 & P2 have same unique id. Combined them and separated by /// into one row (row names need unique ids)

ctrp.drugs <- collapseRows4(ctrp.drugs, which(ctrp.drugs$drugid == "Tipifarnib"), skip.cols=character())

rownames(ctrp.drugs) <- ctrp.drugs$drugid

rownames(ctrp.cells) <- ctrp.cells$cellid



curationTissue <- data.frame(unique.tissueid = cell_all[curationCell$unique.cellid, "unique.tissueid"], "CTRPv2.tissueid" = ctrp.cells$tissueid[match(curationCell$unique.cellid,rownames(ctrp.cells))])
rownames(curationTissue) <- curationCell$unique.cellid


ctrp.cells$tissueid <- curationTissue[ctrp.cells$cellid,"unique.tissueid"]



emptyE <- ExpressionSet()
pData(emptyE)$cellid <- character()
pData(emptyE)$batchid <- character()
fData(emptyE)$BEST <- vector()
fData(emptyE)$Symbol <- character()
annotation(emptyE) <- "CTRP contains no molecular profiles of cell lines. This SE is empty placeholder."

emptySE <- SummarizedExperiment::SummarizedExperiment(
  ## TODO:: Do we want to pass an environment for better memory efficiency?
  assays=S4Vectors::SimpleList(as.list(Biobase::assayData(emptyE))
  ),
  # Switch rearrange columns so that IDs are first, probes second
  rowData=S4Vectors::DataFrame(Biobase::fData(emptyE),
                               rownames=rownames(Biobase::fData(emptyE)) 
  ),
  colData=S4Vectors::DataFrame(Biobase::pData(emptyE),
                               rownames=rownames(Biobase::pData(emptyE))
  ),
  metadata=list("experimentData" = emptyE@experimentData, 
                "annotation" = Biobase::annotation(emptyE), 
                "protocolData" = Biobase::protocolData(emptyE)
  )
)

cellsPresent <- sort(CoreGx::.unionList(sensitivityInfo$cellid))    
ctrp.cells <- ctrp.cells[cellsPresent,]

ctrp.cells$tissueid <- curationTissue[rownames(ctrp.cells), "unique.tissueid"]
ctrp.cells$cellid <- rownames(ctrp.cells)

curationTissue <- curationTissue[rownames(ctrp.cells),]
curationCell <- curationCell[rownames(ctrp.cells),]


drugsPresent <- sort(unique(sensitivityInfo$drugid))
ctrp.drugs <- ctrp.drugs[drugsPresent,]
    
drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))
drug_all <- drug_all[which(!is.na(drug_all[ , "CTRPv2.drugid"])),]
drug_all <- drug_all[ , c("unique.drugid", "CTRPv2.drugid","smiles","inchikey","cid","FDA")]
rownames(drug_all) <- drug_all[ , "unique.drugid"]

drug_all <- drug_all[rownames(ctrp.drugs),]
ctrp.drugs[,c("smiles","inchikey","cid","FDA")] <- drug_all[,c("smiles","inchikey","cid","FDA")]

		 
#add cellosaurus disease type to cell-info
		 

disease <- cell_all$Cellosaurus.Disease.Type[match(ctrp.cells$cellid, cell_all$unique.cellid)]
ctrp.cells$Cellosaurus.Disease.Type <- disease		 
	
#add cellosaurus assession to cell-info
assession <- cell_all$Cellosaurus.Accession.id[match(ctrp.cells$cellid, cell_all$unique.cellid)]
ctrp.cells$Cellosaurus.Accession.id <- assession
		 
#add pharmacodb id to cell-info
pdb <- cell_all$PharmacoDB.id[match(ctrp.cells$cellid, cell_all$unique.cellid)]
ctrp.cells$PharmacoDB.id <- pdb

#add study tissue id to cell_info
study_tissue <- cell_all$unique.tissueid.fromstudies[match(ctrp.cells$cellid, cell_all$unique.cellid)]
ctrp.cells$unique.tissueid.fromstudies <- study_tissue
		 
#add study cell-line type to cell_info
cell_type <- cell_all$CellLine.Type[match(ctrp.cells$cellid, cell_all$unique.cellid)]
ctrp.cells$CellLine.Type <- cell_type
		 
#add metastatic info to cell_info		 
metastatic <- cell_all$Metastatic[match(ctrp.cells$cellid, cell_all$unique.cellid)]
ctrp.cells$Metastatic <- metastatic		 
curationDrug <- curationDrug[rownames(ctrp.drugs),]

if (length(standardize) > 0){

# standardize <- standardizeRawDataConcRange(sens.info = sensitivityInfo, sens.raw = raw.sensitivity)
# sensitivityInfo <- standardize$sens.info
# raw.sensitivity <- standardize$sens.raw

} else {
print("unfiltered PSet")
	
}
                 
CTRPv2 <- PharmacoGx::PharmacoSet(name="CTRPv2", 
 molecularProfiles = list("rna" = emptySE),
 cell=ctrp.cells, 
 drug=ctrp.drugs, 
 sensitivityInfo=sensitivityInfo, 
 sensitivityRaw=raw.sensitivity, 
 sensitivityProfiles=sensitivityProfiles, 
 sensitivityN=NULL,  
 curationCell=curationCell, 
 curationDrug=curationDrug,
 curationTissue=curationTissue, 
 datasetType="sensitivity",
 verify = TRUE)
 
if (length(standardize) > 0){

 noisy_out <- filterNoisyCurves2(CTRPv2)
 print("filter done")
 CTRPv2@sensitivity$profiles[noisy_out$noisy, ] <- NA

} else {
print("unfiltered PSet")
	
}                 

CTRPv2@annotation$version <- 2
saveRDS(CTRPv2,file="/pfs/out/CTRPv2.rds")

dataset <- "CTRPv2"		 
#output ORCESTRA_ID and Pachyderm commit id
write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)				   
pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F) 
 
  
