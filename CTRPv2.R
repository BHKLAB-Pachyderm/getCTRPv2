library(downloader)
library(PharmacoGx)
library(devtools)
library(Biobase)

getCTRPv2P <-
  function (
    verbose=FALSE,
    nthread=1){


options(stringsAsFactors=FALSE)
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

sensitivityProfiles <- data.frame("auc_published" = ctrp.sensitivityPub$area_under_curve, "ec50_published"=ctrp.sensitivityPub$apparent_ec50_umol)

rownames(sensitivityProfiles) <- ctrp.sensitivityPub$experimentIds 

stopifnot(setequal(rownames(sensitivityProfiles), rownames(raw.sensitivity)))

stopifnot(setequal(rownames(sensitivityProfiles), rownames(sensitivityInfo)))


raw.sensitivity <- raw.sensitivity[rownames(sensitivityInfo),,]

sensitivityProfiles <- sensitivityProfiles[rownames(sensitivityInfo),]

#sensitivityProfiles[,c("auc_recomputed", "ic50_recomputed", "Hill-Slope", "E_inf", "EC50")] <- cbind(sensitivityProfiles, recomputed[,"auc_recomputed"][rownames(sensitivityProfiles)], recomputed[,"ic50_recomputed"][rownames(sensitivityProfiles)], recomputed[,"HS"][rownames(sensitivityProfiles)], recomputed[,"E_inf"][rownames(sensitivityProfiles)],recomputed[,"EC50"][rownames(sensitivityProfiles)] )

sensitivityProfiles <- cbind(sensitivityProfiles, recomputed[rownames(sensitivityProfiles),])
sensitivityProfiles[,"aac_recomputed"] <- as.numeric(sensitivityProfiles[,"AAC"])
sensitivityProfiles[,"ic50_recomputed"] <- as.numeric(sensitivityProfiles[,"IC50"])
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



curationTissue <- data.frame(unique.tissueid = cell_all[curationCell$unique.cellid, "unique.tissueid"], "CTRPv2.tissueid" = ctrp.cells$tissueid[curationCell$unique.cellid])
rownames(curationTissue) <- curationCell$unique.cellid


ctrp.cells$tissueid <- curationTissue[ctrp.cells$cellid,]



emptyEset <- ExpressionSet()
annotation(emptyEset) <- "CTRP contains no molecular profiles of cell lines. Please use data from other datasets. This eset is empty placeholder."


CTRPv2 <- PharmacoSet(name="CTRPv2", 
 molecularProfiles = list("rna" = emptyEset),
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
 
 saveRDS(CTRPv2,file="/pfs/out/CTRPv2.rds")
    
 # return (CTRPv2)
 
 }
 
 getCTRPv2P(verbose=FALSE, nthread=1)
 
  
