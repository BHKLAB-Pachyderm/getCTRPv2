library(downloader)
library(PharmacoGx)
library(devtools)

getCTRPv2P <-
  function (
    verbose=FALSE,
    nthread=1){


options(stringsAsFactors=FALSE)
badchars <- "[\xb5]|[]|[ ,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[\\]|[.]|[_]|[ ]|[(]|[)]"

cell_all <- read.csv("/pfs/downAnnotations/cell_annotation_all.csv", na.strings=c("", " ", "NA"))
drug_all <- read.csv("/pfs/downAnnotations/drugs_with_ids.csv", na.strings=c("", " ", "NA"))

curationCell <- cell_all[which(!is.na(cell_all[ , "CTRPv2.cellid"])),]
curationCell <- curationCell[ , c("unique.cellid", "CTRPv2.cellid")]

rownames(curationCell) <- curationCell[ , "unique.cellid"]

curationDrug <- drug_all[which(!is.na(drug_all[ , "CTRPv2.drugid"])),]
curationDrug <- curationDrug[ , c("unique.drugid", "CTRPv2.drugid","smiles")]
rownames(curationDrug) <- curationDrug[ , "unique.drugid"]


matchToIDTableCELL <- function(ids,tbl, column) {
  sapply(ids, function(x) {
    myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
    if(length(myx) > 1){
      stop("Something went wrong in curating cell ids")
    }
    return(tbl[myx, "unique.cellid"])
  })
}


matchToIDTableDRUG <- function(ids,tbl, column) {
    sapply(ids,function(x) {
      myx <- grep(paste0("((///)|^)",x,"((///)|$)"), tbl[,column])
      if(length(myx) > 1){
        stop("Something went wrong in curating drug ids")
      }
      return(tbl[myx, "unique.drugid"])
    })
  }
  
  
  load("/pfs/CTRPv2RawSensitivity/drug_post.RData")
  load("/pfs/ctrpv2ProfilesAssemble/profiles.RData")
  
  recomputed <- res
  
    
    
ctrp.sensitivityInfo <- read.delim("/pfs/ctrpv2raw/v20.meta.per_experiment.txt")

repExps <- unique(ctrp.sensitivityInfo$experiment_id[duplicated(ctrp.sensitivityInfo$experiment_id)])


for(exp in repExps){

  myx <- ctrp.sensitivityInfo$experiment_id == exp
  duplicates <- duplicated(ctrp.sensitivityInfo$experiment_id) & myx
  first <- myx & !duplicates
  
  ctrp.sensitivityInfo[first,] <- apply(ctrp.sensitivityInfo[myx,], 2, function(x) paste(unique(x), collapse="//"))

  ctrp.sensitivityInfo <- ctrp.sensitivityInfo[!duplicates,]
}

ctrp.sensitivityInfo[,"cellid"] <- ctrp.cells$cellid[match(ctrp.sensitivityInfo$master_ccl_id, ctrp.cells$master_ccl_id)]
                                        
                                        
  #######################################################
###### Published AUC and EC50 values 
#######################################################

ctrp.sensitivityPub <- read.delim("/pfs/ctrpv2raw/v20.data.curves_post_qc.txt")


ctrp.sensitivityPub[,"cellid"] <- ctrp.sensitivityInfo[match(ctrp.sensitivityPub$experiment_id, ctrp.sensitivityInfo$experiment_id), "cellid"] ## matches are sometimes duplicated, but the values for cells and drugs are unique either way, so take the first match
ctrp.sensitivityPub[,"drugid"] <- ctrp.drugs$drugid[match(ctrp.sensitivityPub$master_cpd_id, ctrp.drugs$master_cpd_id)] 
ctrp.sensitivityPub[,"culture_media"] <- ctrp.sensitivityInfo[match(ctrp.sensitivityPub$experiment_id, ctrp.sensitivityInfo$experiment_id), "culture_media"] 

experimentIdsPub <- paste(ctrp.sensitivityPub[,"cellid"], ctrp.sensitivityPub[,"drugid"],ctrp.sensitivityPub[,"culture_media"], ctrp.sensitivityPub[,"experiment_id"], sep="_")
ctrp.sensitivityPub$experimentIds <- experimentIdsPub

rownames(ctrp.sensitivityPub) <- ctrp.sensitivityPub$experimentIds 

sensitivityProfiles <- data.frame("auc_published" = ctrp.sensitivityPub$area_under_curve, "ec50_published"=ctrp.sensitivityPub$apparent_ec50_umol)

rownames(sensitivityProfiles) <- ctrp.sensitivityPub$experimentIds 

raw.sensitivity <- raw.sensitivity[rownames(sensitivityInfo),,]

sensitivityProfiles <- sensitivityProfiles[rownames(sensitivityInfo),]

#sensitivityProfiles[,c("auc_recomputed", "ic50_recomputed", "Hill-Slope", "E_inf", "EC50")] <- cbind(sensitivityProfiles, recomputed[,"auc_recomputed"][rownames(sensitivityProfiles)], recomputed[,"ic50_recomputed"][rownames(sensitivityProfiles)], recomputed[,"HS"][rownames(sensitivityProfiles)], recomputed[,"E_inf"][rownames(sensitivityProfiles)],recomputed[,"EC50"][rownames(sensitivityProfiles)] )

sensitivityProfiles <- cbind(sensitivityProfiles, recomputed[rownames(sensitivityProfiles),])
sensitivityProfiles[,"AAC"] <- as.numeric(sensitivityProfiles[,"AAC"])
sensitivityProfiles[,"IC50"] <- as.numeric(sensitivityProfiles[,"IC50"])
sensitivityProfiles[,"HS"] <- as.numeric(sensitivityProfiles[,"HS"])
sensitivityProfiles[,"E_inf"] <- as.numeric(sensitivityProfiles[,"E_inf"])
sensitivityProfiles[,"EC50"] <- as.numeric(sensitivityProfiles[,"EC50"])


#MATCH SENSITIVITY INFO CELL/DRUG ID TO UNIQUE-ID (Both published and recomputed)


#tipifarnib-P1 & P2 have same unique id. Combined them and separated by /// into one row (row names need unique ids)
sensitivityInfo$drugid[sensitivityInfo$drugid == "Tipifarnib-P1"] <- "Tipifarnib"
sensitivityInfo$drugid[sensitivityInfo$drugid == "Tipifarnib-P2"] <- "Tipifarnib"
ctrp.drugs[218,] <- paste0(ctrp.drugs[218,], sep="///" , ctrp.drugs[219,])
ctrp.drugs <- ctrp.drugs[-219,]
drug_match[218] <- "Tipifarnib"
drug_match <- drug_match[-219]
ctrp.drugs$drugid <- drug_match
rownames(ctrp.drugs) <- ctrp.drugs$drugid

rownames(ctrp.cells) <- ctrp.cells$cellid
curationTissue <- cbind("unique.tissueid"=ctrp.cells$tissueid, "CTRPv2.tissueid"=ctrp.cells$tissueid)

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
 verify = FALSE)
 
 save(CTRPv2,file="/pfs/out/CTRPv2.RData")
    
 return (CTRPv2)
 
 }
 
 getCTRPv2P(verbose=FALSE, nthread=1)
 
  
