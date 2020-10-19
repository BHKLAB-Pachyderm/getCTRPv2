library(downloader)
library(PharmacoGx)
library(devtools)
library(Biobase)
library(data.table)
library(reshape2)
library(CoreGx)
library(SummarizedExperiment)


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
 
                 

CTRPv2@annotation$version <- 2
saveRDS(CTRPv2,file="/pfs/out/CTRPv2.rds")

dataset <- "CTRPv2"		 
#output ORCESTRA_ID and Pachyderm commit id
write.table(dataset, file="/pfs/out/dataset.txt", row.names = F ,quote = F, sep = "\t", col.names = F)
write.table(ORCESTRA_ID, file="/pfs/out/orcestra_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F)				   
pach_commit_id <- Sys.getenv("PACH_OUTPUT_COMMIT_ID")
write.table(pach_commit_id, file="/pfs/out/commit_id.txt", row.names = F ,quote = F, sep = "\t", col.names = F) 
 
  
