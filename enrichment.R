#### Enrichment 

library(openxlsx)
library(biomaRt)
library(rjson)
library(httr)


enrich.fishers <- function(exp.paths, obs.paths, fdr=T, pathway.min.support=5){
  # count pathways and store unique entries w/ counts
  exp.counts <- as.data.frame(table(exp.paths))
  colnames(exp.counts) <- c("PATHWAY", "exp.freq")
  obs.counts <- as.data.frame(table(obs.paths))
  colnames(obs.counts) <- c("PATHWAY", "obs.freq")
  
  # merge unique pathway counts and set missing counts to 0
  obs.counts <- merge(exp.counts, obs.counts, by = "PATHWAY", all.x = TRUE)
  #names(obs.counts) <- c("PATHWAY",  "exp.freq" , "obs.freq")
  obs.counts[is.na(obs.counts$obs.freq),"obs.freq"] <- 0
  
  # prepare and empty column for p-vals from Fisher's Exact Test (FET)
  obs.counts$fisher.p <- NA
  
  # iterate through pathways and perform FET on each
  for(p in obs.counts$PATHWAY){
    in.path.mask <- obs.counts$PATHWAY == p
    
    # skip pathways with less than the minimal support
    if (obs.counts$exp.freq[in.path.mask] < pathway.min.support){
      obs.counts[in.path.mask, "fisher.p"] <- NA
      next()
    }

    # contingency table values for the current pathway
    obs.in.path <- obs.counts$obs.freq[in.path.mask]  # observed & in pathway
    obs.no.path <- sum(obs.counts$obs.freq) - obs.in.path         # observed but not in pathway
    no.obs.in.path <- obs.counts$exp.freq[in.path.mask] - obs.in.path # not observed & in pathway
    no.obs.no.path <- sum(obs.counts$exp.freq[!in.path.mask]) - obs.no.path # not observed but not in pathway
    
    # generate contingency table and run FET
    cntg_tab <- matrix(ncol = 2, nrow = 2, c(obs.in.path, obs.no.path, no.obs.in.path, no.obs.no.path))
    obs.counts[in.path.mask, "fisher.p"] <- fisher.test(cntg_tab, simulate.p.value = F, alternative = "greater")$p.value
  }
  
  # order pathways by ascending p-value
  obs.counts <- obs.counts[order(obs.counts$fisher.p),]
  
  # (optional) adjust p-values using FDR
  if (fdr){
    obs.counts$fdr.p <- p.adjust(obs.counts$fisher.p, method = "fdr")
  }
  
  return(obs.counts)
}

gene.ids2symbols <- function(id.list){
  ensembl = useEnsembl(biomart="ensembl",dataset="hsapiens_gene_ensembl", mirror="useast",)
  symbls = getBM(attributes=c("ensembl_gene_id", "hgnc_symbol"), 
                 filters = 'ensembl_gene_id', 
                 values = id.list, 
                 mart = ensembl)
  return(symbls)
}


gene.annotate.from <- function(sym.list, annotation.db) {
  # blank list to hold annotations
  annots <- as.data.frame(sym.list)
  annots$PATHWAY <- NA
  names(annots) <- c("SYMBOL", "PATHWAY")
  
  # open file and iterate through annotations
  con = file(annotation.db, "r")
  on.exit(close(con))
  while ( TRUE ) {
    line = readLines(con, n = 1)
    
    if ( length(line) == 0 ) {
      break
    }
    
    # EnrichR library contents are tab separated
    line.contents = unlist(strsplit(line, "\t"))
    
    # search for matching Entrez symbols
    for (sym in line.contents[3:length(line.contents)]){
      if (sym %in% sym.list){
        if (anyNA(annots[annots$SYMBOL == sym, "PATHWAY"])){
          annots[annots$SYMBOL == sym, "PATHWAY"] <- line.contents[1]
        }else {
          new.row <- data.frame(SYMBOL=sym, PATHWAY=line.contents[1])
          rbind(annots, new.row)
        }
      }
    }
  }
  
  annots[is.na(annots$PATHWAY), "PATHWAY"] <- 'Unannotated'
  
  return(annots)
}

gene.enrichR <- function(exp.list,
                         description = "ClustOmics enrichment analysis", 
                         gene.set.lib='KEGG_2019_Human', export.to=NA){
  ### ADD GENE LIST ###
  ENRICHR_ADD_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
  genes.str = paste(exp.list, "", collapse="\n")
  
  payload = list('list'=c(NULL, genes.str), 'description'=c(NULL, description))
  response.add = POST(ENRICHR_ADD_URL, body=payload, encode = "multipart")
  stop_for_status(response.add)
  
  data.add = fromJSON(content(response.add, "text")) # for user ID
  
  ### PERFORM ENRICHMENT ###
  ENRICHR_ENRICH_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
  query.string = '?userListId=%s&backgroundType=%s'
  user.list.id = data.add$userListId
  
  url = paste0(c(ENRICHR_ENRICH_URL, sprintf(query.string, user.list.id, gene.set.lib)), "", collapse = "")
  response.enrich = GET(url)
  stop_for_status(response.enrich)
  
  data.enrich = fromJSON(content(response.enrich, "text")) # for results
  
  ### EXPORT ENRICHMENT (IF APPLICABLE) ###
  if (!is.na(export.to)){
    ENRICHR_EXPORT_URL = 'http://amp.pharm.mssm.edu/Enrichr/export'
    query.string = '?userListId=%s&filename=%s&backgroundType=%s'
    
    url = paste0(c(ENRICHR_EXPORT_URL, sprintf(query.string, user.list.id, export.to, gene.set.lib)), "", collapse = "")
    response.export = GET(url, stream=T)
    stop_for_status(response.export)
    
    bin <- content(response.export, "raw")
    writeBin(bin, paste0(export.to, ".txt"))
  }
  
  return(data.enrich)
}

met.annotate.from <- function(name.list, annotation.db) {
  # make.names used to force both character vectors to match
  name.list.1 <- make.names(name.list)
  annots <- data.frame(BIOCHEMICAL=make.names(name.list.1))
  
  metadata <- read.delim(file = annotation.db, sep = ",", stringsAsFactors = FALSE)
  metadata$BIOCHEMICAL <- make.names(metadata$BIOCHEMICAL)
  metadata[is.na(metadata$SUB.PATHWAY), 'SUB.PATHWAY'] <- 'Unannotated'
  metadata[is.na(metadata$SUPER.PATHWAY), 'SUPER.PATHWAY'] <- 'Unannotated'
  
  met <- merge(annots, metadata, by.x = "BIOCHEMICAL", by.y = "BIOCHEMICAL") #, all.x = TRUE)
  met <- met[, c("BIOCHEMICAL", "SUPER.PATHWAY", "SUB.PATHWAY")]
  met$PATHWAY <- met$SUPER.PATHWAY
  
  return(met)
}

prot.annotate.from <- function(){

soma_feat <- read.xlsx("./cluster_associations/results_9-3/Feature_Rank_SOMA.xlsx")
prots <- soma$'Protein/Metabolite'


meta <- read.csv("./COPDGene/Data/somalogic/Soma-ProteinNames.csv", row.names = 1)
meta <- as.data.frame(t(meta))



mart <- biomaRt::useDataset(dataset = "hsapiens_gene_ensembl",         
                    mart    = useMart("ENSEMBL_MART_ENSEMBL"))  
att <- listAttributes(mart)
att[grep("prot", listAttributes(mart), ignore.case = T), ]
write.csv(att, file = "./biomart_att.csv", row.names = FALSE)
    
resultTable <- getBM(attributes = c("entrezgene_id", "goslim_goa_description"),       
                      filters    = "entrezgene_id",       
                      values     = as.character(meta[,"EntrezGeneID"]),         
                      mart       = mart)        

#resultTable <- resultTable[!duplicated(resultTable$entrezgene_id),]
unique(resultTable$goslim_goa_description)[grep("extra", unique(resultTable$goslim_goa_description), ignore.case = TRUE)]




res <- aggregate(goslim_goa_description ~ entrezgene_id, data = resultTable, paste)
membrane <- res[grep("membrane", res$goslim_goa_description, ignore.case = TRUE),]
intracellular <- res[grep("intracellular", res$goslim_goa_description, ignore.case = TRUE),]
secreted <- res[grep("secrete", res$goslim_goa_description, ignore.case = TRUE),]


soma.meta <- merge(meta, resultTable, by.x  = "EntrezGeneID", by.y = "entrezgene_id", all.x =T)

soma <- merge(soma, soma.meta, by.x = "Protein/Metabolite", by.y = "TargetFullName", all.x = TRUE)


names(soma)[12] <- "Annotation"
write.xlsx(soma, file = "./cluster_associations/results_9-3/soma_features.xlsx", row.names = FALSE)
}