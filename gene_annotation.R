library(rjson)
library(httr)
library(plyr)

PANTHER.update.tables <- function(table.path="./.PANTHERDB/support_table.Rda") {
  # get supported genomes
  GENOMES_URL <- "http://pantherdb.org/services/oai/pantherdb/supportedgenomes"
  response.genomes = GET(GENOMES_URL)
  supp.genomes <- fromJSON(content(response.genomes, "text"))$search$output$genomes$genome
  
  # create organism taxon-id lookup table
  tid.lookup = list()
  for (l in supp.genomes){
    tid.lookup[tolower(l$name)] = l$taxon_id
    tid.lookup[tolower(l$short_name)] = l$taxon_id
    tid.lookup[tolower(l$long_name)] = l$taxon_id
  }
  
  # get supported annotation databases
  DATASETS_URL= "http://pantherdb.org/services/oai/pantherdb/supportedannotdatasets"
  response.datasets = GET(DATASETS_URL)
  supp.datasets <- fromJSON(content(response.datasets, "text"))$search$annotation_data_sets$annotation_data_type
  
  # create dataset id lookup table
  did.lookup = list()
  for(l in supp.datasets){
    did.lookup[l$label] = l$id
  }
  
  # check lookup table path
  if (TABLE_PATH != table.path){
    TABLE_PATH=table.path
  }
  
  # save updated lookup table
  save(tid.lookup, did.lookup, file=TABLE_PATH)
}

TABLE_PATH <- "./.PANTHERDB/support_table.Rda"
CACHE_PATH <- "./.PANTHERDB/cache.json"
if (file.exists(TABLE_PATH)){
  PANTHER_TABLE <- load(TABLE_PATH)
} else {
  PANTHER.update.tables()
}

PANTHER.annotate <- function(gene.list, organism="human", anno.db="molecular_function", cache=T) {
  # convert list and organism name to PANTHER friendly format
  t.id <- tid.lookup[tolower(organism)]
  
  # generate query url
  ANNOTATION_URL = 'http://pantherdb.org/services/oai/pantherdb/geneinfo'
  query.string = '?geneInputList=%s&organism=%s'
  
  # process in chunks of 999 genes (upper limit of PANTHER)
  anno.df <- data.frame()
  s <- 1
  delta <- 999
  repeat{
    e <- s + delta
    chunk <- as.vector(gene.list[s:e])
    chunk <- chunk[!is.na(chunk)]
    
    gene.comma.list <- paste0(chunk, ",", collapse="")
    url = paste0(c(ANNOTATION_URL, sprintf(query.string, gene.comma.list, t.id)), "", collapse = "")
    
    # POST and await status code
    response.anno <- POST(url)
    stop_for_status(response.anno)
    
    # parse JSON and save if caching
    out.json <- fromJSON(content(response.anno, "text"))
    if (cache){
      jsonlite::write_json(out.json, CACHE_PATH)
    }
    
    # get annotation lookup table
    anno.lookup <- out.json$search$mapped_genes$gene
    
    # assign pathways to genes in gene list
    d.id <- did.lookup[tolower(anno.db)]
    anno.df.chunk <- data.frame()
    i <- 1
    repeat{
      # get list of all DB hits (length = # of DBs containing symbol or ID)
      anno.list <- anno.lookup[[i]]$annotation_type_list$annotation_data_type
      if (length(anno.list) < 1){
        anno.list.filt2 <- list(x = list(name=NA, id=NA))
        names(anno.list.filt2) <- gene.list[i]
      } else {
        # extract only hits from DB of interest
        anno.list.did <- lapply(anno.list, function(x) if (x$content==d.id) x$annotation_list)
        # filter out NA DB hits
        anno.list.filt <- Filter(Negate(is.null), anno.list.did)
        hits <- length(unlist(anno.list.filt)) / 2
        
        # extract annotations from DB (length = number of annotations)
        if (length(anno.list.filt) < 1){
          anno.list.filt2 <- list(x = list(name=NA, id=NA))
          names(anno.list.filt2) <- gene.list[i]
        } else {
          anno.list.filt2 <- lapply(seq(hits), function(i) anno.list.filt[[1]]$annotation[[i]])
          names(anno.list.filt2) <- rep(gene.list[i], hits)
        }
      }
      anno.df.chunk <- rbind(anno.df, ldply(anno.list.filt2, data.frame))
      i = i + 1
      if (i > length(gene.list)){
        break
      }
    }
    
    anno.df <- rbind(anno.df, anno.df.chunk)
    
    s <- e + 1
    if (s > length(gene.list)){
      break
    }
  }
  
  colnames(anno.df) <- c("Gene", "PATHWAY", "PATHWAY.ID")
  return(anno.df)
}
