#Funciones para tratar los ids (conversiones entre ENSEMBL, SYMBOL y ENTREZID)

#1. Convertir IDs
convertIDs <- function(gene_ids, 
                       ori_id, 
                       new_id = "ENTREZID", #por defecto
                       specie = "human",
                       na_remove = T,
                       dup_remove = T) {
  
  #selección de organismo
  org <- if (specie == "human") {
    org.Hs.eg.db} else org.Mm.eg.db
  
  #validar ID orignial
  valid_ids <- c("ENSEMBL", "SYMBOL", "ENTREZID")
  
  if (!ori_id %in% valid_ids) {
    stop(paste("El identificador debe ser:", 
               paste(valid_ids)))
         }
  
  #identificadores
  columns <- unique(c(ori_id, new_id, 
                      "ENSEMBL", "ENTREZID", 
                      "SYMBOL")) #asegurar que se devuelven todas
  ids <- tryCatch({
    AnnotationDbi::select(org,
                          keys = gene_ids,
                          keytype = ori_id,
                          columns = columns)
  }, error = function(e) {
    stop(paste("Error en la conversión de los identificadores", e$message))
  })
  
  #eliminar NAs
  if (na_remove && ori_id %in% colnames(ids)) {
    ids <- ids[!is.na(ids[, new_id]), ]}
  
  #eliminar ids duplicados
  if (dup_remove && ori_id %in% colnames(ids)) {
    ids <- ids[!duplicated(ids[, new_id]), ]}
  
  return(ids)
}

#2. Detectar IDs
detectIDs <- function(gene_ids,
                      specie = "human") {
  
  #detectar patrón ensembl para humano o ratón
  ensembl <- if (specie == "human") {
    "^ENSG\\d{11}"
  } else {
    "^ENSMUSG\\d{11}"
  }
  if (mean(grepl(ensembl, gene_ids)) >= 0.8) { #permitir que haya NAs
    return("ENSEMBL")
  }
  #detectar patrón de entrezid (números)
  if (mean(grepl("^\\d+$", gene_ids))>= 0.8) {
    return("ENTREZID")
  }
  
  return("SYMBOL") #por defecto
}
  