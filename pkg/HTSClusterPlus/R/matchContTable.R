####################################################
## Function of permuting columns in contingency table 
## to compare two clustering results

#' Permute rows/columns in a contingency table comparing two data clusterings
#' 
#' Blah blah blah
#'
#' @param table_1 Partition from a data clustering
#' @param table_2 Partition from a data clustering
#'
#' @return Permuted table 
#' @export
#' 
#' @importFrom e1071 matchClasses
matchContTable <- function(table_1, table_2){
  tab <- table(table_1, table_2)
  ## Put larger clustering in rows if needed, nrow(tab) >= ncol(tab)
  transpose <- FALSE
  if(nrow(tab) < ncol(tab)) transpose <- TRUE;
  if(transpose==TRUE) tab <- t(tab);
  ## Order rows according to largest clusters
  tab <- tab[order(apply(tab,1,max), decreasing=TRUE),]
  ## Match best column with each row of tab
  ## Use unique indices as some columns might map to multiple rows 
  index <- matchClasses(tab, method=ifelse(nrow(tab)==ncol(tab), "exact", "rowmax"))
  tabord <- tab[,unique(index)]
  if(transpose==TRUE) tabord <- t(tabord)
  return(tabord)
}