#' Load all intron counts from rds files outputed by nf-scsajr
#'
#' @param path2rds path to rds folder (nf-scsajr output)
#' @param sids vector of sample ifs
#'
#' @return combined dgCMatrix matrix
#' @export
#'
#' @examples
#' load_introns('path/to/rds',c('sid1','sid2'))
load_introns = function(path2rds,sids){
  ints = lapply(paste0(path2rds,'/',sids,'intron.rds'),readRDS)
  all_introns = unique(unlist(lapply(ints,rownames)))
  result = NULL
  for(i in seq_along(ints)){
    add_introns = setdiff(all_introns,rownames(ints[[i]]))
    tmp =  new("dgCMatrix",
               Dim = c(length(add_introns), ncol(ints[[i]])),
               Dimnames = list(add_introns, colnames(ints[[i]])),
               x = numeric(0),
               i = integer(0),
               p = integer(ncol(ints[[i]]) + 1))
    
    result = cbind(result,rbind(ints[[i]],tmp)[all_introns,])
  }
  rownames(result) = sub('-',':',sub(':-1:',':r:',sub(':1:',':f:',rownames(result))))
  result
}

#' Group intron by shared sites
#' 
#' The feature is splicing site*intron. Exlusion read count is sum of counts for all introns that share same site but not identical to the given.
#'
#' @param int_counts output of \code{\link{load_introns}}
#' @param cells dataframe with cell metadata, for example from cell_meta.rds created by nf-scsajr
#'
#' @return RangedSummarizedExperiment
#' @export
#'
#' @examples
#' ints = load_introns('path/to/rds',c('sid1','sid2'))
#' cells = readRDS('path/to/rds/cell_meta.rds')
#' group_introns(ints,)
group_introns = function(int_counts,cells){
  introns = data.frame(intron=rownames(int_counts))
  introns$left_site = paste0(visutils::splitSub(introns$intron,':',fixed = T,c(1,2,3)),':l')
  introns$right_site = paste0(visutils::splitSub(introns$intron,':',fixed = T,c(1,2,4)),':r')
  sites = c(introns$left_site,introns$right_site)
  sites = table(sites)
  sites = names(sites)[sites>1]
  
  sites_introns_l = introns[introns$left_site %in% sites,c('intron','left_site')]
  sites_introns_r = introns[introns$intron %in% sites,c('intron','left_site')]
  colnames(sites_introns_r) = colnames(sites_introns_l) = c('intron','site')
  sites_introns = rbind(sites_introns_l,sites_introns_r)
  
  rownames(sites_introns) = paste0(sites_introns$intron,':',visutils::splitSub(sites_introns$site,':',4))
  
  i = int_counts[sites_introns$intron,]
  e = t(visutils::calcColSums(t(i),sites_introns$site)[,sites_introns$site])
  e = e - i
  e = e[,rownames(cells)]
  i = i[,rownames(cells)]
  
  rownames(i) = rownames(e) = rownames(sites_introns)
  sitetypes = c('f:l'='d',
                'f:r'='a',
                'r:l'='a',
                'r:r'='d')
  sites_introns$sites = sitetypes[visutils::splitSub(rownames(sites_introns),':',c(2,5))]
  
  sites_introns$start = as.numeric(visutils::splitSub(rownames(sites_introns),':',c(3)))
  sites_introns$stop = as.numeric(visutils::splitSub(rownames(sites_introns),':',c(4)))
  sites_introns$chr_id = visutils::splitSub(rownames(sites_introns),':',c(1))
  sites_introns$strand = c('f'=1,'r'=1)[visutils::splitSub(rownames(sites_introns),':',c(2))]
  
  res = scsajr::make_summarized_experiment(list(seg=sites_introns,i=i,e=e),col_data = cells)
  res
}