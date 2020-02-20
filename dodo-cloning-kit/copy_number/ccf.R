library(Hmisc)

'compute_ccf' <- function(vaf, tt, minor, purity, multiplicity = NULL)
{
  bs = vaf
  cg = tt
  major = tt-minor
  vlen = length(major)
  out = diff1 = rep(list(NA),vlen)
  for(i in 1:vlen) {
    majori = major[i]
    if (!is.na(majori)) {
      mm = 1:major[i]
      out[[i]] = bs[i]*(2*(1-purity[i])+cg[i]*purity[i])/mm/purity[i]
      diff1[[i]] = abs(out[[i]]-1)
    }
  }
  if (is.null(multiplicity)) {
    column.list = lapply(diff1,which.min)
    column.list[sapply(column.list, function(x)length(x)==0)] = NA
    column.list = unlist(column.list)
  } else {
    column.list = multiplicity 
  }
  rawccf = unlist(lapply(1:length(out),function(x)out[[x]][column.list[x]]))
  ccf = rawccf
  ccf[major==0] = NA
  ccf[ccf>1] = 1
  ccf[is.na(minor)] = NA
  return(invisible(list(ccf=ccf, rawccf=rawccf, out=out, multiplicity=column.list)))
}

'conf_ccf' <- function(alt, ref, tt, minor, purity, multiplicity, alpha=0.05)
{
  depth = alt+ref
  ci = binconf(alt, depth, alpha, method=c("exact"))
  lowerci = ci[,2]
  upperci = ci[,3]
  lower = compute_ccf(vaf=lowerci, tt, minor, purity, multiplicity)$ccf
  upper = compute_ccf(vaf=upperci, tt, minor, purity, multiplicity)$ccf
  return(invisible(list(lower=lower, upper=upper)))
}
