# ---------------------------------------------------------- #
# -------------------- Global functions -------------------- #
# ---------------------------------------------------------- #

formatting_fun <- function(idx,or,ci.lb,ci.ub){
  if(!is.na(or[idx])){
    return(paste0(sprintf(or[idx], fmt="%.2f")," (",sprintf(ci.lb[idx], fmt="%.2f"),",",sprintf(ci.ub[idx], fmt="%.2f"),")"))
  } else {
    return("N/A")
  }
}