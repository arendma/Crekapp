expandcol <- function(df, cols) {
  #takes a data frame and expands it by truning columns in key value pairs
  if(length(cols)<dim(df)[2]) {
    oldcols <- df[,-cols]
  }
  resdf <- data.frame()
  for (i in cols) {
    newdf=data.frame(oldcols,x1=df[,i], x2=rep(colnames(df)[i], dim(df)[1]))
    resdf=rbind(resdf,newdf)
  }
  return(resdf)
}