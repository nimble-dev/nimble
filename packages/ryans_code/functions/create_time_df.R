create_time_df<-function(suite_output_time,nmcmcs){
  time=suite_output_time[1:nmcmcs]
  timing=data.frame(names(time),time)
  names(timing)=c('method','time')
  rownames(timing)=c()
  return(timing)
}