#' Start monitor ram usage in linux
#'
#' Create a log file in the working directory that then is elaborated by get_max_ram_used
#'
#' @param log_name Character string that will be given to the log file saved on the disk
#' @return Log file on the disk with the ram usage
#' @export
start_monitor_ram=function(log_name){
  #One core one chunk of pathways to process ----
  if(.Platform$OS.type == "unix") {
    message("Starting to monitor the RAM usage")
  }else{
    stop("Sorry but I can monitor the RAM usage only in Linux");
  }

  report_memory_path="tmp_memory_report_%s.txt"
  report_memory_path=sprintf(report_memory_path,log_name)

  command="free -h -s 5 >> %s"
  command=sprintf(command,report_memory_path)

  system(command=command,wait=F)
  return(report_memory_path)
}

#' End monitor ram usage in linux and get maximum value
#'
#' End monitor ram usage in linux and get maximum value
#'
#' @param log_name Character string that has been given to the log file saved on the disk with start_monitor_ram
#' @return Maximum ram used and collected in the log file
#' @importFrom data.table rbindlist
#' @export
get_max_ram_used=function(log_path){
  message("Determining the max RAM used")
  final_log_path=gsub("tmp_","final_",log_path)
  command1="cp %s %s"
  command1=sprintf(command1,log_path,final_log_path)

  command2="rm %s"
  command2=sprintf(command2,log_path)

  system(command=command1,wait=T)
  system(command=command2,wait=T)

  x=readLines(final_log_path)
  x=x[grep("Mem: ",x)]

  for(ri in 1:length(x)){
    y=x[ri]
    y=gsub(",.","",y)
    for(n_rep in seq(10,1)){
      spaces=paste(rep(" ",n_rep),collapse = "")
      y=gsub(spaces,",",y)
    }
    y=gsub("Mem:,,","",y)
    x[ri]=y
  }

  xl=strsplit(x,split=",",fixed = T)
  xdf=lapply(xl,function(x){
    df=as.data.frame(t(as.data.frame(x)))
    return(df)
  })
  df=as.data.frame(rbindlist(xdf))
  colnames(df)=c("total","used","free","shared","cache","available")
  used_ram=df$used
  used_ram=as.numeric(gsub("G","",used_ram))
  start_used=used_ram[1]
  max_used=max(used_ram)-start_used

  return(max_used)
}

#' Convert minutes time value in hours
#'
#' @param raw_min Double value
#' @return New double value formatted in hours
#' @export
convert_min2h=function(raw_min){
  if(raw_min>=60){
    h=floor(raw_min/60)
    mins=round(raw_min-(60*h))
    read_run_time=as.numeric(paste(h,mins,sep="."))
  }else{
    h=0
    mins=round(raw_min)
    if(mins>9){
      mins=paste(h,mins,sep=".")
    }else{
      mins=paste(h,".0",mins,sep="")
    }
    read_run_time=as.numeric(mins)
  }
  return(read_run_time)
}

