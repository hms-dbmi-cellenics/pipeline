generate_plotuuid <- function(sample_uuid, task_name, plot_idx) {

  if(sample_id != "") {
    return(paste(sample_uuid, task_name, plot_idx, sep="-"))
  }

  return(paste(task_name, plot_idx, sep="-"))

}