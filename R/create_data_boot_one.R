#'This function create one data for bootsrapping
#' @title create one bootsraping data
#' @param boot_count the count number of bootstrapping
#' @param Indata the data that you used
#' @export create_data_boot_one
create_data_boot_one=function(boot_count,Indata)
{
	#j.boot<-boot_count*7+j #print(j.boot) #set.seed(j.boot)
    data.boot<-as.data.frame(Indata[sample(1:nrow(Indata), replace=TRUE),])
  return(data.boot)
}
