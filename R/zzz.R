.onAttach<-function(...){
  # uses packageStartupMessage which can then be
  # surpressed

  version <- utils::packageVersion("mmds")
  built <- utils::packageDescription("mmds",fields="Built")

  hello <- paste("This is mmds ",version,"\nBuilt: ",built,sep="")
  packageStartupMessage(hello)
}
