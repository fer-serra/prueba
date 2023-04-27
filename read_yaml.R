if (!installed.packages()[, 1] == "yaml") {
  install.packages("yaml")
}
if (!installed.packages()[, 1] == "jsonlite") {
  install.packages("jsonlite")
}

require(jsonlite, warn.conflicts = F)
require(yaml, warn.conflicts = F)

yaml_read <- function (file) {
	out <- read_yaml (file=file)	
	out <- toJSON(out)
	out <- fromJSON(out)
	return(out)
}


