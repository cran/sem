# last modified 7 March 2006 by J. Fox

specify.model <- function(file=""){
    ram <- scan(file=file, what=list(path="", par="", start=1, dump=""), sep=",", 
        strip.white=TRUE, comment.char="#", fill=TRUE) 
            # dump permits comma at line end
    ram <- cbind(ram$path, ram$par, ram$start)
    class(ram) <- "mod"
    ram
    }
    
print.mod <- function(x, ...){
    path <- x[,1]
    parameter <- x[,2]
    parameter[is.na(parameter)] <- "<fixed>"
    startvalue <- as.numeric(x[,3])
    startvalue[is.na(startvalue)] <- " "
    if (all(startvalue == " "))  print(data.frame(Path=path, Parameter=parameter),
        right=FALSE)
    else print(data.frame(Path=path, Parameter=parameter, StartValue=startvalue),
        right=FALSE)
    invisible(x)
    }
    
