# last modified 31 May 2005 by J. Fox

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
    startvalue <- as.numeric(x[,3])
    print(data.frame(Path=path, Parameter=parameter, StartValue=startvalue),
        right=FALSE)
    invisible(x)
    }
    
