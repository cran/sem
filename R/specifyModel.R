# last modified 2011-08-07

specify.model <- function(...){
	.Deprecated("specifyModel", package="sem")
	specifyModel(...)
}

specifyModel <- function(file="", exog.variances=FALSE, endog.variances=TRUE, covs){
	add.variances <- function () {
		variables <- need.variance()
		nvars <- length(variables)
		if (nvars == 0) return(model)
		message("NOTE: adding ", nvars, " variances to the model")
		paths <- character(nvars)
		par.names <- character(nvars)
		for (i in 1:nvars) {
			paths[i] <- paste(variables[i], "<->", variables[i])
			par.names[i] <- paste("V[", variables[i], "]", sep = "")
		}
		model.2 <- cbind(c(model[, 1], paths), c(model[, 2], par.names), 
			c(model[, 3], rep(NA, length(paths))))
		class(model.2) <- "semmod"
		model.2
	}
	need.variance <- function () {
		all.vars <- classifyVariables(model)
		exo.vars <- all.vars$exogenous
		end.vars <- all.vars$endogenous
		variables <- logical(0)
		for (paths in model[, 1]) {
			vars <- strip.white(paths)
			vars <- sub("-*>", "->", sub("<-*", "<-", vars))
			vars <- sub("<->|<-", "->", vars)
			vars <- strsplit(vars, "->")[[1]]
			if (vars[1] != vars[2]) {
				for (a.variable in vars) {
					if (is.na(variables[a.variable])) variables[a.variable] <- TRUE
				}
			}
			else {
				variables[vars[1]] <- FALSE
			}
		}
		if (!exog.variances && length(exo.vars) > 0) variables[exo.vars] <- FALSE
		if (!endog.variances && length(end.vars) > 0) variables[end.vars] <- FALSE
		names(variables)[variables]
	}
    model <- scan(file=file, what=list(path="", par="", start=1, dump=""), sep=",", 
        strip.white=TRUE, comment.char="#", fill=TRUE) 
            # dump permits comma at line end
	model$par[model$par == ""] <- NA
    model <- cbind(model$path, model$par, model$start)
	if (!(missing(covs))){
		for (cov in covs){
			vars <- strsplit(cov, "[ ,]+")[[1]]
			nvar <- length(vars)
			for (i in 1:nvar){
				for (j in i:nvar){
					row <- c(paste(vars[i], "<->", vars[j]), 
						if (i == j) paste("V[", vars[i], "]", sep="") else paste("C[", vars[i], ",", vars[j], "]", sep=""),
						NA)
					if (row[2] %in% model[,2]) next
					model <- rbind(model, row)
				}
			}
		}
	}
	model <- removeRedundantPaths(model, warn=FALSE)
	add.variances()
    }
    
print.semmod <- function(x, ...){
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
	
classifyVariables <- function(model) {
	variables <- logical(0)
	for (paths in model[, 1]) {
		vars <- strip.white(paths)
		vars <- sub("-*>", "->", sub("<-*", "<-", vars))
		if (grepl("<->", vars)){
			vars <- strsplit(vars, "<->")[[1]]
			if (is.na(variables[vars[1]])) variables[vars[1]] <- FALSE
			if (is.na(variables[vars[2]])) variables[vars[2]] <- FALSE
		}
		else if (grepl("->", vars)){
			vars <- strsplit(vars, "->")[[1]]
			if (is.na(variables[vars[1]])) variables[vars[1]] <- FALSE
			variables[vars[2]] <- TRUE
		}
		else if (grepl("<-", vars)){
			vars <- strsplit(vars, "<-")[[1]]
			if (is.na(variables[vars[1]])) variables[vars[2]] <- FALSE
			variables[vars[1]] <- TRUE
		}
		else stop("incorrectly specified model")
	}
	list(endogenous=names(variables[variables]), exogenous=names(variables[!variables]))
}

strip.white <- function(x) gsub(' ', '', x)

removeRedundantPaths <- function(model, warn=TRUE){
	paths <- model[, 1]
	paths <- strip.white(paths)
	paths <- sub("-*>", "->", sub("<-*", "<-", paths))
	start <- regexpr("<->|<-|->", paths)
	end <- start + attr(start, "match.length") - 1
	arrows <- substr(paths, start, end)
	vars <- matrix(unlist(strsplit(paths, "<->|<-|->")), ncol=2, byrow=TRUE)
	for (i in 1:length(arrows)){
		if (arrows[i] == "<-"){
			arrows[i] <- "->"
			vars[i, ] <- vars[i, 2:1]
		}
	}
	vars <- cbind(vars, arrows)
	dupl.paths <- duplicated(vars)
	if (warn && any(dupl.paths)){
		warning("the following duplicated paths were removed: ", paste(model[dupl.paths, 1], collapse=", "))
	}
	model <- model[!dupl.paths, ]
	class(model) <- "semmod"
	model
}


