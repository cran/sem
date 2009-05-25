# last modified 24 May 2009 by Adam Kramer (originally by J. Fox)
# last modified 25 May by J. Fox

path.diagram <- function(model, ...){
    UseMethod("path.diagram")
    }

path.diagram.sem <- function(model, out.file, min.rank=NULL, max.rank=NULL, same.rank=NULL,
	variables=model$var.names, parameters=rownames(model$ram), ignore.double=TRUE,
	edge.labels=c("names", "values", "both"), size=c(8,8), node.font=c("Helvetica", 14),
	edge.font=c("Helvetica", 10), rank.direction=c("LR", "TB"), digits=2, standardize=FALSE, ...){
	if(!missing(out.file)){
		handle <- file(out.file, "w")
		on.exit(close(handle))
	}
	else handle <- stdout()
	edge.labels <- match.arg(edge.labels)
	rank.direction <- match.arg(rank.direction)
	cat(file=handle, paste('digraph "', deparse(substitute(model)),
			'" {\n', sep=""))
	cat(file=handle, paste('  rankdir=', rank.direction, ';\n', sep=""))
	cat(file=handle, paste('  size="',size[1],',',size[2],'";\n', sep=""))
	cat(file=handle, paste('  node [fontname="', node.font[1],
			'" fontsize=', node.font[2], ' shape=box];\n', sep=""))
	cat(file=handle, paste('  edge [fontname="', edge.font[1],
			'" fontsize=', edge.font[2], '];\n', sep=""))
	cat(file=handle, '  center=1;\n')
	if (!is.null(min.rank)){
		min.rank <- paste('"', min.rank, '"', sep="")
		min.rank <- gsub(',', '" "', gsub(' ', '', min.rank))
		cat(file=handle, paste('  {rank=min ', min.rank, '}\n', sep=""))
	}
	if (!is.null(max.rank)){
		max.rank <- paste('"', max.rank, '"', sep="")
		max.rank <- gsub(',', '" "', gsub(' ', '', max.rank))
		cat(file=handle, paste('  {rank=max ', max.rank,'}\n', sep=""))
	}
	if (!is.null(same.rank)){
		for (s in 1:length(same.rank)){
			same <- paste('"', same.rank[s], '"', sep="")
			same <- gsub(',', '" "', gsub(' ', '', same))
			cat(file=handle, paste('  {rank=same ', same,'}\n', sep=""))
		}
	}
	latent <- variables[-(1:model$n)]
	for (lat in latent){
		cat(file=handle, paste('  "', lat, '" [shape=ellipse]\n', sep=""))
	}
	ram <- model$ram
	if (standardize) ram[, 5] <-  std.coef(model)[,2] 
	else ram[names(model$coeff), 5] <- model$coeff
	rownames(ram)[model$fixed] <- ram[model$fixed,5]
	values <- round(ram[,5], digits)
	heads <- ram[,1]
	to <- ram[,2]
	from <- ram[,3]
	labels <- if (edge.labels == "names") parameters else if (edge.labels == "values") values else paste(parameters,values,sep="=")
	direction <- ifelse((heads == 2), ' dir=both', "")
	for (par in 1:nrow(ram)){
		if((!ignore.double) || (heads[par] == 1))
			cat(file=handle, paste('  "',
					variables[from[par]], '" -> "', variables[to[par]],
					'" [label="', labels[par], '"', direction[par], '];\n', sep=""))
	}
	cat(file=handle, '}\n')
}
	
#path.diagram.sem <- function(model, out.file, min.rank=NULL, max.rank=NULL, same.rank=NULL,
#    variables=model$var.names, parameters=rownames(model$ram), ignore.double=TRUE,
#    edge.labels=c("names", "values"), size=c(8,8), node.font=c("Helvetica", 14),
#    edge.font=c("Helvetica", 10), rank.direction=c("LR", "TB"), digits=2, ...){
#    if(!missing(out.file)){
#        handle <- file(out.file, "w")
#        on.exit(close(handle))
#        }
#        else handle <- stdout()
#    edge.labels <- match.arg(edge.labels)
#    rank.direction <- match.arg(rank.direction)
#    cat(file=handle, paste('digraph "', deparse(substitute(model)),
#        '" {\n', sep=""))
#    cat(file=handle, paste('  rankdir=', rank.direction, ';\n', sep=""))
#    cat(file=handle, paste('  size="',size[1],',',size[2],'";\n', sep=""))
#    cat(file=handle, paste('  node [fontname="', node.font[1], 
#        '" fontsize=', node.font[2], ' shape=box];\n', sep=""))
#    cat(file=handle, paste('  edge [fontname="', edge.font[1],
#        '" fontsize=', edge.font[2], '];\n', sep=""))
#    cat(file=handle, '  center=1;\n')
#    if (!is.null(min.rank)){
#        min.rank <- paste('"', min.rank, '"', sep="")
#        min.rank <- gsub(',', '" "', gsub(' ', '', min.rank))
#        cat(file=handle, paste('  {rank=min ', min.rank, '}\n', sep=""))
#        }
#    if (!is.null(max.rank)){
#        max.rank <- paste('"', max.rank, '"', sep="")
#        max.rank <- gsub(',', '" "', gsub(' ', '', max.rank))
#        cat(file=handle, paste('  {rank=max ', max.rank,'}\n', sep=""))
#        }
#    if (!is.null(same.rank)){
#        for (s in 1:length(same.rank)){
#            same <- paste('"', same.rank[s], '"', sep="")
#            same <- gsub(',', '" "', gsub(' ', '', same))
#            cat(file=handle, paste('  {rank=same ', same,'}\n', sep=""))
#            }
#        }
#    latent <- variables[-(1:model$n)]
#    for (lat in latent){
#        cat(file=handle, paste('  "', lat, '" [shape=ellipse]\n', sep=""))
#        }
#    ram <- model$ram
#    ram[names(model$coeff), 5] <- model$coeff
#    rownames(ram)[model$fixed] <- ram[model$fixed,5]
#    values <- round(ram[,5], digits)
#    heads <- ram[,1]
#    to <- ram[,2]
#    from <- ram[,3]
#    labels <- if (edge.labels == "names") parameters else values
#    direction <- ifelse((heads == 2), ' dir=both', "")
#    for (par in 1:nrow(ram)){
#        if((!ignore.double) || (heads[par] == 1))
#        cat(file=handle, paste('  "', 
#            variables[from[par]], '" -> "', variables[to[par]], 
#            '" [label="', labels[par], '"', direction[par], '];\n', sep=""))
#        }
#    cat(file=handle, '}\n')
#    }
