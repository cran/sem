# with contributions by Adam Kramer and Michael Friendly (originally by J. Fox)
# last modified 2012-11-01 by J. Fox

path.diagram <- function(...){
    .Deprecated("pathDiagram", package="sem")
    pathDiagram(...)
}

pathDiagram <- function (model, ...)
{
    UseMethod("pathDiagram")
}

pathDiagram.semmod <- function(model, obs.variables, ...){    
    parse.path <- function(path) {                                           
        path.1 <- gsub("-", "", gsub(" ","", path))
        direction <- if (regexpr("<>", path.1) > 0) 2 
        else if (regexpr("<", path.1) > 0) -1
        else if (regexpr(">", path.1) > 0) 1
        else stop(paste("ill-formed path:", path))
        path.1 <- strsplit(path.1, "[<>]")[[1]]
        list(first=path.1[1], second=path.1[length(path.1)], direction=direction)
    }
    if ((!is.matrix(model)) | ncol(model) != 3) stop ("model argument must be a 3-column matrix")
    startvalues <- as.numeric(model[,3])
    par.names <- model[,2]
    n.paths <- length(par.names)
    heads <- from <- to <- rep(0, n.paths)
    for (p in 1:n.paths){
        path <- parse.path(model[p,1])
        heads[p] <- abs(path$direction)
        to[p] <- path$second
        from[p] <- path$first
        if (path$direction == -1) {
            to[p] <- path$first
            from[p] <- path$second
        }
    }
    ram <- matrix(0, p, 5)
    all.vars <- unique(c(to, from))
    latent.vars <- setdiff(all.vars, obs.variables)
    vars <- c(obs.variables, latent.vars)
    pars <- na.omit(unique(par.names))
    ram[,1] <- heads
    ram[,2] <- apply(outer(vars, to, "=="), 2, which)
    ram[,3] <- apply(outer(vars, from, "=="), 2, which)   
    par.nos <- apply(outer(pars, par.names, "=="), 2, which)
    if (length(par.nos) > 0)
        ram[,4] <- unlist(lapply(par.nos, function(x) if (length(x) == 0) 0 else x))
    ram[,5]<- startvalues
    colnames(ram) <- c("heads", "to", "from", "parameter", "start value")
    pars <- unique(na.omit(par.names))
    coeff <- rep(0, length(pars))
    names(coeff) <- pars
    fake.sem <- list(ram=ram, n=length(obs.variables), var.names=vars, coeff=coeff)
    class(fake.sem) <- "sem"
    pathDiagram(fake.sem, ...)
}


pathDiagram.sem <- function (model, file, min.rank = NULL, max.rank = NULL,
                             same.rank = NULL, variables = model$var.names, parameters, # = rownames(model$ram),
                             ignore.double = TRUE, edge.labels = c("names", "values",
                                                                   "both"), size = c(8, 8), node.font = c("Helvetica", 14),
                             edge.font = c("Helvetica", 10), rank.direction = c("LR",
                                                                                "TB"), digits = 2, standardize = FALSE, output.type=c("graphics", "dot"),
                             graphics.fmt="pdf", dot.options=NULL, ...)
{
    output.type <- match.arg(output.type)
    if (!missing(file)) {
        dot.file <- paste(file, ".dot", sep="")
        handle <- file(dot.file, "w")
        on.exit(close(handle))
        if (output.type == "graphics") graph.file <- paste(file, ".", graphics.fmt, sep="")
    }
    else handle <- stdout()
    edge.labels <- match.arg(edge.labels)
    rank.direction <- match.arg(rank.direction)
    cat(file = handle, paste("digraph \"", deparse(substitute(model)),
                             "\" {\n", sep = ""))
    cat(file = handle, paste("  rankdir=", rank.direction, ";\n",
                             sep = ""))
    cat(file = handle, paste("  size=\"", size[1], ",", size[2],
                             "\";\n", sep = ""))
    cat(file = handle, paste("  node [fontname=\"", node.font[1],
                             "\" fontsize=", node.font[2], " shape=box];\n", sep = ""))
    cat(file = handle, paste("  edge [fontname=\"", edge.font[1],
                             "\" fontsize=", edge.font[2], "];\n", sep = ""))
    cat(file = handle, "  center=1;\n")
    if (!is.null(min.rank)) {
        min.rank <- paste("\"", min.rank, "\"", sep = "")
        min.rank <- gsub(",", "\" \"", gsub(" ", "", min.rank))
        cat(file = handle, paste("  {rank=min ", min.rank, "}\n",
                                 sep = ""))
    }
    if (!is.null(max.rank)) {
        max.rank <- paste("\"", max.rank, "\"", sep = "")
        max.rank <- gsub(",", "\" \"", gsub(" ", "", max.rank))
        cat(file = handle, paste("  {rank=max ", max.rank, "}\n",
                                 sep = ""))
    }
    if (!is.null(same.rank)) {
        for (s in 1:length(same.rank)) {
            same <- paste("\"", same.rank[s], "\"", sep = "")
            same <- gsub(",", "\" \"", gsub(" ", "", same))
            cat(file = handle, paste("  {rank=same ", same, "}\n",
                                     sep = ""))
        }
    }
    latent <- variables[-(1:model$n)]
    for (lat in latent) {
        cat(file = handle, paste("  \"", lat, "\" [shape=ellipse]\n",
                                 sep = ""))
    }
    ram <- model$ram
    if (missing(parameters)){
        par.names <- names(coef(model))
        rownames(ram)[ram[, "parameter"] != 0] <- par.names[ram[, "parameter"]]
        rownames(ram)[ram[, "parameter"] == 0] <- ram[ram[, "parameter"] == 0, "start value"]
        parameters <- rownames(ram)
    }
    if (standardize)
        ram[, 5] <- std.coef(model)[, 2]
    else ram[names(model$coeff), 5] <- model$coeff
    #    rownames(ram)[model$fixed] <- ram[model$fixed, 5]
    values <- round(ram[, 5], digits)
    heads <- ram[, 1]
    to <- ram[, 2]
    from <- ram[, 3]
    labels <- if (edge.labels == "names")
        parameters
    else if (edge.labels == "values")
        values
    else paste(parameters, values, sep = "=")
    direction <- ifelse((heads == 2), " dir=both", "")
    for (par in 1:nrow(ram)) {
        if ((!ignore.double) || (heads[par] == 1))
            cat(file = handle, paste("  \"", variables[from[par]],
                                     "\" -> \"", variables[to[par]], "\" [label=\"",
                                     labels[par], "\"", direction[par], "];\n", sep = ""))
    }
    cat(file = handle, "}\n")
    if (output.type == "graphics" && !missing(file)){
        cmd <- paste("dot -T", graphics.fmt, " -o ", graph.file, " ", dot.options, " ", dot.file, sep="")
        cat("Running ", cmd, "\n")
        result <- try(system(cmd))
    }
    invisible(NULL)
}
