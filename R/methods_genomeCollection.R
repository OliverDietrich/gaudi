### Methods for the genomeCollection class

#-------------------------------------------------------------------------------
# constructor function for genomeCollection class
#-------------------------------------------------------------------------------

#' Create a \code{genomeCollection}
#'
#' Create a \code{genomeCollection} object from a data.frame with defined components.
#' 
#' @param ids Unique identifiers for each row in meta.data. Can be a column name or vector.
#' @param path Character of length 1, path to a 
#' @param meta.data A data.frame containing meta data associated to ids
#'
#' @return A \code{\link{genomeCollection}} object
#'
#' @rdname create_genomeCollection
#' @export
#'
genomeCollection <- function(ids, path, meta.data=NULL,
                             path.key = 'path', genome.key = 'genome', genes.key = 'genes',
                             cds.key = 'CDS', proteins.key = 'proteins'
                            ) {

    # Minimal check
    stopifnot(
        length(path) == 1,
        typeof(path) == 'character',
        length(ids) != 0,
        typeof(ids) == 'character'
    )

    # Check metadata
    if (is.data.frame(meta.data) & length(ids) == 1 & all(ids %in% names(meta.data))) {
        msg <- paste('Found column',ids,'in meta data. Reading IDs...')
        message(msg)
        id.vector <- meta.data[[ids]]
        meta.data[[ids]] <- NULL
        ids <- id.vector
    }
    if (any(duplicated(ids))) {
        msg <- "All entries must have unique IDs."
        stop(msg)
    }

    # Path
    if (!endsWith(path,'/')) {
        path <- paste0(path,'/')
    }
    if (!dir.exists(path)) {
        dir.create(path)
    }

    # Add meta data
    metadata <- data.frame(row.names = ids)
    meta.cols <- names(meta.data)
    ## Path
    if (path.key %in% meta.cols) {
        metadata[['path']] <- meta.data[[path.key]]
        meta.data[[path.key]] <- NULL
    } else {
        metadata[['path']] <- paste0(path,'assemblies/',ids,'/') 
        # TODO:
        # Error prone: pasting incomplete PATHs in the following steps...
        # Needs NA handling. 
    }
    ## Genome
    if (genome.key %in% meta.cols) {
        metadata[['genome']] <- meta.data[[genome.key]]
        meta.data[[genome.key]] <- NULL
    } else {
        metadata[['genome']] <- paste0(metadata$path,'genome.fna')
    }
    ## Genes
    if (genes.key %in% meta.cols) {
        metadata[['genes']] <- meta.data[[genes.key]]
        meta.data[[genes.key]] <- NULL
    } else {
        metadata[['genes']] <- paste0(metadata$path,'genes.gff')
    }
    ## CDS
    if (cds.key %in% meta.cols) {
        metadata[['CDS']] <- meta.data[[cds.key]]
        meta.data[[cds.key]] <- NULL
    } else {
        metadata[['CDS']] <- paste0(metadata$path,'cds.fna')
    }
    ## Proteins
    if (proteins.key %in% meta.cols) {
        metadata[['proteins']] <- meta.data[[proteins.key]]
        meta.data[[proteins.key]] <- NULL
    } else {
        metadata[['proteins']] <- paste0(metadata$path,'proteins.faa')
    }

    # Add remaining columns
    for (i in names(meta.data)) {
        metadata[[i]] <- meta.data[[i]]
    }

    # Create collection
    object <- methods::new("genomeCollection", 
                 index = ids,
                 path = path,
                 meta.data = metadata,
                 reads = list(),
                 assembly = list(), 
                 annotation = list(), 
                 distances = list(),
                 misc = list()
                )
    return(object)
}

#-------------------------------------------------------------------------------
# validity check for genomeCollection class object
#-------------------------------------------------------------------------------

.gc_validity <- function(object) {
    msg <- NULL

    # PATHs
    if (length(path(object)) > 1) {
      msg <- c(msg, "path(object) must be a single location.")
    }
    if (!dir.exists(path(object))) {
      msg <- c(msg, "path(object) must exist.")
    }

    # Meta.data
    if ( nrow(metadata(object)) != length(index(object)) ) {
        msg <- c(msg, "nrow of meta.data is not equal to length of index")
    }
    mandatory.cols <- c('path','genome','genes','CDS','proteins')
    cols.exist <- mandatory.cols %in% names(metadata(object))
    if (!all(cols.exist)) {
        cols.names <- paste(mandatory.cols[!cols.exist], collapse=', ')
        msg <- c(msg, paste(cols.names,"not part of meta.data."))
    }

    if (length(msg)) { return(msg) }
    return(TRUE)
}

methods::setValidity("genomeCollection", .gc_validity)

#-------------------------------------------------------------------------------
# show
#-------------------------------------------------------------------------------

.simple_slot_summary <- function(object, name) {

    # Check input
    data <- slot(object, name)
    if (is.null(data)) {
        msg <- paste0("Slot '",name,"' does not exist.")
        stop(msg)
    }
    
    # Reproduce str(object, max.level=0)
    slot.str <- paste(stringr::str_to_title(class(data)),'of',length(names(data)))

    # Print names
    slot.names <- names(data)
    if (length(slot.names) > 3) {
        slot.print <- paste(paste(head(slot.names,3), collapse=' '), "...")
    } else {
        slot.print <- paste(slot.names, collapse=' ')
    }

    paste(slot.str, slot.print, sep='   ')
}

.show.genomeCollection <- function(object) {

  # Samples
  if (length(index(object)) > 6) {
      index.print <- paste(paste(head(index(object),3), collapse=' '), "...", paste(tail(index(object),3), collapse=' '))
  } else {
      index.print <- paste(index(object), collapse=' ')
  }

  # Meta.data
  meta.cols <- names(object@meta.data)
  mandatory.cols <- c('path','genome','genes','CDS','proteins')
  if (length(meta.cols) > 6) {
      meta.print <- paste(paste(head(meta.cols,3), collapse=' '), "...", paste(tail(meta.cols,3), collapse=' '))
  } else {
      meta.print <- paste(meta.cols, collapse=' ')
  }
  meta.summary <- sapply(lapply(object@meta.data[, mandatory.cols], file.exists), sum)

  # Print
    cat(
        is(object),"\n","containing", length(index(object)), "entries:", index.print, "\n",
        "stored in", path(object), "\n",
        "\n",
        "meta.data: Total",ncol(object@meta.data),"columns, ",meta.print,"\n"
    )
    cat("Files exist?",paste(paste(names(meta.summary), meta.summary, sep=': '), collapse=', '))
    cat(
        "\n",
        "reads:",.simple_slot_summary(object, 'reads'),"\n",
        "assembly:",.simple_slot_summary(object, 'assembly'),"\n",
        "annotation: ",.simple_slot_summary(object, 'annotation'),"\n",
        "distances: ",.simple_slot_summary(object, 'distances'),"\n"
    )
}

#' @export 
setMethod("show", "genomeCollection", .show.genomeCollection)

#-------------------------------------------------------------------------------
# Accessors
#-------------------------------------------------------------------------------

#' Get the number of genomes in a genomeCollection
setMethod("nrow", "genomeCollection", function(x) nrow(x@meta.data))

          
#' Accessors for the 'index' element of an genomeCollection object.
#'
#' @description 
#' The \code{index} slot in an genomeCollection object holds
#' a character vector containing unique genome names.
#'
#' @author Oliver Dietrich
#' @export
#'
setMethod("index", "genomeCollection", function(x) x@index)

          
#' Accessors for the 'path' element of an genomeCollection object.
#'
#' @description 
#' The \code{path} slot in an genomeCollection object holds
#' a character vector containing file paths.
#'
#' @author Oliver Dietrich
#' @export
#'
setMethod("path", "genomeCollection", function(x) x@path)

          
#' Accessors for the 'meta.data' element of an genomeCollection object.
#'
#' @description 
#' The \code{meta.data} slot in an genomeCollection object holds
#' a data.frame containing unstructured information associated to each genome.
#'
#' @author Oliver Dietrich
#' @export
#'
setMethod("metadata", "genomeCollection", function(x) x@meta.data)

setMethod("metadata<-", "genomeCollection", function(x, value) {
  x@meta.data <- value
  validObject(x)
  x
})

          
#' Accessors to Reads objects contained in a genomeCollection object.
#'
setMethod("Reads", "genomeCollection", function(x, name=NULL) x@reads[[name]])

setMethod("Reads<-", "genomeCollection", function(x, name=NULL, value) {
  x@reads[[name]] <- value
  validObject(x)
  x
})

setMethod("ReadsNames", "genomeCollection", function(x) names(x@reads))


#' Accessors to Assembly objects contained in a genomeCollection object.
#'
setMethod("Assembly", "genomeCollection", function(x, name=NULL) x@assembly[[name]])

setMethod("Assembly<-", "genomeCollection", function(x, name=NULL, value) {
  x@assembly[[name]] <- value
  validObject(x)
  x
})

setMethod("Assemblies", "genomeCollection", function(x) names(x@assembly))

#' Accessors to Assembly objects contained in a genomeCollection object.
#' 
setMethod("Annotation", "genomeCollection", function(x, name=NULL) x@annotation[[name]])

setMethod("Annotation<-", "genomeCollection", function(x, name=NULL, value) {
  x@annotation[[name]] <- value
  validObject(x)
  x
})

setMethod("Annotations", "genomeCollection", function(x) names(x@annotation))

          
#-------------------------------------------------------------------------------
# Subsetting
#-------------------------------------------------------------------------------

setMethod("[", "genomeCollection", function(x, i, j, ..., drop=TRUE) {
  # ...
  validObject(x)
  x
})

setMethod("subset", "genomeCollection", function(x, i, j, ..., drop=TRUE) {
  # ...
  validObject(x)
  x
})

#-------------------------------------------------------------------------------
# Subobjects and metadata
#-------------------------------------------------------------------------------

#' Dollar-sign autocompletion
#'
#' @importFrom utils .DollarNames
#'
#' @export
.DollarNames.genomeCollection <- function(x, pattern = "") {
    grep(pattern, names(metadata(x)), value=TRUE)
}

#' Metadata
#'
#' Get and set genome meta data
#'
#' @param x A genomeCollection object
#' @param name Name of genome meta data
#'
setMethod("$", "genomeCollection", function(x, name) metadata(x)[[name]])

#' @param value A vector to add as genome meta data
#'
setReplaceMethod("$", "genomeCollection", function(x, name, value) {
  metadata(x)[[name]] <- value
  validObject(x)
  x
})

#' Subobjects and genome meta data
#'
#' The \code{[[} operator pulls either subobjects
#' (e.g. assembly, annotation, distances, ...) or
#' genome meta data from a genomeCollection object
#'
#' @param x A genomeCollection object
#' @param i Name of genome meta data
#' @param j
#'
setMethod("[[", c("genomeCollection","ANY","missing"), function(x, i, j, ...) {
  md <- slot(object = x, name = 'meta.data')
  if (rlang::is_missing(i)) {
    return(md)
  } else 
  if (is.null(i)) {
    return(NULL)
  } else 
  if (i %in% names(md)) {
    metadata(x)[[i, ...]]
  } else
  if (i %in% ReadsNames(x)) {
    Reads(x, i)
  } else
  if (i %in% Assemblies(x)) {
    Assembly(x, i)
  } else
  if (i %in% Annotations(x)) {
    Annotation(x, i)
  } else {
    return(NULL)
  }
  # Only works for some slots!
  # Check out https://github.com/satijalab/seurat-object/blob/main/R/seurat.R for implementation of subobjects ...
})

#' @param value A vector to add as subobject or genome meta data
#'
setReplaceMethod("[[", c("genomeCollection","ANY","missing"), function(x, i, j, ..., value) {
  metadata(x)[[i, ...]] <- value
  validObject(x)
  x
})