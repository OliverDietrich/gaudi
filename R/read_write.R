#' Extract attributes from GFF3
#'
#' Extract values from GFF3 attributes and add columns to object
#'
#' @param object Data.frame formatted as GFF3
#' @param attribute_key String appended to columns to indicate they are extracted from 'attributes'
#'
#' @returns data.frame
#'
#' @export
extract_gff3_attributes <- function(object=NULL, attribute_key='attribute_') {

    stopifnot(
        !is.null(object),
        class(object) == 'data.frame'
    )

    # Checks
    if (!'attributes' %in% names(object)) {
        msg <- 'Column attributes not present in GFF3 file. Aborting.'
        stop(msg)
    }

    # Extract attributes
    x <- str_split(object[['attributes']], ';')

    # Format attributes
    ## For each list entry
    attr_names <- c()
    for (i in 1:length(x)) {
        v <- x[[i]]
        n <- length(v)
        v_values <- character(length = n)
        v_names <- character(length = n)
        # For each vector item
        for (j in 1:length(v)) {
            vj <- v[[j]]
            ss <- str_split(vj,'=')
            v_values[[j]] <- ss[[1]][[2]] # Righthand side becomes value
            v_names[[j]] <- ss[[1]][[1]] # Lefthand side becomes name
        }
        v <- v_values
        names(v) <- v_names
    
        # Store unique entry names
        ind <- which(!v_names %in% attr_names)
        attr_names <- c(attr_names, v_names[ind])
    
        # Return formatted vector
        x[[i]] <- v
    }

    # Check for duplicated column names
    if (any(attr_names %in% names(object))) {
        msg <- 'Conflicting attribute names. Appending attribute_ to each new column'
        warning(msg)
        stop('Not yet implemented.')
    }

    # Return attributes to object
    ## For each attribute
    for (i in attr_names) {
        key <- paste0(attribute_key,i)
        object[[key]] <- NA
        # For each row
        for (j in 1:length(x)) {
            if (i %in% names(x[[j]])) {
                object[j,key] <- x[[j]][[i]]
            }
        }
    }

    return(object)
}

#' Read GFF3 file
#'
#' Read file formatted as GFF3 (https://gmod.org/wiki/GFF3#GFF3_Annotation_Section)
#'
#' @param file Path to file
#' @param keep_fasta_sequence_as_attributes Boolean value to indicate whether to keep
#' appended FASTA sequences as attributes to the returned data.frame
#' @param verbose Boolean, whether to display file header and object structure
#'
#' @returns data.frame
#'
#' @export
read_gff3 <- function(file=NULL, keep_fasta_sequences_as_attributes=FALSE, verbose=FALSE) {

    stopifnot(
        !is.null(file),
        file.exists(file)
    )

    # Variables
    gff3_column_names <- c('seqid','source','type','start','end','score','strand','phase','attributes')

    # Read flat file
    flat <- readLines(file)

    ## Check for version
    version <- flat[[1]]
    if (version == '##gff-version 3') {
        if (verbose) {message('Reading gff-version 3.')}
    } else {
        msg <- paste0('Unknown header: ', print(version),'. Aborting.')
        stop(msg)
    }
    
    # Detect FASTA sequences
    fasta <- stringr::str_detect(flat, '^##FASTA')
    if (any(fasta)) {
        if (keep_fasta_sequences_as_attributes) {
            msg <- 'FASTA sequences present. Will be returned as attributes.'
            if (verbose) {warning(msg)}
        } else {
            msg <- 'FASTA sequences present. Will be removed.'
            if (verbose) {warning(msg)}
        }
        fa_start <- which(fasta)
        gff_start <- 1
        gff_stop <- fa_start-1
        flat <- flat[gff_start:gff_stop]
    }

    # Remove sequence regions
    seq_region <- stringr::str_which(flat, '^##sequence-region')
    flat <- flat[-seq_region]
    
    # Extract and print header
    header_region <- stringr::str_which(flat, '^#')
    header <- flat[header_region]
    header <- paste(header, collapse='\n')
    if (verbose) {
        cat('\n', header, '\n\n')
    }
    
    # Create object
    object <- read.table(text=flat, header = FALSE, sep = '\t', comment.char = '#', col.names = gff3_column_names)
    if (any(fasta) & keep_fasta_sequences_as_attributes) {
        attr(object, 'fasta') <- Biostrings::readDNAStringSet(file, seek.first.rec = TRUE)
    }

    # Extract attributes
    object <- extract_gff3_attributes(object)

    # View
    if (verbose) {
        str(object, max.level = 1)
    }

    return(object)
}

#' Write GFF3 file
#'
#' Write data.frame to GFF3 formatted file
#'
#' @param object Data.frame
#' @param file Character, file path
#'
#' @export
write_gff3 <- function() {

    
}

#' Extract FASTA headers
#'
#' @param input.fasta Character, input FASTA file
#'
#' @export
fasta_headers <- function(input.fasta=NULL) {

    stopifnot(
        !is.null(input.fasta)
    )

    # Check input
    is_file(input.fasta, suffix=c('.fasta','.fna'))

    # Read file
    object <- readLines(input.fasta)

    # Detect headers
    ind <- str_which(object,'>')
    object <- object[ind]
    
    # Format
    object <- str_remove(object, '>')
    object <- str_split(object,'_', simplify=TRUE) %>% as.data.frame()

    # Exit
    return(object)
}

#' Read assembly graph
#'
#' Read assembly graph from .gfa file. File specification is available at https://github.com/GFA-spec/GFA-spec.
#'
#' @param input.graph Character, path to input file
#'
#' @export
read_assembly_graph <- function(input.graph=NULL, remove.sequences=TRUE) {

    # Minimal check
    stopifnot(
        is_file(input.graph, suffix='gfa')
    )

    # Read file
    object <- readLines(input.graph)

    # List components
    components <- list(
        header = stringr::str_which(object,'^H'),
        segments = stringr::str_which(object,'^S'),
        links = stringr::str_which(object,'^L'),
        jumps = stringr::str_which(object,'^J'),
        containments = stringr::str_which(object,'^C'),
        paths = stringr::str_which(object,'^P'),
        walks = stringr::str_which(object,'^W'),
        comments = stringr::str_which(object,'^#')
    )

    # Filter components
    ind <- components %>% lapply(length) %>% unlist()
    ind <- ind[which(ind > 0)]
    components <- components[names(ind)]

    # Filter object
    for (i in names(components)) {
        components[[i]] <- object[components[[i]]]
        components[[i]] <- stringr::str_split(components[[i]], '\t', simplify=TRUE)
        components[[i]] <- as.data.frame(components[[i]])
    }

    # Format
    ## Header
    if ('header' %in% names(components)) {
        names(components[['header']])[[1]] <- 'RecordType'
        ### VersionNumber
        index <- unlist(lapply(lapply(components[['header']], stringr::str_detect, pattern='VN:Z:'), all))
        if (sum(index) == 1) {
            names(components[['header']])[index] <- 'VersionNumber'
            components[['header']][['VersionNumber']] <- stringr::str_remove(components[['header']][['VersionNumber']], 'VN:Z:')
        }        
        ### ProgramVersion
        if (all(stringr::str_detect(components[['header']][['V3']], 'sp:Z:'))) { # Requires position 3!
            names(components[['header']])[[3]] <- 'ProgramVersion'
            components[['header']][['ProgramVersion']] <- stringr::str_remove(components[['header']][['ProgramVersion']], 'sp:Z:')
        }
    }
    ## Segments
    if ('segments' %in% names(components)) {
        names(components[['segments']])[1:3] <- c('RecordType','Name','Sequence')
        ### K-mer counts
        index <- unlist(lapply(lapply(components[['segments']], stringr::str_detect, pattern='KC:i:'), all))
        if (sum(index) == 1) {
            names(components[['segments']])[index] <- 'KmerCount'
            components[['segments']][['KmerCount']] <- stringr::str_remove(components[['segments']][['KmerCount']], 'KC:i:')
            components[['segments']][['KmerCount']] <- as.numeric(components[['segments']][['KmerCount']])
        }
        ### Coverage
        index <- unlist(lapply(lapply(components[['segments']], stringr::str_detect, pattern='DP:f:'), all))
        if (sum(index) == 1) {
            names(components[['segments']])[index] <- 'Coverage'
            components[['segments']][['Coverage']] <- stringr::str_remove(components[['segments']][['Coverage']], 'DP:f:')
            components[['segments']][['Coverage']] <- as.numeric(components[['segments']][['Coverage']])
        }
        ### Length
        components[['segments']][['Length']] <- stringr::str_length(components[['segments']][['Sequence']])
        ### Order
        ind <- order(components[['segments']][['Length']], decreasing = TRUE)
        components[['segments']] <- components[['segments']][ind, ]
        ### Rownames
        row.names(components[['segments']]) <- 1:nrow(components[['segments']])
        components[['segments']][['Contig']] <- 1:nrow(components[['segments']])
        ### Remove sequence column (takes space, bad for printing)
        if (remove.sequences) {
            components[['segments']][['Sequence']] <- NULL
        }
    }
    ## Links
    if ('links' %in% names(components)) {
        names(components[['links']])[1:6] <- c('RecordType','From','FromOrient','To','ToOrient','Overlap')
        ### Map names to contigs
        lookup <- lookup <- setNames(components[['segments']][['Contig']], components[['segments']][['Name']])
        components[['links']][['FromContig']] <- lookup[components[['links']][['From']]]
        components[['links']][['ToContig']] <- lookup[components[['links']][['To']]]
    }
    ## Jumps
    ## Containments
    ## Paths
    ## Walks
    ## Comments

    # Create graph
    graph <- igraph::graph_from_data_frame(
        d = components$links[,c('FromContig','ToContig','From','To','Overlap','FromOrient','ToOrient')], 
        vertices = components$segments[,c('Contig','Name','Coverage','Length','KmerCount')], 
        directed = FALSE
    )

    # Transfer 
    components$segments$group <- components(graph)$membership
    components$graph <- graph
    
    # Exit
    return(components)
}