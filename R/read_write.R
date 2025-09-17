#' Extract tag-value pairs by pattern
#'
#' @param x A character string containing tag-value pairs, sepearated by pattern
#' @param pattern A string separating tags from values (e.g. ID=myname)
#'
#' @importFrom stringr str_detect str_split
#'
#' @export
extract_tag_value_pairs <- function(x, pattern='=', verbose=FALSE) {

    # Minimal check
    stopifnot(
        is.character(x)
    )

    # Input
    index <- stringr::str_detect(x, pattern)
    if (!all(index)) {
        msg <- paste0('Missing pattern "', pattern, '" in some tag-value pairs. Removing...')
        if (verbose) warning(msg)
        #removed <- x[!index] # Originally thought to append as column, but it seems more clean to remove...
        x <- x[index]
    }

    # Main
    pairs <- str_split(x, pattern, simplify=TRUE)
    pairs <- setNames(pairs[,2], pairs[,1])

    # Exit
    return(pairs)
}

#' Extract attributes from GFF3
#'
#' Extract values from GFF3 attributes and add columns to object
#'
#' @param object Data.frame formatted as GFF3
#' @param attribute_key String appended to columns to indicate they are extracted from 'attributes'
#'
#' @returns data.frame
#'
#' @importFrom dplyr bind_rows
#' @importFrom stringr str_split
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
    x <- object[['attributes']]

    # Split by ';'
    x <- str_split(x, ';')

    # Split by '='
    x <- lapply(x, extract_tag_value_pairs)
    x <- bind_rows(x)

    # Adjust duplicated names
    ind <- names(x) %in% names(object)
    if (length(names)) {
        names(x)[ind] <- paste0('attribute_',names(x)[ind])
    }

    # Add back to object
    object <- cbind(object, x)
    
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
#' @importFrom readr read_table col_character col_integer col_number
#' @importFrom stringr str_detect str_which str_split
#'
#' @export
read_gff3 <- function(file, 
                      extract.attributes = TRUE,
                      keep.sequences = TRUE, 
                      verbose = FALSE
                     ) {

    # Minimal check
    stopifnot(
        file.exists(file)
    )

    # Variables
    gff3_columns <- list(
        'seqid' = col_character(),
        'source' = col_character(),
        'type' = col_character(),
        'start' = col_integer(),
        'end' = col_integer(),
        'score' = col_character(),
        'strand' = col_character(),
        'phase' = col_character(),
        'attributes' = col_character()
    )

    # Read flat file
    flat <- readLines(file)

    # Detect empty files
    if (!length(flat)) {
        msg <- paste('File',file,'is empty. Returning NULL!')
        warning(msg)
        return(NULL)
    }
    
    ## Check for version
    ind <- str_detect(flat, '##gff-version')
    if (sum(ind) >= 1) {
        version_header <- flat[ind][1]
        version <- str_split(version_header, pattern = '\\s+')[[1]][[2]]
        message('Reading gff version ', version)
    } else {
        msg <- paste0('Unknown header: Trying anyway...')
        warning(msg)
    }
    
    # Detect FASTA sequences
    fasta <- str_detect(flat, '^##FASTA')
    if (any(fasta)) {
        message('FASTA sequences present.')
        fa_start <- which(fasta)
    } else {
        fa_start <- Inf
    }

    # Remove sequence regions
    seq_region <- stringr::str_which(flat, '^##sequence-region')
    
    # Extract and print header
    header_region <- stringr::str_detect(flat, '^#')
    header <- flat[header_region]
    header <- paste(header, collapse='\n')
    if (verbose) cat('\n', header, '\n\n')
    
    # Create object
    n <- fa_start - sum(header_region)
    object <- read_tsv(file, col_types = readr::as.col_spec(gff3_columns), col_names = names(gff3_columns), comment = '#', n_max = n)
    object <- as.data.frame(object)

    # Extract GFF3 column 'attributes'
    object <- if (extract.attributes) extract_gff3_attributes(object) else object

    # Set object attributes
    if (keep.sequences & any(fasta)) {
        attr(object, 'sequence') <- Biostrings::readDNAStringSet(file, seek.first.rec = TRUE)
    }

    # Exit
    if (verbose) str(object, max.level = 1)
    return(object)
}

#' Enforce GFF3 formatting
#'
#' Format a data.frame according to GFF3 structure.
#' 
#' @param x Data.frame
#' @param replace.attributes Whether to replace the attributes column by gathering all non-standard columns
#' @param remove.duplicates Whether to remove rows that share the same seqid, type, start, and end to a row higher in the file.
#'
#' @importFrom dplyr group_by mutate
#'
#' @export
#'
format_gff3 <- function(x, replace.attributes = FALSE, remove.duplicates = FALSE) {

    # Check
    stopifnot(
        is.data.frame(x)
    )

    # Variables
    gff3_columns <- c('seqid','source','type','start','end','score','strand','phase','attributes')
    mandatory <- c('seqid','source','type','start','end')
    x_cols <- names(x)

    # Attributes
    index <- which(!x_cols %in% gff3_columns)
    if (length(index)) {
        attr <- x[, index]
        attr <- apply(attr, 1, function(y) setNames(as.list(y), names(attr)))
        attr <- lapply(attr, na.omit)
        attr <- lapply(attr, function(y) paste(names(y), y, sep='='))
        attr <- lapply(attr, paste, collapse=';')
        attr <- unlist(attr)
    } else {
        attr <- rep('', nrow(x))
    }

    # Main
    if (!all(mandatory %in% x_cols)) {
        msg <- paste0('The mandatory columns (', paste(mandatory), ') are not present.')
        stop(msg)
    }
    GFF <- data.frame(
        'seqid' = x[['seqid']],
        'source' = x[['source']],
        'type' = x[['type']],
        'start' = x[['start']],
        'end' = x[['end']],
        'score' = if ('score' %in% names(x)) x[['score']] else 0,
        'strand' = if ('score' %in% names(x)) x[['strand']] else '+',
        'phase' = if ('score' %in% names(x)) x[['phase']] else '.',
        'attributes' = if ('attributes' %in% x_cols & !replace.attributes) x[['attributes']] else attr
    )

    # Order IDs
    lvls <- unique(GFF$seqid)
    GFF$seqid <- factor(GFF$seqid, lvls)
        
    # Order source
    lvls <- unique(GFF$source)
    GFF$source <- factor(GFF$source, lvls)
        
    # Order data.frame
    GFF <- GFF[order(GFF$seqid, GFF$start, GFF$end, GFF$source, decreasing = c(FALSE, FALSE, TRUE, FALSE)), ]

    # Remove duplicated entries
    x <- group_by(GFF, seqid, type, start)
    x <- mutate(x, dup = duplicated(end))
    index <- !x$dup
    GFF <- if (remove.duplicates) GFF[index, ] else GFF

    # Exit
    return(GFF)
}

#' Write GFF3 file
#'
#' Write data.frame to GFF3 formatted file
#'
#' @param x Data.frame containing GFF3 formatted columns
#' @param filename File path
#' @param replace.attributes Whether to replace the attributes column 
#' with tag-value pairs gather from non-GFF columns.
#'
#' @export
write_gff3 <- function(x, filename, replace.attributes=FALSE) {

    # Get sequences
    fasta <- attr(x, 'sequence')

    # Format attributes
    x <- format_gff3(x, replace.attributes = replace.attributes)
    
    # Write header
    cat('##gff-version 3\n', file = filename)

    # Write main
    x <- split(x, x$seqid)
    names(x) <- paste0('##sequence-region ',names(x),'\n')
    for (i in names(x)) {
        cat(i, file = filename, append = TRUE)
        write.table(x[[i]], filename, sep = '\t', append = TRUE, row.names = FALSE, col.names=FALSE, quote=FALSE)
    }
    
    # Append FASTA
    if (length(fasta)) {
        cat('##FASTA\n', file = filename, append = TRUE)
        Biostrings::writeXStringSet(fasta, filename, append = TRUE, format = 'fasta')
    }
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
read_assembly_graph <- function(input.graph=NULL, remove.seqs = TRUE) {

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
        index <- unlist(lapply(lapply(components[['segments']], stringr::str_detect, pattern='dp:f:'), all))
        if (sum(index) == 1) {
            names(components[['segments']])[index] <- 'Coverage'
            components[['segments']][['Coverage']] <- stringr::str_remove(components[['segments']][['Coverage']], 'dp:f:')
            components[['segments']][['Coverage']] <- as.numeric(components[['segments']][['Coverage']])
        }
        ### Length
        index <- unlist(lapply(lapply(components[['segments']], stringr::str_detect, pattern='LN:i:'), all))
        if (sum(index) == 1) {
            names(components[['segments']])[index] <- 'Length'
            components[['segments']][['Length']] <- stringr::str_remove(components[['segments']][['Length']], 'LN:i:')
            components[['segments']][['Length']] <- as.numeric(components[['segments']][['Length']])
        } else {
            components[['segments']][['Length']] <- stringr::str_length(components[['segments']][['Sequence']])
        }
        ### Order
        ind <- order(components[['segments']][['Length']], decreasing = TRUE)
        components[['segments']] <- components[['segments']][ind, ]
        ### Rownames
        row.names(components[['segments']]) <- 1:nrow(components[['segments']])
        components[['segments']][['Contig']] <- 1:nrow(components[['segments']])

        ### Sequences
        if (remove.seqs) {
            components$sequence <- Biostrings::DNAStringSet(components[['segments']][['Sequence']])
            components[['segments']][['Sequence']] <- NULL
        }
    }
    
    ## Links
    if ('links' %in% names(components)) {
        names(components[['links']])[1:6] <- c('RecordType','From','FromOrient','To','ToOrient','Overlap')
        ### Map names to contigs
        lookup <- setNames(
            components[['segments']][['Contig']], 
            components[['segments']][['Name']]
        )
        components[['links']][['FromContig']] <- lookup[components[['links']][['From']]]
        components[['links']][['ToContig']] <- lookup[components[['links']][['To']]]        
    }
    ## Jumps
    ## Containments
    ## Paths
    ## Walks
    ## Comments

    # Create graph
    link_cols <- c('FromContig','ToContig','From','To','Overlap','FromOrient','ToOrient')
    link_cols <- link_cols[link_cols %in% names(components$links)]
    segment_cols <- c('Contig','Name','Coverage','Length','KmerCount')
    segment_cols <- segment_cols[segment_cols %in% names(components$segments)]
    graph <- igraph::graph_from_data_frame(
        d = components$links[, link_cols],
        vertices = components$segments[, segment_cols], 
        directed = FALSE
    )
    components$graph <- graph

    # Summarize outgoing links of segments
    lookup <- components$links
    lookup$N <- 1
    lookup <- dplyr::group_by(lookup, FromContig)
    lookup <- dplyr::summarize(lookup, ToContig = paste(ToContig, collapse = ','), LinkCount = sum(N))
    #lookup <- setNames(lookup$ToContig, lookup$FromContig)
    ind <- match(components$segments$Contig, lookup$FromContig)
    components$segments$ToContig <- lookup$ToContig[ind]
    components$segments$LinkCount <- lookup$LinkCount[ind]
    
    # Annotate segments
    components$segments$Group <- components(graph)$membership
    components$segments$N <- 1
    components$segments <- dplyr::mutate(dplyr::group_by(components$segments, Group), 
                                         Group_Length = sum(Length), 
                                         Group_Coverage = sum(Length * Coverage) / Group_Length,
                                         Group_Size = sum(N)
                                        )
    components$segments$Circular <- unlist(Map(grepl, components$segments$Contig, components$segments$ToContig))
    
    # Exit
    return(components)
}