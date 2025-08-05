#' Calculate plaque forming units (PFU) from spot assay
#'
#' Calculate PFU from plaque counts observed for a dilution of phage stock. Requires the spot volume.
#'
#' @param dil.factor Dilution factor of the droplet
#' @param plaque.count Number of plaques counted from the droplet
#' @param droplet.volume Volume of the deposited droplet (mL)
#' 
#' @export
calc_pfu <- function(
    dil.factor = NULL,
    plaque.count = 10,
    droplet.volume = 2.5e-3
) {
    stopifnot(
        !is.null(dil.factor),
        is.numeric(dil.factor),
        is.numeric(plaque.count),
        is.numeric(droplet.volume)
    )
    
    # Combine data
    object <- data.frame(
        'dilution' = dil.factor,
        'droplet_volume' = droplet.volume,
        'plaque_count' = plaque.count
    )
    
    # Correct errors
    
    ## Dilution 0 (should be 1)
    ind <- which(object$dilution == 0)
    object$dilution[ind] <- 1
    
    ## NAs
    ind <- which(is.na(object$dilution))
    object$dilution[ind] <- 1
    
    ## Dilution 1 (implies no plaques)
    ind <- which(object$dilution == 1)
    object$plaque_count[ind] <- 0
    
    # Calculate PFU
    object <- object %>% mutate(PFU = (plaque_count/droplet_volume) * (1/dilution))
    
    # Exit 0
    return(object[['PFU']])
}

#' Calculate efficiency of plating/plaquing (EOP) from plaque forming units (PFU)
#'
#' @param test.pfu Numeric, vector of PFU units
#' @param iso.pfu Numeric, vector of PFU units
#'
#' @export
calc_eop <- function(
    test.pfu = NULL,
    iso.pfu = NULL
) {
    
    # Minimal check
    stopifnot(
        !is.null(test.pfu),
        !is.null(iso.pfu)
    )
    
    # Check iso.pfu
    if (any(iso.pfu == 0)) {
        msg <- 'Isolation PFU must not be 0, aborting...'
        stop(msg)
    }
    
    # Combine data
    object <- data.frame(
        'iso' = iso.pfu,
        'test' = test.pfu
    )
    
    # Calculate EOP
    object <- object %>% mutate(eop = log10(test/iso))
    
    # Correct technical artifacts
    ind <- which(object$test == 0)
    if (length(ind) > 0) {
        msg <- 'Some EOP values have been assumed from isolation PFU.'
        warning(msg)
        object$eop[ind] <- -log10(object$iso[ind] * 0.003)
    }
    
    # Exit 0
    return(object[['eop']])
}

#' Find NCBI database entries for a topic per year
#' 
#' @param db NCBI database to query
#' @param year Year to search entries for
#' @param term Search term
#' 
#' @export 
entrez_search_count_by_year <- function(db='nuccore', year=NULL, term=NULL){

    stopifnot(
        !is.null(year),
        !is.null(term)
    )
    
    # Run NCBI query
    query <- paste(term, "AND (", year, "[PDAT])")
    result <- rentrez::entrez_search(db=db, term=query, retmax=0)$count

    return(result)
}

#' Navigate to project root
#'
#' Navigate to the root of a project directory. The directory itself can be anywhere in a filesystem
#' but if e.g. Jupyter starts from a subdirectory this function navigates back to enable setting relative PATHs.
#' 
#' @param project.name Name of the project, needs to be part of the PATH
#' @param subdirs Character, Name of common subdirectories that can be recognized without project.name
#'
#' @export
navigate_to_project_root <- function(project.name = NULL,
                                     subdirs = c('/bin','/analysis','/data','/scripts')
                                    ) {

    # Get current directory
    wd <- getwd()

    # Define cases
    subdir_present <- endsWith(wd, subdirs)
    if (is.null(project.name)) {
        name_present <- FALSE
    } else {
        name_present <- stringr::str_detect(wd, project.name)
        dir_names <- str_split(wd,project.name)[[1]]
        subdir <- dir_names[2]
    }

    # Warning
    if (!is.null(project.name) & !name_present) {
        msg <- paste('Project name',project.name,'not found in',wd)
        warning(msg)
    }

    # Subdirs
    if (any(subdir_present)) {
        subdir <- subdirs[which(subdir_present)]
        wd <- stringr::str_remove(wd, subdir)
        setwd(wd)
        msg <- paste('Found subdirectory',subdir,'and navigated back to',wd)
        message(msg)
    } else if (name_present & subdir != '') {
        wd <- paste0(dir_names[1], project.name)
        setwd(wd)
        msg <- paste('Found project name and changed directory to',wd)
        message(msg)
    } else {
        msg <- paste('No subdirectory identified. No changes made to',wd)
        message(msg)
    }
}

#' Invoke a System Command
#'
#' system3 invokes the OS command specified by command.
#'
#' @param command the system command to be invoked, as a character string.
#' 
#' @export
system3 <- function(command) {
    
    # Redirect stderr
    command <- paste(command, '2>&1')

    # Call system
    stdout <- system(command, intern=TRUE)

    # Re-format output
    stdout <- paste(stdout, collapse='\n')

    # Print
    cat(stdout)
}

#' Activate conda environment
#' 
#' Add conda environment to PATH to call programs not available via R
#'
#' @param name Name of the conda environment
#' @param prefix Path to the environment
#' @param conda_home Path to the miniconda installation
#'
#' @export
activate_conda_env <- function(name=NULL, prefix=NULL, conda_home='~/miniconda3') {

    # Minimal check
    stopifnot(
        !is.null(conda_home) & is.character(conda_home) & length(conda_home) == 1
    )

    # Variables
    if (!endsWith(conda_home,'/')) {
        conda_home <- paste0(conda_home,'/')
    }
    conda_env_path <- paste0(conda_home,'envs/')
    all_envs <- list.files(conda_env_path)
    pwd <- getwd()
    
    # Checks
    if (is.null(name) & is.null(prefix)) {
        msg <- 'No name or prefix given. Please specify...'
        cmd <- paste0(conda_home,'bin/conda ','env list')
        system3(cmd)
        stop(msg)
    }
    if (!is.null(name) & !is.null(prefix)) {
        msg <- 'You can not specify a name AND prefix. Aborting...'
        stop(msg)
    }
    
    # By name
    if (!is.null(name)) {
        if (!endsWith(name,'/')) {
            name <- paste0(name,'/')
        }
        name_present <- name %in% all_envs
        if (!name_present) {
            msg <- paste0(name,' not present in ',conda_env_path,'. Treating it like a prefix...')
            warning(msg)
            conda_env <- paste0(pwd,'/',name,'bin')
        } else {
            conda_env <- paste0(conda_env_path,'/',name,'bin')
        }
    }

    # By prefix
    if (!is.null(prefix)) {
        if (!endsWith(prefix,'/')) {
            prefix <- paste0(prefix,'/')
        }
        conda_env <- paste0(pwd,'/',prefix,'bin')
    }

    if (!dir.exists(conda_env)) {
        msg <- paste('Directory',conda_env,'not found. Aborting...')
        stop(msg)
    }
    
    # Add to PATH
    PATH <- Sys.getenv('PATH')
    path_elements <- str_split(PATH, ':')[[1]]
    if (conda_env %in% path_elements) {
        msg <- paste(conda_env,'already part of PATH. Moving to top ...')
        warning(msg)
        path_elements <- c(conda_env, path_elements[which(path_elements != conda_env)])
        PATH <- paste(path_elements, collapse=':')
        Sys.setenv('PATH' = PATH)
    } else {
        msg <- paste('Adding',conda_env,'to PATH...')
        message(msg)
        PATH <- paste(conda_env, PATH, sep = ':')
        Sys.setenv('PATH' = PATH)
    }

    cat(PATH)
}