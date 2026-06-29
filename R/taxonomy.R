#' Infer taxonomic ranks from Genome Taxonomy Database (GTDB)
#'
#' GTDB-Tk v2.7.0 requires ~140G of external data which needs to be downloaded
#'    and extracted. This can be done automatically, or manually.
#'
#'    Automatic:
#'
#'        1. Run the command "download-db.sh" to automatically download and extract to:
#'            /projects/health_sciences/bms/microbiology/fineranlab/adaptable-phage-solutions/.envs/assembly/share/gtdbtk-2.7.0/db/
#'
#'    Manual:
#'
#'        1. Manually download the latest reference data:
#'            wget https://data.gtdb.aau.ecogenomic.org//releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz
#'
#'        2. Extract the archive to a target directory:
#'            tar -xvzf gtdbtk_r226_data.tar.gz -C "/path/to/target/db" --strip 1 > /dev/null
#'            rm gtdbtk_r226_data.tar.gz
#'
#'        3. Set the GTDBTK_DATA_PATH environment variable by running:
#'            conda env config vars set GTDBTK_DATA_PATH="/path/to/target/db"
#'
assign_bacterial_taxonomy_by_GTDB <- function() {
    
}

#' Infer taxonomic ranks for phage genomes from International Committee for Taxonomy of Viruses (ICTV) database
#'
assign_phage_taxonomy_by_ICTV <- function() {
    
}