#' Minimap2
#'
#' Minimap2 alignment of FASTQ reads to FASTA reference.
#'
#' @param reference FASTA file
#' @param out.dir Path to output directory
#' @param short.1 FASTQ file of first short reads in each pair (R1)
#' @param short.2 FASTQ file of second short reads in each pair (R2)
#' @param unpaired FASTQ file of unpaired short reads (S)
#' @param long FASTQ or FASTA file of long reads
#' @param tool Name of alignment tool (minimap2 or bowtie2)
#' @param preset Preset for tuning the minimap 2 params. One of 'map-ont', 'map-pb', or 'map-ilcr' for long reads. 
#' Will default to 'sr' for short reads.
#' @param threads Number of threads
#'
#' @export
#'
map_reads_to_reference <- function(reference, out.dir,
                                   short.1=NULL, short.2=NULL, unpaired=NULL, long=NULL,
                                   tool = 'minimap2',
                                   preset = 'map-ont',
                                   file.prefix = NULL,
                                   remove.intermediates = TRUE,
                                   recompute = FALSE,
                                   threads = n_proc()
                                  ) {

    # Input
    dir.create(out.dir, recursive=TRUE, showWarnings=FALSE)
    read_1_exists <- if (is.null(short.1)) FALSE else file.exists(short.1)
    read_2_exists <- if (is.null(short.2)) FALSE else file.exists(short.2)
    paired <- if (read_1_exists & read_2_exists) TRUE else FALSE
    read_S_exists <- if (is.null(unpaired)) FALSE else file.exists(unpaired)
    read_L_exists <- if (is.null(long)) FALSE else file.exists(long)
    index <- c('paired'=paired, 'unpaired'=read_S_exists,'long'=read_L_exists)
    if (sum(index) == 0) {
        stop('No reads (pairs) found. Aborting...')
    } else 
    if (sum(index) == 1) {
        prefix <- names(index)[index]
        msg <- paste('Found',prefix,'reads. Mapping to',reference,'...')
        message(msg)
    } else if (sum(index) > 1) {
        prefix <- paste(names(index)[index], collapse=' and ')
        msg <- paste('Found', prefix, 'reads. Please specify a single type.')
        stop(msg)
    }

    # Variables
    preset <- if (prefix %in% c('paired','unpaired')) 'sr' else preset
    prefix <- if (!is.null(file.prefix)) file.prefix else prefix
    out.dir <- if(endsWith(out.dir,'/')) out.dir else paste0(out.dir,'/')
    base.name <- paste0(out.dir,prefix,'-',tool)
    out.sam <- paste0(base.name,'.','sam')
    out.bam <- paste0(base.name,'.','bam')
    out_sorted.bam <- paste0(base.name,'_sorted','.','bam')
    out.bai <- paste0(base.name,'_indexed','.','bai')
    out_mapped.bam <- paste0(base.name,'_mapped','.','bam')
    out.report <- paste0(base.name,'_report','.','txt')
    out.bcf <- paste0(base.name,'_variants','.','bcf')
    out.vcf <- paste0(base.name,'_variants','.','vcf')
    out.coverage <- paste0(base.name,'_coverage','.','tsv')

    # Alignment
    if (tool == 'minimap2') {
        fastq <- if (paired) paste(short.1, short.2) else if (unpaired) unpaired else if (long) long
        cmd <- paste('minimap2','-a','-x',preset,reference,fastq,'-o',out.sam)
    } else 
    if (tool == 'bowtie2') {
        fastq <- if (paired) paste('-1',short.1,'-2',short.2) else 
            if (unpaired) unpaired else 
            if (long) stop('Bowtie for long reads is not supported.')
        cmd <- paste('bowtie2','-x',out.dir,fastq)
        stop('No implemented yet.')
    } else {
        stop('Unknown alignment tool selected. Please specify minimap2 or bowtie2.')
    }
    if (!file.exists(out.sam) | recompute) {
        message('Running alignment...')
        system3(cmd)
    } else {
        message('Alignment file (SAM) already exists.')
        cat(cmd)
    }

    # Convert to BAM
    # samtools view -bT reference.fa test.sam > test.bam #(if header absent)
    # samtools view -bS test.sam > test.bam #(if header present)
    cmd <- paste('samtools','view','-bT',reference,out.sam,'>',out.bam)
    if (!file.exists(out.bam) | recompute) {
        message('Converting SAM to BAM...')
        system3(cmd)
    } else {
        message('Alignment file (BAM) already exists.')
        cat(cmd)
    }

    # Sort BAM
    # samtools sort test.bam test_sorted
    cmd <- paste('samtools','sort',out.bam,'-o',out_sorted.bam)
    if (!file.exists(out_sorted.bam) | recompute) {
        message('Sorting BAM file...')
        system3(cmd)
    } else {
        message('BAM file has already been sorted.')
        cat(cmd)
    }

    # Index BAM
    # samtools index test_sorted.bam test_sorted.bai
    cmd <- paste('samtools','index',out_sorted.bam,'-o',out.bai)
    if (!file.exists(out.bai) | recompute) {
        message('Indexing BAM file...')
        system3(cmd)
    } else {
        message('BAM file has already been indexed.')
        cat(cmd)
    }

    # Filter BAM
    # samtools view -h -F 4 blah.bam > blah_only_mapped.sam
    cmd <- paste('samtools','view','-h','-F','4',out_sorted.bam,'-o',out_mapped.bam)
    if (!file.exists(out_mapped.bam) | recompute) {
        system3(cmd)
    } else {
        message('BAM file has already been filtered.')
        cat(cmd)
    }    

    # Create report
    cmd <- paste('samtools','stats',out_sorted.bam,'--reference',reference,'--threads',threads,'>',out.report)
    if (!file.exists(out.report) | recompute) {
        system3(cmd, include.errors = FALSE)
    } else {
        message('Alignment report already exists.')
        cat(cmd)
    }
    
    # Call variants
    # bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
    cmd <- paste('bcftools','mpileup','-Ou','-f',reference,out_mapped.bam,'|','bcftools','call','-mv','-Ob','-o',out.bcf)
    if (!file.exists(out.bcf) | recompute) {
        message('Calling variants...')
        system3(cmd)
    } else {
        message('Variant calling file (VCF) present.')
        cat(cmd)
    }

    # Convert BCF -> VCF
    cmd <- paste('bcftools','view',out.bcf,'>',out.vcf)
    if (!file.exists(out.vcf) | recompute) {
        message('Converting BCF to VCF...')
        system3(cmd, include.errors = FALSE)
    } else {
        message('VCF file present.')
        cat(cmd)
    }

    # Compute coverage
    cmd <- paste('samtools','depth',out_sorted.bam,'-o',out.coverage)
    if (!file.exists(out.coverage) | recompute) {
        message('Computing coverage...')
        system3(cmd)
    } else {
        message('Coverage already computed.')
        cat(cmd)
    }

    # Cleanup
    remove.intermediates <- FALSE # During DEVELOPMENT
    if (remove.intermediates) {
        message('Removing (large) intermediates...')
        unlink(out.sam)
        unlink(out.bam)
        unlink(out_sorted.bam)
        unlink(out.bai)
    }

    # Exit
    message('Done.')
}

# MMseqs

# MUMer