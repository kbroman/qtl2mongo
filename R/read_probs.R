#' Read genotype probabilities from database
#'
#' Read genotype probabilities from database
#'
#' @param db Name of database
#' @param chr Single chromosome to read.
#' @param pos Optional length-2 vector specifying an interval of positions to read.
#' @param markers Optional vector of marker names to read. If
#'     \code{pos} is provided, \code{marker} is ignored. The markers
#'     must all be on the same chromosome.
#' @param url URL for mongo server
#'
#' @return Genotype probabilities, as an object of class \code{"calc_genoprob"}.
#'
#' @examples
#' \dontrun{
#' pr <- read_probs("iron_probs", chr="5", url="mongodb://localhost")
#'
#' pr <- read_probs("iron_probs", chr="5", pos=c(50, 75))
#'
#' pr <- read_probs("iron_probs", markers=c("D5Mit11", "c5.loc18", "c5.loc18.5"))
#' }
#'
#' @seealso \code{\link{write_probs}}
#'
#' @importFrom mongolite mongo
#' @export
read_probs <-
    function(db, chr=NULL, pos=NULL, markers=NULL, url="mongodb://localhost")
{

    if(is.null(chr) && is.null(markers))
        stop("Must specify either chr or markers")

    # read chromosome table
    m_chr <- mongolite::mongo("chr", db, url)
    chr_tab <- m_chr$find()

    if(!is.null(chr)) {
        if(length(chr) != 1) {
            chr <- chr[1]
            warning('Argument "chr" must be a single character string.')
        }

        if(!(chr %in% unlist(chr_tab$chr)))
            stop("Chromosome ", chr, " not in database")
    }
    is_x_chr <- as.logical(unlist(chr_tab$is_x_chr))
    names(is_x_chr) <- unlist(chr_tab$chr)

    # attributes
    m_meta <- mongolite::mongo("meta", db, url)
    attrib <- m_meta$find()
    attrib <- lapply(attrib, '[[', 1) # strip off list business

    alleles <- attrib$alleles
    crosstype <- attrib$crosstype
    alleleprobs <- attrib$alleleprobs
    ind <- attrib$ind[[1]]

    m_probs <- mongolite::mongo("probs", db, url)

    if(!is.null(pos)) {
        if(!("pos" %in% names( m_probs$find(fields='{"pos":1}', limit=1) )))
            stop("Database doesn't contain position information")

        if(!is.null(markers))
            warning('Argument "markers" ignored if "pos" is provided')

        if(length(pos) != 2 || pos[1] > pos[2])
            stop('Argument "pos" should have length 2 with pos[1] <= pos[2]')

        pr <- m_probs$find(query=paste0('{"chr":"',chr,'","pos":{"$gte":', pos[1],'},',
                                        '"pos":{"$lte":', pos[2], '}}'),
                           sort='{"marker_index":1}')

        if(is.null(pr) || length(pr$marker) == 0) {
            warning("No markers found in that interval")
            return(NULL)
        }
    }
    else if(is.null(markers)) {
        pr <- m_probs$find(query=paste0('{"chr":"',chr,'"}'),
                           sort='{"marker_index":1}')

        return(pr)
        if(is.null(pr) || length(pr$marker) == 0) {
            warning("No markers found on that chromosome")
            return(NULL)
        }
    }
    else {
        pr <- m_probs$find(paste0('{"marker":{"$in":[',
                                  paste0('"', markers, '"', collapse=','),
                                  ']}}'),
                           sort='{"marker_index":1}')

        chr <- unlist(pr$chr)
        if(length(unique(chr)) > 1)
            stop("markers on multiple chromosomes")
    }

    # turn into array and add attributes
    markers <- unlist(pr$marker)
    chr <- pr$chr[[1]]
    geno <- chr_tab$geno[[ which(unlist(chr_tab$chr) == chr) ]]

    pr <- list("1"=array(unlist(pr$probs), dim=c(length(ind), length(geno), length(markers))))
    dimnames(pr[[1]]) <- list(ind, geno, markers)
    names(pr) <- chr

    attr(pr, "crosstype") <- crosstype
    attr(pr, "is_x_chr") <- is_x_chr[chr]
    attr(pr, "alleles") <- alleles
    attr(pr, "alleleprobs") <- alleleprobs
    class(pr) <- c("calc_genoprob", "list")

    pr

}
