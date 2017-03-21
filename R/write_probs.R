#' Write genotype probabilities to database
#'
#' Write genotype probabilities to database
#'
#' @param db Name of database
#' @param probs Genotype probabilities, as calculated by \code{\link[qtl2geno]{calc_genoprob}}.
#' @param map Optional map of marker positions (a list of vectors of positions).
#' @param url URL for mongo server
#' @param quiet If FALSE, print some information about progress.
#'
#' @return None. (See Details.)
#'
#' @details
#' The genotype probabilities are written to a single Mongo database in a pair of collections.
#'
#' The data are placed in a series of tables.
#' \itemize{

#' \item \code{probs} - the probabilities, with fields \code{marker},
#'     \code{chr}, \code{marker_index}, \code{probs}, and potentially
#'     \code{pos}.
#' \item \code{chr} - chromosome information, with fields \code{chr}, \code{is_x_chr}, \code{geno}
#' \item \code{meta} - Other information, with fields \code{ind},
#'     \code{alleles}, \code{crosstype}, \code{alleleprobs}
#' }
#'
#' @examples
#' \dontrun{
#' library(qtl2geno)
#' iron <- read_cross2(system.file("extdata", "iron.zip", package="qtl2geno"))
#' map <- insert_pseudomarkers(iron, map, error_prob=0.002)
#' pr <- calc_genoprob(iron, map)
#'
#' write_probs(db="iron_probs", pr, map)
#' }
#'
#' @seealso \code{\link{read_probs}}
#'
#' @importFrom mongolite mongo
#' @export
write_probs <-
    function(db, probs, map=NULL, url="mongodb://localhost", quiet=TRUE)
{
    # check inputs
    if(!is.null(map)) match_probs_map(probs, map)
    nind <- vapply(probs, nrow, 1)
    if(length(unique(nind)) != 1)
        stop("probs has different numbers of individuals on different chromosomes")
    ind <- lapply(probs, rownames)
    if(length(ind) > 1) { # more than one chromosome
        for(i in seq(along=ind)[-1]) {
            if(any(ind[[1]] != ind[[i]]))
                stop("probs has different row names on different chromosomes")
        }
    }

    # insert meta information
    m_meta <- mongolite::mongo("meta", db, url)
    if(m_meta$count() > 0) m_meta$drop() # drop what's there
    m_meta$insert( list(ind=ind,
                   alleles=attr(probs, "alleles"),
                   crosstype=attr(probs, "crosstype"),
                   alleleprobs=attr(probs, "alleleprobs")) )

    # is_x_chr attribute
    is_x_chr <- attr(probs, "is_x_chr")
    if(is.null(is_x_chr)) {
        warning("Missing is_x_chr attribute; assuming all are autosomes")
        is_x_chr <- rep(FALSE, length(probs))
    }

    # insert chr information
    m_chr <- mongolite::mongo("chr", db, url)
    if(m_chr$count() > 0) m_chr$drop() # drop what's there
    lapply(seq(along=probs), function(chr)
        m_chr$insert(list(chr=names(probs)[chr],
                          is_x_chr=is_x_chr[chr],
                          geno=colnames(probs[[chr]]),
                          chr_index=chr)))

    # insert chr information
    m_probs <- mongolite::mongo("probs", db, url)
    if(m_probs$count() > 0) m_probs$drop() # drop what's there
    for(chr in seq(along=probs)) {
        chr_nam <- names(probs)[chr]

        if(!quiet) message("Writing probs for chr ", chr_nam)

        mnames <- dimnames(probs[[chr]])[[3]]
        if(!is.null(map)) pos <- as.list(map[[chr]])
        else pos <- lapply(seq(along=mnames), function(a) NULL)

        lapply(seq(dim(probs[[chr]])[[3]]), function(i)
            m_probs$insert( list(marker=mnames[i],
                                 marker_index=i,
                                 chr=chr_nam,
                                 pos=pos[[i]],
                                 probs=probs[[chr]][,,i]) ))
    }

    # add indexes
    m_probs$index('{"marker":1}')
    m_probs$index('{"chr":1}')
    m_probs$index('{"marker_index":1}')

    invisible(NULL)
}


# check that probs and map conform
match_probs_map <-
    function(probs, map, map_name="map")
{
    nmar_probs <- vapply(probs, function(a) dim(a)[[3]], 1)

    if(!all(nmar_probs == vapply(map, length, 1)))
        stop("probs and ", map_name, " have different numbers of markers")
    if(!all(unlist(lapply(probs, function(a) dimnames(a)[[3]])) ==
            unlist(lapply(map, names))))
        stop("probs and ", map_name, " have different marker names")

    TRUE
}
