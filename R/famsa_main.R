#' @export
famsa <- function(stringset = NULL, help = FALSE, verbose = FALSE,
                  very_verbose = FALSE, tree = FALSE, distance_matrix = FALSE, advanced_settings = "") {
  # stringset is an object of class 'XStringSet': DNAStringSet, RNAStringSet, or AAStringSet.

  args <- strsplit(advanced_settings, " ")[[1]]
  argv <- NULL
  if (help) {
    argv <- append(argv, c("-help"))
  } else {
    cl <- class(stringset)
    if (length(grep("StringSet", cl)) == 0) {
      stop("Input must be an object of class XStringSet: DNAStringSet, RNAStringSet, or AAStringSet!\nSee Biostrings: https://bioconductor.org/packages/release/bioc/html/Biostrings.html\n")
    }

    # Get names and sequence order for re-ordering post-algorithm.
    argv <- append(argv, c("famsa"))
    if (very_verbose) {
      argv <- append(argv, c("-vv"))
    } else if (verbose) {
      argv <- append(argv, c("-v"))
    }


    sequences_names <- names(stringset)
    if (!tree && !distance_matrix) {
      names(stringset) <- 1:length(stringset)
      tempOut <- tempfile(fileext = ".fasta")
    } else if (tree) {
      argv <- append(argv, c("-gt_export"))
      tempOut <- tempfile()
    } else if (distance_matrix) {
      argv <- append(argv, c("-dist_export"))
      tempOut <- tempfile(fileext = ".csv")
    } else {
      print("Error - both tree and distance matrix have been selected as the output!\nSet help=TRUE, to print instruction.")
    }

    if (length(args) > 0) {
      for (i in 1:length(args)) {
        if (is.logical(args[[i]])) {
          # arg is a flag.
          argv <- append(argv, as.character(names(args)[i]))
        } else {
          # arg is an option.
          argv <- append(argv, c(as.character(names(args)[i]), as.character(args[[i]])))
        }
      }
    }

    tempIn <- tempfile(fileext = ".fasta")
    Biostrings::writeXStringSet(stringset, filepath = tempIn, format = "fasta")

    # argv <- append(argv,c("famsa", "-gt_export", tempIn, tempOut))
    # argv <- append(argv,c("famsa", "-dist_export", tempIn, tempOut))
    # argv <- append(argv,c("famsa", "-dist_export", tempIn, tempOut))
    # argv <- append(argv,c("famsa", "-help", tempIn, tempOut))
    argv <- append(argv, c(tempIn, tempOut))
  }


  nargs <- as.integer(length(argv))
  # print(nargs)
  # print(argv)
  .C(.famsaCPP, nargs, as.character(argv))
  if (help) {
    return(NULL)
  }
  if (!tree && !distance_matrix) {
    ret <- switch(
            class(stringset),
            AAStringSet = Biostrings::readAAStringSet(tempOut, format = "fasta"),
            DNAStringSet = Biostrings::readDNAStringSet(tempOut, format = "fasta"),
            RNAStringSet = Biostrings::readRNAStringSet(tempOut, format = "fasta")
        )

    # Re-order and re-name output (alphabetical).
    ret <- ret[order(as.integer(names(ret)))]
    names(ret) <- sequences_names

    ret <- switch(
            class(stringset),
            AAStringSet = Biostrings::AAMultipleAlignment(ret),
            DNAStringSet = Biostrings::DNAMultipleAlignment(ret),
            RNAStringSet = Biostrings::RNAMultipleAlignment(ret)
        )

  } else if (tree) {
    ret <- phytools::read.newick(tempOut)
  } else {
    ret <- as.dist(
      data.matrix(
        read.csv(file=tempOut, row.names = 1, header = F, fill=T, col.names=append(sequences_names,NA,))
      )
    )
  }
  file.remove(tempIn, tempOut)

  return(ret)
}
