
# blast-parser -----------------------------------------------------------

#' @importFrom Biostrings BString
#' @importFrom XML xmlInternalTreeParse xmlRoot xpathApply xmlDoc xmlValue
NULL

#' Parse NCBI BLAST XML files.
#' 
#' @param x BLAST output in XML format (as file path or character string).
#' @param asText If \code{TRUE} the XML BLAST output is passed in as a string
#' instead of a file path.
#' @return A \code{\linkS4class{BlastReport}} object.
#' @export
#' @examples 
#' hla.xml <- system.file("extdata", "HLA-B-UTR.blast.xml", package = "blastr")
#' x <- blastReport(hla.xml)
#' 
#' ## look at the Header
#' getHeader(x)
#' 
#' ## look at the BLAST parameters
#' getParams(x)
#' 
#' ## look at the results of the first query
#' getQuery(x, 1)
#' 
#' ## or easier
#' x[[1]]
#' 
#' ## access the first two hits in the second query
#' getHit(getQuery(x, 2), 1:2)
#' 
#' ## for better readability use piping
#' \dontrun{
#' library("magrittr")
#' x %>% getQuery(2) %>% getHit(1:2)
#' }
#' 
#' ## or simple subsetting
#' x[[2]][1:2]
#' 
#' ## it's also possible to access the two best hits for all queries at once
#' getHit(x, 1:2)
#' 
#' ## access the hit ranges
#' getHitRange(x[[2]][1:2])
#' 
#' ## access evalues
#' getEvalue(x[[2]][1:2])
#' 
blastReport <- function(x, asText = FALSE) {
  doc <- xmlRoot(xmlInternalTreeParse(x, asText = asText))
  dbr <- parseBlastDatabaseReport(doc)
  query_elems <- xpathApply(doc, '/BlastOutput/BlastOutput_iterations/Iteration')
  new_BlastReport(
    header = dbr$header,
    parameters = dbr$params,
    queries = parseQueries(query_elems)
  )
}

parseBlastDatabaseReport  <- function(doc) {
  list(
    header = new_BlastHeader(
      version = xvalue(doc, '//BlastOutput/BlastOutput_version'),
      reference = xvalue(doc, '//BlastOutput/BlastOutput_reference'),
      database = xvalue(doc, '/BlastOutput/BlastOutput_db')
    ),
    params = new_BlastParameters(
      program = xvalue(doc, '/BlastOutput/BlastOutput_program'),
      matrix = xvalue(doc, '/BlastOutput/BlastOutput_param/Parameters/Parameters_matrix'),
      expect = xvalue(doc, '/BlastOutput/BlastOutput_param/Parameters/Parameters_expect', as = 'numeric'), 
      penalties = {
        open <- xvalue(doc, '/BlastOutput/BlastOutput_param/Parameters/Parameters_gap-open', as = 'integer')
        extend <- xvalue(doc, '/BlastOutput/BlastOutput_param/Parameters/Parameters_gap-extend', as = 'integer')
        c(open = open, extend = extend)
      },
      filter = xvalue(doc, '/BlastOutput/BlastOutput_param/Parameters/Parameters_filter'),
      sc_match = xvalue(doc, '/BlastOutput/BlastOutput_param/Parameters/Parameters_sc-match', as = 'integer'),
      sc_mismatch = xvalue(doc, '/BlastOutput/BlastOutput_param/Parameters/Parameters_sc-mismatch', as = 'integer'),
      num_sequences = xvalue(doc, '//Iteration[position()=1]/Iteration_stat/Statistics/Statistics_db-num'),
      num_letters = xvalue(doc, '//Iteration[position()=1]/Iteration_stat/Statistics/Statistics_db-len'),
      hsp_length = xvalue(doc, '//Iteration[position()=1]/Iteration_stat/Statistics/Statistics_hsp-len', as = 'numeric'),
      effective_space = xvalue(doc, '//Iteration[position()=1]/Iteration_stat/Statistics/Statistics_eff-space', as = 'numeric'), 
      ka_params = {
        k <- xvalue(doc, '//Iteration[position()=1]/Iteration_stat/Statistics/Statistics_kappa', as = 'numeric')
        lambda <- xvalue(doc, '//Iteration[position()=1]/Iteration_stat/Statistics/Statistics_lambda', as = 'numeric')
        h <- xvalue(doc, '//Iteration[position()=1]/Iteration_stat/Statistics/Statistics_entropy', as = 'numeric')
        c(k = k, lambda = lambda, h = h)
      }
    )
  )
}

parseQueries <- function(query_elems) {
  QueryList(
    lapply(query_elems, function(elem) {
      doc <- xmlDoc(elem)
      hit_elems <- xpathApply(doc, '/Iteration/Iteration_hits/Hit')
      query_env <- new.env(parent=emptyenv())
      query_env[['iter_num']] <- xvalue(doc, '/Iteration/Iteration_iter-num', as='integer')
      query_env[['query_id']] <- xvalue(doc, '/Iteration/Iteration_query-ID')
      query_env[['query_def']] <- xvalue(doc, '/Iteration/Iteration_query-def')
      query_env[['query_len']] <- xvalue(doc, '/Iteration/Iteration_query-len', as='integer')
      message <- xvalue(doc, '/Iteration/Iteration_message')
      new_Query(
        iter_num = query_env[['iter_num']] , query_id = query_env[['query_id']],
        query_def = query_env[['query_def']], query_len = query_env[['query_len']],
        hits = message %|na|% parseHits(hit_elems, query_env),
        query_env = query_env
      )
    })
  )
}

parseHits <- function(hit_elems, query_env) {
  HitList(
    lapply(hit_elems, function(elem) {
      doc <- xmlDoc(elem)
      hsp_elems <- xpathApply(doc, '/Hit/Hit_hsps/Hsp')
      hit_len <- xvalue(doc, '/Hit/Hit_len', as='integer')
      query_env[['hit_len']] <- hit_len
      new_Hit(
        hit_num = xvalue(doc, '/Hit/Hit_num', as='integer'),
        hit_def = Deflines(x=paste(xvalue(doc, '/Hit/Hit_id'),
                                   xvalue(doc, '/Hit/Hit_def'))
        ),
        hit_acc = xvalue(doc, '/Hit/Hit_accession'),
        hit_len = hit_len,
        hsps = parseHsps(hsp_elems, query_env),
        query_env = query_env
      )
    }),
    query_env = query_env
  )
}

parseHsps <- function(hsp_elems, query_env) {
  HspList(
    lapply(hsp_elems, function(elem) {
      doc <- xmlDoc(elem)
      new_Hsp(
        hsp_num = xvalue(doc, '/Hsp/Hsp_num', as='integer'),
        bit_score = xvalue(doc, '/Hsp/Hsp_bit-score', as='numeric'),
        score = xvalue(doc, "/Hsp/Hsp_score", as='numeric'),
        evalue = xvalue(doc, "/Hsp/Hsp_evalue", as='numeric'),
        query_from = xvalue(doc, "/Hsp/Hsp_query-from", as='integer'),
        query_to = xvalue(doc, "/Hsp/Hsp_query-to", as='integer'),
        hit_from = xvalue(doc, "/Hsp/Hsp_hit-from", as='integer'),
        hit_to = xvalue(doc, "/Hsp/Hsp_hit-to", as='integer'),
        query_frame = xvalue(doc, "/Hsp/Hsp_query-frame", as='integer'),
        hit_frame = xvalue(doc, "/Hsp/Hsp_hit-frame", as='integer'),
        identity = xvalue(doc, "/Hsp/Hsp_identity", as='integer'),
        positive = xvalue(doc, "/Hsp/Hsp_positive", as='integer'),
        gaps = xvalue(doc, "/Hsp/Hsp_gaps", as='integer'),
        align_len = xvalue(doc, "/Hsp/Hsp_align-len", as='integer'),
        qseq = BString(xvalue(doc, "/Hsp/Hsp_qseq")), 
        hseq = BString(xvalue(doc, "/Hsp/Hsp_hseq")),
        match = BString(xvalue(doc, "//Hsp_midline")),
        query_env = query_env
      )
    }),
    query_env = query_env
  )
}

#' Parse NCBI BLAST table files into \linkS4class{BlastTable} objects.
#' 
#' @param x BLAST output in tabular format (file path or character string)
#' @return A \code{\linkS4class{blastTable}} object.
#' @rdname blastTable
#' @export
#' @examples 
#' ##
blastTable <- function(x) {
  ## split multiple query results
  if (is.character(x) && !is.na(x[2]) && grepl("^# Query:.*$", x[2])) {
    start <- grep("^# Query:.*$", x) - 1
    end <- unique(c(grep("^# Query:.*$", x)[-1] - 2, length(x)))
    x <- Map(function(i) x[i], i = Map(seq, from = start, to = end))  
  }
  res <- list()
  for (blout in x) {
    if (is.character(blout) && file.exists(blout)) {
      cf <- count.fields(blout, sep = "\t", comment.char = "#")
      file_path <- file(blout, open = "r")
    } else {
      cf <- count.fields(textConnection(as.character(blout)),
                         sep = "\t", comment.char = "#")
      file_path <- textConnection(as.character(blout))
    }
    if (!all(cf == 12)) {
      stop(sprintf("Malformed blast table. %s columns in row %s.\n",
                   sQuote(cf[cf > 12]), sQuote(which(cf > 12))))
    }
    
    hasComment <- TRUE
    comStr <- NULL
    while (hasComment) {
      line <- readLines(file_path, n = 1)
      if (hasComment <- str_detect(line, "^#")) {
        comStr <- c(comStr, str_match(line, "[^(# *)].*"))
      }
    }
    pushBack(line, connection = file_path)  
    hit_table <- 
      read.table(file_path, header = FALSE, sep = "\t",
                 as.is = TRUE, nrows = length(cf),
                 col.names = c("qid", "sid", "pident",
                               "length", "mismatch", "gapopen",
                               "qstart", "qend", "sstart",
                               "send","evalue","bit_score"),
                 colClasses = c("character", "character", "numeric",
                                "integer", "integer", "integer",
                                "integer", "integer", "integer",
                                "integer", "numeric", "numeric"))
    
    ## parse subjectIds in 'hit_table'
    all_ids <- with(hit_table, strsplit(sid, "\\|"))
    gi  <- vapply(all_ids, '[', 2, FUN.VALUE = '')
    source_tag <- vapply(all_ids, '[', 3, FUN.VALUE = '')
    accn <- ifelse(grepl("pdb", source_tag),
                   paste0(sapply(all_ids, '[', 4), "_", sapply(all_ids, '[', 5)),
                   sapply(all_ids, '[', 4))
    if (all(is.na(accn))) {
      accn <- hit_table[["sid"]]
    }
    neg_log_evalue <- with(hit_table, -log(as.numeric(evalue)))
    neg_log_evalue[is.infinite(neg_log_evalue)] <- -log(1e-323)
    if (!is.null(comStr) && length(comStr) == 5) {
      program <- comStr[1]
      query <- str_split_fixed(comStr[2], "Query: ", 2)[, 2]
      database <- str_split_fixed(comStr[3], "Database: ", 2)[, 2]
    } else {
      program <- query <- database <- NA_character_
    }
    
    close(file_path)
    
    res <- c(res, new('BlastTable', program = program, database = database,
                      query = query,
                      bit_score = as.numeric(hit_table[["bit_score"]]),
                      evalue = as.numeric(hit_table[["evalue"]]),
                      mlog.evalue = neg_log_evalue,
                      geneid = gi,
                      accession = accn,
                      table = hit_table))
  }
  
  if (length(res) == 1) {
    res[[1]]
  } else res
}