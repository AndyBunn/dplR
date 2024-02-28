read.fh <- function(fname, BC_correction = FALSE) {
  inp <- readLines(fname, ok=TRUE, warn=FALSE)
  
  ## Get start and end positions of headers and data blocks
  header.begin <- grep("^HEADER:$", inp)
  header.end <- grep("^DATA:(Tree|Single)$", inp)
  n <- length(header.end)
  if(n == 0) {
    stop('file has no data in "Tree" or "Single" formats')
  }
  ## For each data block in one of the supported formats, find the
  ## corresponding header block
  header.taken <- logical(length(header.begin))
  for (i in seq_len(n)) {
    n.preceding <- sum(header.begin < header.end[i] - 1)
    if (n.preceding == 0 || header.taken[n.preceding]) {
      stop("invalid file: HEADER and DATA don't match")
    } else {
      header.taken[n.preceding] <- TRUE
    }
  }
  if (!all(header.taken)) {
    warning("more HEADER blocks than DATA blocks in supported formats")
  }
  ## For each data block in one of the supported formats, find the
  ## following header block (or end of file)
  data.end <- numeric(n)
  for (i in seq_len(n-1)) {
    tmp <- header.begin[header.begin > header.end[i]]
    data.end[i] <- tmp[1]
  }
  tmp <- header.begin[header.begin > header.end[n]]
  if (length(tmp) > 0) {
    data.end[n] <- tmp[1]
  } else {
    data.end[n] <- length(inp) + 1
  }
  ## Forget headers that are not used by the data blocks
  header.begin <- header.begin[header.taken]
  
  ## Get essential metadata from headers
  keycodes <- character(n)
  lengths <- numeric(n)
  end.years <- numeric(n)
  start.years <- numeric(n)
  multipliers <- rep(1, n)
  divisors <- rep(100, n)
  site.code <- rep(NA_character_, n)
  tree.vec <- rep(NA_real_, n)
  core.vec <- rep(NA_real_, n)
  radius.vec <- rep(NA_real_, n)
  stemdisk.vec <- rep(NA_real_, n)
  pith.offset <- rep(NA_real_, n)
  for (i in seq_len(n)) {
    this.header <- inp[(header.begin[i]+1):(header.end[i]-1)]
    ## get keycode (= series id)
    this.keycode <- sub("KeyCode=", "", fixed=TRUE,
                        x=grep("^KeyCode=|Keycode=", this.header, value=TRUE))
    if (length(this.keycode) != 1) {
      string2 <- gettext('number of "KeyCode" lines is not 1',
                         domain="R-dplR")
      stop(gettextf("in series %s: ", as.character(i), domain="R-dplR"),
           string2, domain=NA)
    } else {
      keycodes[i] <- this.keycode
    }
    ## get length
    this.length <- sub("Length=", "", fixed=TRUE,
                         x=grep("^Length=", this.header, value=TRUE))
    if (length(this.length) != 1) {
      string2 <- gettext('number of "Length" lines is not 1',
                         domain="R-dplR")
      stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
           string2, domain=NA)
    } else {
      lengths[i] <- as.numeric(this.length)
    }
    ## get end year
    this.end.year <- sub("DateEnd=", "", fixed=TRUE,
                         x=grep("^DateEnd=", this.header, value=TRUE))
    if (length(this.end.year) != 1) {
      string2 <- gettext('number of "DateEnd" lines is not 1',
                         domain="R-dplR")
      stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
           string2, domain=NA)
    } else {
      end.years[i] <- as.numeric(this.end.year)
    }
    ## get start year
    this.start.year <- sub("DateBegin=", "", fixed=TRUE,
                           x=grep("^DateBegin=", this.header, value=TRUE))
    if (length(this.start.year) != 1) {
      if(length(this.end.year) == 1) {
        start.years[i] <- end.years[i]- lengths[i] +1
      } else {
        string2 <- gettext('number of "DateBegin" lines is not 1',
                         domain="R-dplR")
        stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
           string2, domain=NA) }
      } else {
        start.years[i] <- as.numeric(this.start.year)
        # check for BC - AD dates and add not existing year 0
        if (start.years[i] <0 && end.years[i]>0 && BC_correction == TRUE) {
          start.years[i] <- start.years[i]+1 }
      }
    ## correct BC dates to +1, because rwl-format uses non-existing year 0 and this is not necessary in the FH-format
    if (start.years[i]<0 && end.years[i]<0 && BC_correction == TRUE) {
      start.years[i] <- start.years[i] + 1
      end.years[i] <- end.years[i] + 1
    }
    ## get unit (by default, divide by 100)
    this.unit <- sub("Unit=", "", fixed=TRUE,
                     x=grep("^Unit=", this.header, value=TRUE))
    if (length(this.unit) == 1) {
      this.unit <- sub("mm", "", this.unit, fixed=TRUE)
      div.loc <- regexpr("/", this.unit, fixed=TRUE)
      if (div.loc > 0) {
        multipliers[i] <- as.numeric(substr(this.unit, 1, div.loc-1))
        divisors[i] <- as.numeric(substr(this.unit, div.loc+1,
                                         nchar(this.unit)))
      } else {
        multipliers[i] <- as.numeric(this.unit)
        divisors[i] <- 1
      }
      if (is.na(multipliers[i]) || is.na(divisors[i])) {
        string2 <- gettext('cannot interpret "Unit" line',
                           domain="R-dplR")
        stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
             string2, domain=NA)
      }
    } else if (length(this.unit) > 1) {
      string2 <- gettext('number of "Unit" lines is > 1',
                         domain="R-dplR")
      stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
           string2, domain=NA)
    }
    ## get site code
    this.site <- sub("SiteCode=", "", fixed=TRUE,
                     x=grep("^SiteCode=", this.header, value=TRUE))
    if (length(this.site) == 1) {
      site.code[i] <- this.site
    }
    ## get tree number
    this.tree <- sub("TreeNo=", "", fixed=TRUE,
                     x=grep("^TreeNo=", this.header, value=TRUE))
    if (length(this.tree) == 1) {
      tmp <- suppressWarnings(as.numeric(this.tree))
      if (identical(tmp, round(tmp))) {
        tree.vec[i] <- tmp
      }
    }
    ## get core number
    this.core <- sub("CoreNo=", "", fixed=TRUE,
                     x=grep("^CoreNo=", this.header, value=TRUE))
    if (length(this.core) == 1) {
      tmp <- suppressWarnings(as.numeric(this.core))
      if (identical(tmp, round(tmp))) {
        core.vec[i] <- tmp
      }
    }
    ## get radius number
    this.radius <- sub("RadiusNo=", "", fixed=TRUE,
                       x=grep("^RadiusNo=", this.header, value=TRUE))
    if (length(this.radius) == 1) {
      tmp <- suppressWarnings(as.numeric(this.radius))
      if (identical(tmp, round(tmp))) {
        radius.vec[i] <- tmp
      }
    }
    ## get stem disk number
    this.stemdisk <- sub("StemDiskNo=", "", fixed=TRUE,
                         x=grep("^StemDiskNo=", this.header, value=TRUE))
    if (length(this.stemdisk) == 1) {
      tmp <- suppressWarnings(as.numeric(this.stemdisk))
      if (identical(tmp, round(tmp))) {
        stemdisk.vec[i] <- tmp
      }
    }
    ## get pith offset (missing rings before start of series)
    this.missing <-
      sub("MissingRingsBefore=", "", fixed=TRUE,
          x=grep("^MissingRingsBefore=", this.header, value=TRUE))
    if (length(this.missing) == 1) {
      tmp <- suppressWarnings(as.numeric(this.missing))
      if (identical(tmp, round(tmp)) && tmp >= 0) {
        pith.offset[i] <- tmp + 1
      }
    }
  }
  
  ## calculate time span for data.frame
  min.year <- min(start.years)
  r.off <- min.year - 1
  max.year <- max(end.years)
  span <- min.year:max.year
  dendro.matrix <- matrix(NA, ncol = n, nrow = length(span))
  colnames(dendro.matrix) <- keycodes
  rownames(dendro.matrix) <- span
  ## get rid of comments (if any)
  strip.comment <- function(x) {
    strsplit(x, ";")[[1]][1]
  }
  for (i in seq_len(n)) { # loop through data blocks
    portion.start <- header.end[i] + 1
    portion.end <- data.end[i] - 1
    n.expected <- end.years[i] - start.years[i] + 1
    if (n.expected < 0) {
      stop(gettext("Series ", keycodes[i], " has an end date before the start date"))
    }
    if (portion.end < portion.start) {
      stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
           gettextf("too few values (expected %d, got %d)",
                    n.expected, 0, domain="R-dplR"), domain=NA)
    }
    portion <- inp[portion.start:portion.end]
    if (nchar(portion[1]) < 60 ||
        grepl(";", portion[1], fixed=TRUE)) { # data is in column format
      data <- as.numeric(vapply(portion, strip.comment, "foo"))
    } else { # data is in block format
      data <- numeric(length(portion) * 10)
      for (j in seq_along(portion)) {
        row.fwf <- substring(portion[j],
                             seq(from=1, by=6, length=10),
                             seq(from=6, by=6, length=10))
        row.numeric <- as.numeric(row.fwf)
        data[(j * 10 - 9):(j * 10)] <- row.numeric
      }
      ## Remove trailing zeros
      zeros <- which(data == 0)
      if (length(zeros) > 0) {
        nonzeros <- setdiff(zeros[1]:length(data), zeros)
        if (length(nonzeros) > 0) {
          zeros <- zeros[zeros > max(nonzeros)]
          if (length(zeros) > 0) {
            data <- data[-zeros]
          }
        } else {
          data <- data[-zeros]
        }
      }
    }
    data <- data * multipliers[i] / divisors[i]
    n.true <- length(data)
    if (n.true == n.expected) {
      ## write data into matrix
      dendro.matrix[(start.years[i]-r.off):(end.years[i]-r.off), i] <-
        data
    } else if (n.true < n.expected) {
      stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
           gettextf("too few values (expected %d, got %d)",
                    n.expected, n.true, domain="R-dplR"), domain=NA)
    } else if (all(is.na(data[(n.expected+1):n.true]))) {
      dendro.matrix[(start.years[i]-r.off):(end.years[i]-r.off), i] <-
        data[seq_len(n.expected)]
    } else {
      stop(gettextf("in series %s: ", keycodes[i], domain="R-dplR"),
           gettextf("too many values (expected %d, got %d)",
                    n.expected, n.true, domain="R-dplR"), domain=NA)
    }
  }
  cat(sprintf(ngettext(n,
                       "There is %d series\n",
                       "There are %d series\n",
                       domain="R-dplR"),
              n))
  start.years.char <- format(start.years, scientific=FALSE, trim=TRUE)
  end.years.char <- format(end.years, scientific=FALSE, trim=TRUE)
  seq.series.char <- format(seq_len(n), scientific=FALSE, trim=TRUE)
  cat(paste0(format(seq.series.char, width=5), "\t",
             format(keycodes, width=8), "\t",
             format(start.years.char, width=5, justify="right"), "\t",
             format(end.years.char, width=5, justify="right"), "\t",
             format(multipliers/divisors,
                    scientific=FALSE, drop0trailing=TRUE),"\n"), sep="")
  rwl <- as.data.frame(dendro.matrix) # return data.frame
  ## Create data.frame for site, tree, core, radius, stem disk IDs
  all.have.treeID <- !any(is.na(tree.vec))
  na.core <- is.na(core.vec)
  all.have.coreID <- !any(na.core)
  ## Try to find implicit core IDs (tree ID occurs once)
  if (all.have.treeID && !all.have.coreID) {
    foo <- table(tree.vec)
    measured.once <- as.numeric(names(foo)[foo == 1])
    core.vec[na.core & tree.vec %in% measured.once] <- 1
    all.have.coreID <- !any(is.na(core.vec))
  }
  ## Only include "ids" data.frame if all tree and core IDs are known
  if (all.have.treeID && all.have.coreID) {
    unique.sites <- unique(site.code)
    n.unique <- length(unique.sites)
    if (n.unique > 1) {
      site.vec <- match(site.code, unique.sites)
      tree.vec2 <- complex(n, NA_real_, NA_real_)
      total.dupl <- 0
      for (i in seq_len(n.unique)) {
        idx <- which(site.vec == i)
        ut <- unique(tree.vec[idx])
        for (this.tree in ut) {
          idx2 <- idx[tree.vec[idx] == this.tree]
          if (this.tree %in% tree.vec2) {
            tree.vec2[idx2] <- 1i * (total.dupl + 1)
            total.dupl <- total.dupl + 1
          } else {
            tree.vec2[idx2] <- this.tree
          }
        }
      }
      if (total.dupl > 0) {
        dont.change <- Im(tree.vec2) == 0
        existing <- unique(Re(tree.vec2[dont.change]))
        max.existing <- max(existing)
        if (max.existing < 1) {
          free.ids <- 1:total.dupl
        } else {
          free.ids <- which(!(1:max.existing %in% existing))
          free.ids <-
            c(free.ids,
              seq(from=max.existing+1, by=1,
                  length.out=max(0, total.dupl-length(free.ids))))
        }
        tree.vec2[!dont.change] <-
          free.ids[Im(tree.vec2[!dont.change])]
      }
      tree.vec2 <- Re(tree.vec2)
      adf <- data.frame(tree=tree.vec2, core=core.vec, site=site.vec,
                        row.names=keycodes)
    } else {
      adf <- data.frame(tree=tree.vec, core=core.vec, row.names=keycodes)
    }
    if (any(!is.na(radius.vec))) {
      adf <- cbind(adf, radius=radius.vec)
    }
    if (any(!is.na(stemdisk.vec))) {
      adf <- cbind(adf, stemDisk=stemdisk.vec)
    }
    attr(rwl, "ids") <- adf
    cat(gettext('Tree and core IDs were found. See attribute "ids".\n',
                domain="R-dplR"))
  }
  ## Include pith offset data.frame if some pith offsets are known
  na.po <- is.na(pith.offset)
  if (any(!na.po)) {
    attr(rwl, "po") <- data.frame(series=keycodes, pith.offset=pith.offset)
    if (any(na.po)) {
      cat(gettext('Pith offsets were found (some missing values). See attribute "po".\n',
                  domain="R-dplR"))
    } else {
      cat(gettext('Pith offsets were found (no missing values). See attribute "po".\n',
                  domain="R-dplR"))
    }
  }
  class(rwl) <- c("rwl", "data.frame")
  rwl
}
