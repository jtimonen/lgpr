library(lgpr)

# Parsing function
parse_df <- function(df, dict) {
  n <- dim(df)[1]
  DF <- matrix(0, n, 4)
  for (i in 1:n) {
    meas <- as.character(df$nam[i])
    parts <- strsplit(meas, split = "_")[[1]]
    grp <- substr(parts[2], 1, 4)
    if (grp == "Case") {
      group <- 1
    } else if (grp == "Cont") {
      group <- 0
    } else {
      print(grp)
      stop("Unknown group!")
    }
    if (length(parts) == 4) {
      age <- parts[3]
    } else {
      age <- parts[4]
    }
    DF[i, 1] <- gsub("[^0-9]", "", parts[1])
    DF[i, 2] <- which(dict == parts[2])
    DF[i, 3] <- gsub("[^0-9]", "", age)
    DF[i, 4] <- group
  }
  num <- as.numeric(DF)
  mat <- matrix(num, n, 4)
  mat <- cbind(mat, df$y, df$nf)
  colnames(mat) <- c("idx", "id", "age", "group", "y", "norm.fac")
  DF <- as.data.frame(mat)
  DF <- DF[order(DF$id, DF$age), ]
  DF$id <- as.factor(DF$id)
  DF$group <- as.factor(DF$group)
  return(DF)
}

# Data creation function
create_df <- function(idx, dat) {
  y <- dat$x[idx, ]
  nf <- dat$norm_fac
  gnam <- dat$gene_names[idx]
  ensg <- rownames(dat$x)[idx]
  cat(paste0("Getting data for ", ensg, "(", gnam, ")\n"))

  nam <- names(y)
  y <- as.numeric(y)
  df <- data.frame(nam, y, nf)

  # ID of individual is index in this sequence
  dict <- c(
    "Case1", "Case2", "Case3", "Case5", "Case9", "Case10", "Case11",
    "Control1", "Control2", "Control3", "Control5", "Control9",
    "Control10", "Control11"
  )
  DF <- parse_df(df, dict)
  DF <- data.frame(DF[, 2:6])

  # Add diseaseAges
  t_init <- c(12, 12, 18, 24, 18, 12, 18) # seroconversion time in months
  names(t_init) <- c(1:7)
  DF <- add_dis_age(DF, t_init)

  # Add categorical covariates
  sex <- rep(c(0, 1, 1, 0, 0, 0, 1), 2) # 0 = Female, 1 = Male
  loc <- rep(c(0, 0, 0, 0, 0, 1, 1), 2) # 0 = Finland, 1 = Estonia
  names(sex) <- c(1:14)
  names(loc) <- c(1:14)

  DF <- add_factor(DF, sex)
  DF <- add_factor(DF, loc)
  DF <- DF[, c(1, 2, 6, 7, 3, 4, 5)]
  norm_fac <- DF$norm.fac
  DF <- DF[, 1:6]
  out <- list(df = DF, norm_fac = norm_fac)
  return(out)
}
