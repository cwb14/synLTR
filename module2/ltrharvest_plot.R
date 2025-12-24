#!/usr/bin/env Rscript


Rscript ltrharvest_plot.R --scn Bdact_ltr.ltrtools.stitched2.scn --gff=Bdact_ltr.work/Bdact_ltr.ltrtools.internals.fa.rexdb-plant.dom.gff3 --tsv=Bdact_ltr_kmer2ltr_dedup -o Bdact_plots
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--tsv"), type="character", default=NULL,
              help="LTR-RT TSV (col1=name, col2=LTRlen1, col3=LTRlen2)"),
  make_option(c("--scn"), type="character", default=NULL,
              help="SCN file (needs col1 s(ret), col2 e(ret), col6 l(lLTR), col12 chrom)"),
  make_option(c("--gff"), type="character", default=NULL,
              help="Domain HMM GFF (TEsorter HMM gff; needs cols 4-5 and col9 attributes/ID)"),
  make_option(c("-o","--output"), type="character", default="plot",
              help="Output prefix [default: %default]"),
  make_option(c("--max_individual"), type="integer", default=10000,
              help="Max elements to draw in an individual-elements plot per lineage [default: %default]")
)

parser <- OptionParser(option_list=option_list, usage="Rscript plot.R --tsv ltrrt.tsv --scn file.scn --gff hmm.gff -o outprefix")
opt <- parse_args(parser)

if (is.null(opt$tsv) || is.null(opt$scn) || is.null(opt$gff)) {
  print_help(parser)
  stop("Mandatory: --tsv, --scn, --gff")
}

# -----------------------------
# Helpers: parsing
# -----------------------------
trim <- function(x) gsub("^\\s+|\\s+$","", x)

# Column1 format: chrom:start-end#class/superfamily/family (example: Bdact_chr1:1689985-1695563#LTR/Gypsy/Reina)
parse_tsv_name <- function(x){
  # returns list(key, chrom, s, e, class, superfamily, family)
  parts <- strsplit(x, "#", fixed=TRUE)[[1]]
  key <- parts[1]
  tax <- if (length(parts) >= 2) parts[2] else ""
  tax_parts <- strsplit(tax, "/", fixed=TRUE)[[1]]
  class <- if (length(tax_parts) >= 1) tax_parts[1] else NA_character_
  superfamily <- if (length(tax_parts) >= 2) tax_parts[2] else NA_character_
  family <- if (length(tax_parts) >= 3) tax_parts[3] else NA_character_

  # parse key: chrom:start-end
  m <- regexec("^(.+):([0-9]+)-([0-9]+)$", key)
  r <- regmatches(key, m)[[1]]
  if (length(r) != 4) {
    return(list(key=key, chrom=NA, s=NA_integer_, e=NA_integer_,
                class=class, superfamily=superfamily, family=family))
  }
  list(
    key = key,
    chrom = r[2],
    s = as.integer(r[3]),
    e = as.integer(r[4]),
    class = class,
    superfamily = superfamily,
    family = family
  )
}

# Extract from GFF attributes column:
# - TE key: before "|" (Bdact_chr1:1689985-1695563)
# - gene/domain: from gene=... if present, else from Name=... (e.g., Reina-RT), else last token after ":" in ID
parse_gff_attrs <- function(attr){
  attr <- as.character(attr)
  # TE key tends to be in ID=... and before a '|' :
  # ID=Bdact_chr1:1689985-1695563|Class_I/LTR/...:Ty3-RT;...
  te_key <- NA_character_
  m <- regexec("ID=([^;]+)", attr)
  r <- regmatches(attr, m)[[1]]
  if (length(r) == 2) {
    idval <- r[2]
    te_key <- strsplit(idval, "\\|")[[1]][1]
  }

  # gene=GAG; etc
  gene <- NA_character_
  mg <- regexec("gene=([^;]+)", attr)
  rg <- regmatches(attr, mg)[[1]]
  if (length(rg) == 2) gene <- rg[2]

  # fallback: Name=Reina-RT
  if (is.na(gene)) {
    mn <- regexec("Name=([^;]+)", attr)
    rn <- regmatches(attr, mn)[[1]]
    if (length(rn) == 2) {
      nm <- rn[2]
      # pull last hyphen token if looks like Reina-RT, Ty3-RT, etc
      if (grepl("-", nm)) {
        gene <- tail(strsplit(nm, "-", fixed=TRUE)[[1]], 1)
      } else {
        gene <- nm
      }
    }
  }

  # Normalize gene labels to your plotting set
  # Common ones: GAG, PROT, RT, RH, INT, CH/CHD
  gene_norm <- toupper(gene)
  gene_norm <- gsub("^TY3-","", gene_norm)
  gene_norm <- gsub("^TY1-","", gene_norm)
  if (gene_norm %in% c("CH")) gene_norm <- "CH"
  if (gene_norm %in% c("CHD")) gene_norm <- "CH"
  if (gene_norm %in% c("IN")) gene_norm <- "INT"
  if (gene_norm %in% c("INTEGRASE")) gene_norm <- "INT"
  if (!gene_norm %in% c("GAG","PROT","RT","RH","INT","CH")) {
    # keep as-is (some lineages may have extra labels); will still plot with grey if not in color map
    gene_norm <- gene_norm
  }

  list(te_key=te_key, gene=gene_norm)
}

# -----------------------------
# Read inputs
# -----------------------------
# (1) TSV
tsv <- read.table(opt$tsv, sep="\t", header=FALSE, quote="", comment.char="", stringsAsFactors=FALSE)
if (ncol(tsv) < 3) stop("TSV must have at least 3 columns (name, LTRlen1, LTRlen2).")
colnames(tsv)[1:3] <- c("name","ltr1","ltr2")

tsv_parsed <- lapply(tsv$name, parse_tsv_name)
tsv$key <- vapply(tsv_parsed, `[[`, character(1), "key")
tsv$chrom <- vapply(tsv_parsed, `[[`, character(1), "chrom")
tsv$s <- vapply(tsv_parsed, `[[`, integer(1), "s")
tsv$e <- vapply(tsv_parsed, `[[`, integer(1), "e")
tsv$class <- vapply(tsv_parsed, `[[`, character(1), "class")
tsv$superfamily <- vapply(tsv_parsed, `[[`, character(1), "superfamily")
tsv$family <- vapply(tsv_parsed, `[[`, character(1), "family")

tsv$ltr_len <- pmax(as.integer(tsv$ltr1), as.integer(tsv$ltr2), na.rm=TRUE)
tsv$full_len <- as.integer(tsv$e - tsv$s + 1)

# lineage label for grouping in plots (Gypsy/Reina)
tsv$lineage <- ifelse(!is.na(tsv$superfamily) & !is.na(tsv$family),
                      paste0(tsv$superfamily, "/", tsv$family),
                      ifelse(!is.na(tsv$family), tsv$family, "Unknown"))

# (2) SCN
# (2) SCN (robust parser: tolerate ragged lines)
scn_lines <- readLines(opt$scn, warn = FALSE)

# drop empty/comment-ish lines
scn_lines <- scn_lines[trimws(scn_lines) != ""]
scn_lines <- scn_lines[!grepl("^\\s*#", scn_lines)]

# split on whitespace
scn_split <- strsplit(scn_lines, "\\s+")

# keep only lines with at least 12 fields
keep <- vapply(scn_split, length, integer(1)) >= 12L
if (!all(keep)) {
  bad_n <- sum(!keep)
  warning(sprintf("SCN: dropped %d malformed line(s) with <12 columns", bad_n))
}
scn_split <- scn_split[keep]

# build data.frame from the needed columns
scn <- data.frame(
  sret  = as.integer(vapply(scn_split, `[`, character(1), 1)),
  eret  = as.integer(vapply(scn_split, `[`, character(1), 2)),
  lLTR  = as.integer(vapply(scn_split, `[`, character(1), 6)),
  chrom =           vapply(scn_split, `[`, character(1), 12),
  stringsAsFactors = FALSE
)

# build key and shift map
scn$key <- paste0(scn$chrom, ":", scn$sret, "-", scn$eret)
scn_shift <- scn$lLTR
names(scn_shift) <- scn$key

# (3) GFF (domain hits)
gff <- read.table(opt$gff, sep="\t", header=FALSE, quote="", comment.char="#", stringsAsFactors=FALSE)
if (ncol(gff) < 9) stop("GFF must have 9 columns.")
colnames(gff)[1:9] <- c("seqid","source","type","start","end","score","strand","phase","attr")
gff$start <- as.integer(gff$start)
gff$end <- as.integer(gff$end)
gff$strand <- as.character(gff$strand)

gff_info <- lapply(gff$attr, parse_gff_attrs)
gff$te_key <- vapply(gff_info, `[[`, character(1), "te_key")
gff$gene <- vapply(gff_info, `[[`, character(1), "gene")

# keep only rows that have a te_key
gff <- gff[!is.na(gff$te_key) & gff$te_key != "", , drop=FALSE]

# -----------------------------
# Build element feature table (relative coordinates)
# -----------------------------
# Map TE coordinates from TSV (source of element boundaries)
te_meta <- tsv[, c("key","chrom","s","e","full_len","ltr_len","lineage")]
rownames(te_meta) <- te_meta$key

# For each GFF domain, shift by SCN lLTR and convert to relative coords within element
# If missing SCN entry, shift=0 (but you said SCN used for true shift; better to warn)
shift_for <- function(key){
  if (!is.na(scn_shift[key])) return(as.integer(scn_shift[key]))
  return(0L)
}

# Create long-format features: LTR5, domains, LTR3
features_long <- list()

# Add LTRs from TSV
ltr_rows <- data.frame(
  te_key = te_meta$key,
  lineage = te_meta$lineage,
  feature = "LTR5",
  rel_start = 1L,
  rel_end = as.integer(te_meta$ltr_len),
  full_len = as.integer(te_meta$full_len),
  stringsAsFactors=FALSE
)
features_long[[length(features_long)+1]] <- ltr_rows

ltr3_start <- as.integer(te_meta$full_len - te_meta$ltr_len + 1)
ltr3_rows <- data.frame(
  te_key = te_meta$key,
  lineage = te_meta$lineage,
  feature = "LTR3",
  rel_start = ltr3_start,
  rel_end = as.integer(te_meta$full_len),
  full_len = as.integer(te_meta$full_len),
  stringsAsFactors=FALSE
)
features_long[[length(features_long)+1]] <- ltr3_rows

# Add domains from GFF
dom_rows <- gff
# attach element boundaries (s,e) and lineage
dom_rows$elem_s <- te_meta[dom_rows$te_key, "s"]
dom_rows$elem_e <- te_meta[dom_rows$te_key, "e"]
dom_rows$lineage <- te_meta[dom_rows$te_key, "lineage"]
dom_rows$full_len <- te_meta[dom_rows$te_key, "full_len"]

# apply SCN shift
dom_rows$shift <- vapply(dom_rows$te_key, shift_for, integer(1))
dom_rows$start_shifted <- dom_rows$start + dom_rows$shift
dom_rows$end_shifted <- dom_rows$end + dom_rows$shift

# relative coordinates (respect strand of the element as given by the domain hit)
dom_rows$rel_start <- ifelse(dom_rows$strand == "+",
                            dom_rows$start_shifted - dom_rows$elem_s + 1,
                            dom_rows$elem_e - dom_rows$end_shifted + 1)
dom_rows$rel_end <- ifelse(dom_rows$strand == "+",
                           dom_rows$end_shifted - dom_rows$elem_s + 1,
                           dom_rows$elem_e - dom_rows$start_shifted + 1)

# ensure rel_start <= rel_end
swap <- dom_rows$rel_start > dom_rows$rel_end
if (any(swap, na.rm=TRUE)) {
  tmp <- dom_rows$rel_start[swap]
  dom_rows$rel_start[swap] <- dom_rows$rel_end[swap]
  dom_rows$rel_end[swap] <- tmp
}

dom_keep <- dom_rows[, c("te_key","lineage","gene","rel_start","rel_end","full_len")]
colnames(dom_keep)[3] <- "feature"
features_long[[length(features_long)+1]] <- dom_keep

features_long <- do.call(rbind, features_long)

# Drop domains that fall outside the element after shifting (keep only sensible ones)
features_long <- features_long[
  !is.na(features_long$rel_start) & !is.na(features_long$rel_end) &
    features_long$rel_start >= 1 & features_long$rel_end <= features_long$full_len &
    features_long$rel_end >= features_long$rel_start,
  , drop=FALSE
]

# If multiple hits for same TE+feature, keep the best "widest" one (simple heuristic)
features_long$width <- features_long$rel_end - features_long$rel_start + 1
features_long <- features_long[order(features_long$te_key, features_long$feature, -features_long$width),]
features_long <- features_long[!duplicated(features_long[,c("te_key","feature")]),]

# -----------------------------
# Plotting
# -----------------------------
feature_colors <- c(
  "LTR5"="#333333",
  "LTR3"="#333333",
  "GAG" ="#fc850e",
  "PROT"="#f80bfb",
  "RT"  ="#0808f7",
  "RH"  ="#fc1413",
  "INT" ="#05fc09",
  "CH"  ="#7a001a"
)

get_col <- function(f){
  if (!is.na(feature_colors[f])) return(feature_colors[f])
  "#AAAAAA"
}

# Compute "average element" coordinates for a lineage by median starts/ends
compute_average_coords <- function(df_lineage){
  # df_lineage: rows are features of elements (te_key, feature, rel_start, rel_end, full_len)
  feats <- unique(df_lineage$feature)
  out <- lapply(feats, function(f){
    sub <- df_lineage[df_lineage$feature == f, , drop=FALSE]
    data.frame(
      feature=f,
      start=as.numeric(stats::median(sub$rel_start, na.rm=TRUE)),
      end=as.numeric(stats::median(sub$rel_end, na.rm=TRUE)),
      stringsAsFactors=FALSE
    )
  })
  out <- do.call(rbind, out)

  # ensure in canonical order where possible
  ord <- c("LTR5","GAG","PROT","RT","RH","INT","CH","LTR3")
  out$order <- match(out$feature, ord)
  out$order[is.na(out$order)] <- 999
  out <- out[order(out$order, out$start), , drop=FALSE]

  avg_full_len <- as.numeric(stats::median(unique(df_lineage[,c("te_key","full_len")])$full_len, na.rm=TRUE))
  list(coords=out[,c("feature","start","end"),drop=FALSE], avg_full_len=avg_full_len)
}

plot_average_element <- function(df_lineage, lineage_name, file_png){
  avg <- compute_average_coords(df_lineage)
  coords <- avg$coords
  max_x <- max(coords$end, avg$avg_full_len, na.rm=TRUE)

  png(file_png, width=1400, height=450, pointsize=18)
  par(mar=c(5,5,4,2))
  plot(0,0, type="n", xlim=c(0, max_x*1.05), ylim=c(0, 10),
       xlab="Position [bp]", ylab="", axes=FALSE,
       main=paste0(lineage_name, " : Structure of average element"))
  axis(1, las=1)
  rect(0, 5, max_x, 7, col="#DDDDDD", border=NA)

  for (i in seq_len(nrow(coords))){
    f <- coords$feature[i]
    xs <- coords$start[i]
    xe <- coords$end[i]
    rect(xs, 4.8, xe, 7.2, col=get_col(f), border=NA)
    lab <- if (f=="LTR5") "5'LTR" else if (f=="LTR3") "3'LTR" else f
    text((xs+xe)/2, 3.5, labels=lab, cex=1)
  }

  # legend
  leg_feats <- intersect(c("LTR5","LTR3","GAG","PROT","RT","RH","INT","CH"), unique(coords$feature))
  legend("topright",
         legend=gsub("LTR5","LTR", gsub("LTR3","LTR", leg_feats)),
         fill=sapply(leg_feats, get_col),
         bty="n", cex=1.1)
  dev.off()
}

plot_individual_elements <- function(df_lineage, lineage_name, file_png, maxN=10000){
  # Build per-element rows of rectangles (one row per TE)
  # Optionally subsample to maxN elements
  tes <- unique(df_lineage$te_key)
  if (length(tes) == 0) return()

  if (length(tes) > maxN) {
    set.seed(1)
    tes <- sample(tes, maxN)
    df_lineage <- df_lineage[df_lineage$te_key %in% tes, , drop=FALSE]
  }

  # Order elements using hclust on feature centers (roughly like your original approach)
  all_features <- c("LTR5","GAG","PROT","RT","RH","INT","CH","LTR3")
  centers <- matrix(NA_real_, nrow=length(all_features), ncol=length(tes),
                    dimnames=list(all_features, tes))
  widths <- matrix(NA_real_, nrow=length(all_features), ncol=length(tes),
                   dimnames=list(all_features, tes))

  for (ti in tes){
    sub <- df_lineage[df_lineage$te_key==ti, , drop=FALSE]
    for (f in intersect(sub$feature, all_features)){
      row <- sub[sub$feature==f,][1,]
      centers[f, ti] <- (row$rel_start + row$rel_end)/2
      widths[f, ti] <- (row$rel_end - row$rel_start + 1)
    }
  }

  # Replace NA centers with -1 for distance calculation (keeps “missing” informative)
  centers2 <- centers
  centers2[is.na(centers2)] <- -1
  ord <- tryCatch({
    hclust(dist(t(centers2), method="manhattan"), method="average")$order
  }, error=function(e) seq_along(tes))
  tes_ord <- tes[ord]

  max_len <- max(df_lineage$full_len, na.rm=TRUE)

  png(file_png, width=1200, height=1000, pointsize=18)
  par(mar=c(5,5,4,2))
  plot(0,0, type="n", xlim=c(0, max_len*1.05), ylim=c(0, length(tes_ord)+1),
       xlab="Position [bp]", ylab="", axes=FALSE,
       main=paste0(lineage_name, " : Structure of individual elements (N=", length(tes_ord), ")"))
  axis(1, las=1)

  for (i in seq_along(tes_ord)){
    te <- tes_ord[i]
    yb <- i - 0.4
    yt <- i + 0.4

    # baseline full element bar (light grey)
    fl <- unique(df_lineage[df_lineage$te_key==te, "full_len"])[1]
    rect(0, yb, fl, yt, col="#EEEEEE", border=NA)

    sub <- df_lineage[df_lineage$te_key==te, , drop=FALSE]
    # draw features
    for (j in seq_len(nrow(sub))){
      f <- sub$feature[j]
      xs <- sub$rel_start[j]
      xe <- sub$rel_end[j]
      rect(xs, yb, xe, yt, col=get_col(f), border=NA)
    }
  }

  # legend
  leg_feats <- intersect(all_features, unique(df_lineage$feature))
  legend("topright",
         legend=gsub("LTR5","LTR", gsub("LTR3","LTR", leg_feats)),
         fill=sapply(leg_feats, get_col),
         bty="n", cex=1.1)
  dev.off()
}

plot_all_elements <- function(df_all, file_png){
  lineages <- sort(unique(df_all$lineage))
  if (length(lineages) == 0) return()

  # precompute averages per lineage
  avgs <- lapply(lineages, function(LN){
    compute_average_coords(df_all[df_all$lineage==LN, , drop=FALSE])
  })
  names(avgs) <- lineages

  # global xlim based on the *max* lineage avg length (not feature ends)
  max_x <- max(vapply(avgs, function(x) x$avg_full_len, numeric(1)), na.rm=TRUE)

  png(file_png, width=1500, height=max(600, 80*length(lineages)), pointsize=18)
  par(mar=c(5,10,4,2))

  plot(0,0, type="n",
       xlim=c(0, max_x*1.05),
       ylim=c(0, length(lineages)+1),
       xlab="Position [bp]", ylab="", axes=FALSE,
       main="all_elements.png : Average structure per lineage")

  axis(1, las=1)

  for (i in seq_along(lineages)){
    LN <- lineages[i]
    coords <- avgs[[LN]]$coords
    avg_len <- avgs[[LN]]$avg_full_len

    y <- length(lineages) - i + 1
    yb <- y - 0.25
    yt <- y + 0.25

    # optional faint full-width row background (just a canvas)
    rect(0, yb, max_x, yt, col="#FBFBFB", border=NA)

    # THIS is the real element-length bar (varies by lineage!)
    rect(0, yb, avg_len, yt, col="#DDDDDD", border=NA)

    # draw average features, clipped to avg_len
    for (j in seq_len(nrow(coords))){
      f <- coords$feature[j]
      xs <- coords$start[j]
      xe <- coords$end[j]
      if (is.na(xs) || is.na(xe)) next

      xs <- max(1, xs)
      xe <- min(avg_len, xe)
      if (xe < xs) next

      rect(xs, yb, xe, yt, col=get_col(f), border=NA)
    }

    text(-max_x*0.02, y, labels=LN, adj=1, xpd=NA, cex=1.1)
  }

  leg_feats <- intersect(c("LTR5","LTR3","GAG","PROT","RT","RH","INT","CH"),
                         unique(df_all$feature))
  legend("topright",
         legend=gsub("LTR5","LTR", gsub("LTR3","LTR", leg_feats)),
         fill=sapply(leg_feats, get_col),
         bty="n", cex=1.1)

  dev.off()
}

# -----------------------------
# Output
# -----------------------------
out_dir <- paste0(opt$output, "_plots")
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

# Per-lineage plots
lineages <- sort(unique(features_long$lineage))
for (LN in lineages){
  dfL <- features_long[features_long$lineage == LN, , drop=FALSE]
  safe <- gsub("[^[:alnum:]_\\-\\.]", "_", LN)

  avg_png <- file.path(out_dir, paste0(safe, "_average.png"))
  ind_png <- file.path(out_dir, paste0(safe, "_individual.png"))

  plot_average_element(dfL, LN, avg_png)
  plot_individual_elements(dfL, LN, ind_png, maxN=opt$max_individual)
}

# Global all_elements.png
plot_all_elements(features_long, file.path(out_dir, "all_elements.png"))

message("Done.")
message("Outputs in: ", out_dir)
