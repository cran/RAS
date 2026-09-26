# ---------------------------------------------------------------------------
# Summary-informed RAS: use an EXTERNAL GWAS's effect estimates as the Stage-1
# weights, and scan on ALL individuals of the target cohort.
#
# The internal pipeline (ras_scan_original()/ras_scan()) splits the target cohort
# 50/50, estimates SNP effects on one half, scans the other, and repeats. That
# split exists only because the same individuals cannot both estimate an effect
# and be tested on it. When the effects come from an INDEPENDENT cohort the
# conflict is gone: the weights are fixed, no split is drawn, and every target
# individual contributes to the regional scan.
#
#   ras_harmonize_sumstats()   external file + local map -> weights + a QC table
#   ras_weights_from_sumstats()  the matcher/aligner underneath it
#   ras_scan_external()   weights -> RAS profile on the full cohort (C path)
#   ras_scan_external_original()        the same, in memory, for small cohorts
#
# Harmonisation contract:
#   * `match_by` is EXPLICIT -- there is no silent fall-through from ID to
#     position matching. The default c("id","pos") uses both on purpose: an
#     rsID is build-invariant but not unique (multi-allelic sites share one),
#     chr:pos:alleles is unique but build-dependent, and not every variant has
#     an rsID. The ID-matched subset is also the only thing the build gate can
#     verify positions on, so it is what licenses the position matching.
#   * `build_check` compares, for every SNP matched by ID, the position the
#     summary statistics claim against the position the local map claims. On
#     the same build the difference is identically zero; across builds it is
#     tens of kb, and a PIECEWISE offset (the dangerous case, because it
#     truncates the weights spatially rather than uniformly) shows up as a
#     fraction of non-zero differences strictly between 0 and 1.
#   * `match_min_prop` refuses to return a weight vector that matched too
#     little of the smaller dataset (bigsnpr's guard, default 20%).
#   * Strand-ambiguous (A/T, C/G) variants are REMOVED, not rescued. Rescuing
#     them needs a frequency rule, and a frequency rule is a modelling decision
#     this layer deliberately does not make.
#   * `blocks` reports the match rate along the chromosome. RAS's statistic is
#     REGIONAL, so "65% overall, uniformly" and "100% then 30%" are two
#     completely different inputs -- a distinction a PRS does not care about.
#
# Genome build / naming are two separate problems, often confused:
#   (i)  NAMING  -- both sides on the same build, but one says "rs1558902" and
#                   the other says "chr16:53769311:A:G". Fixed by a dictionary
#                   LOOKUP. No coordinate is altered.
#   (ii) BUILD   -- the same variant carries different coordinates on each side.
#                   Fixed by a coordinate CONVERSION, or sidestepped by a
#                   dictionary, since looking rs1558902 up in a GRCh38 table
#                   RETURNS its GRCh38 position without mapping anything.
# Three routes, chosen from the data rather than from preference:
#   B0  ras_map_ids_to_rsid()  rename the LOCAL map's variants to rsIDs using a
#       dictionary on the map's own build. Emits a PLINK --update-name file.
#       PREFERRED when the map is ours.
#   B1  dictionary translation of the SUMSTATS into the map's naming and build.
#       A table lookup only ever FAILS to find; it cannot find the wrong thing,
#       which is why it beats coordinate mapping.
#   B2  ras_liftover_sumstats()  last resort: liftOver has real error rates
#       (~1.7% GRCh37->38 position discordance) so this route always round-trips
#       (from -> to -> from) and drops anything that fails.
# PLINK itself does not convert builds -- it only applies a conversion someone
# else computed (--update-map / --update-name / --ref-from-fa).
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Column-name synonyms seen in the wild: PLINK1 .assoc.linear / .qassoc,
# PLINK2 .glm.linear, GWAS-SSF (the GWAS Catalog standard), and the usual public
# conventions. The `hm_*` names come FIRST on purpose: in a GWAS Catalog
# harmonised file those columns are the authoritative GRCh38 / forward-strand
# values, and the un-prefixed columns are the submitter's originals.
.SYN <- list(
  snp    = c("hm_rsid","hm_variant_id",
             "SNP","ID","rsid","RSID","rs_id","MarkerName","variant_id","snpid"),
  a1     = c("hm_effect_allele",
             "A1","ALT","EA","effect_allele","EFFECT_ALLELE","Allele1","A1_EFFECT",
             "Tested_Allele"),
  a2     = c("hm_other_allele",
             "A2","REF","OA","NEA","other_allele","OTHER_ALLELE","Allele2","AX",
             "Other_Allele"),
  beta   = c("hm_beta",
             "BETA","B","beta","EFFECT","effect","Effect","log_OR","LOG_OR"),
  or     = c("hm_odds_ratio","OR","odds_ratio","OddsRatio"),
  se     = c("SE","StdErr","standard_error","se"),
  p      = c("P","PVAL","pval","P_BOLT_LMM","p_value","P.value"),
  chr    = c("hm_chrom","CHR","chr","chrom","CHROM","#CHROM","chromosome","Chr"),
  pos    = c("hm_pos","POS","BP","pos","bp","position","base_pair_location","Position"),
  eaf    = c("hm_effect_allele_frequency","effect_allele_frequency","EAF","FRQ",
             "freq","A1FREQ","Freq_Tested_Allele_in_HRS","MAF")
)

# Dictionary column synonyms (dbSNP slice, All of Us Variant Annotation Table, VCF).
.DSYN <- list(
  rsid = c("RSID","rsid","rs_id","ID","SNP","dbsnp_rsid","variant_rsid","rs"),
  chr  = c("CHR","chr","CHROM","#CHROM","chrom","chromosome","contig"),
  pos  = c("BP","POS","pos","position","start","location"),
  ref  = c("REF","ref","reference_allele","A2","other_allele"),
  alt  = c("ALT","alt","alternate_allele","alternate_bases","A1","effect_allele")
)

.RE_RSID <- "^rs[0-9]+$"

# Normalise a chromosome label: drop any "chr" prefix, map X/Y/MT to numbers so
# that "chr16", "16" and 16 all collide on the same key.
.norm_chr <- function(x) {
  x <- toupper(trimws(as.character(x))); x <- sub("^CHR", "", x)
  x[x == "X"] <- "23"; x[x == "Y"] <- "24"; x[x %in% c("MT","M")] <- "26"; x
}

# Unordered allele-pair key, so A/G and G/A collide (used to disambiguate
# multi-allelic sites when matching on position rather than on ID).
.allele_key <- function(a1, a2) {
  a1 <- toupper(a1); a2 <- toupper(a2)
  paste(pmin(a1, a2), pmax(a1, a2), sep = "/")
}

.pick  <- function(nms, key) { hit <- intersect(.SYN[[key]], nms); if (length(hit)) hit[1] else NA_character_ }
.flip  <- function(a) c(A = "T", T = "A", C = "G", G = "C")[a]
.ambig <- function(a1, a2) (a1 == .flip(a2))          # A/T or C/G

.read_delim_table <- function(x) {
  utils::read.table(x, header = TRUE, stringsAsFactors = FALSE,
                    comment.char = "", check.names = FALSE)
}

# PLINK --linear/--logistic with --covar emits one row per TEST term (ADD, then
# each covariate). Keep only the additive SNP term.
.keep_add_rows <- function(ss, verbose = FALSE) {
  if (!("TEST" %in% names(ss))) return(ss)
  n0 <- nrow(ss)
  ss <- ss[trimws(as.character(ss$TEST)) == "ADD", , drop = FALSE]
  if (verbose)
    message(sprintf("[sumstats] TEST column found: kept %d ADD rows of %d", nrow(ss), n0))
  ss
}

# --- build-consistency gate -------------------------------------------------
# Only meaningful for SNPs matched by ID: those carry an independent position
# claim from each side. Returns a verdict; the caller decides how loud to be,
# because a build mismatch is harmless when positions are never used for
# matching and fatal when they are.
.build_verdict <- function(ss_bp, map_bp, min_n = 20L, tol_frac = 0.01) {
  ok <- is.finite(ss_bp) & is.finite(map_bp)
  n  <- sum(ok)
  if (n < min_n)
    return(list(evaluable = FALSE, n = n, reason = sprintf(
      "only %d ID-matched SNPs carry positions on both sides (need %d)", n, min_n)))
  d <- ss_bp[ok] - map_bp[ok]
  med <- stats::median(d); frac <- mean(d != 0)
  list(evaluable = TRUE, n = n, median_dbp = med, frac_nonzero = frac,
       consistent = (med == 0 && frac <= tol_frac),
       # a fraction strictly inside (tol, 1) means part of the chromosome agrees
       # and part does not: the piecewise-offset signature.
       piecewise = (frac > tol_frac && frac < 0.99),
       quantiles = stats::quantile(d, c(0, .25, .5, .75, 1)))
}

# --- match rate along the chromosome ---------------------------------------
.match_blocks <- function(usable, bp, block_size = 500L) {
  n <- length(usable)
  if (n == 0L) return(NULL)
  grp <- ((seq_len(n) - 1L) %/% block_size) + 1L
  data.frame(
    block     = as.integer(sort(unique(grp))),
    idx_start = as.integer(tapply(seq_len(n), grp, min)),
    idx_end   = as.integer(tapply(seq_len(n), grp, max)),
    bp_start  = if (!is.null(bp)) as.numeric(tapply(bp, grp, function(z) suppressWarnings(min(z, na.rm = TRUE)))) else NA_real_,
    bp_end    = if (!is.null(bp)) as.numeric(tapply(bp, grp, function(z) suppressWarnings(max(z, na.rm = TRUE)))) else NA_real_,
    n_snps    = as.integer(tapply(usable, grp, length)),
    rate      = as.numeric(tapply(usable, grp, mean)),
    row.names = NULL, stringsAsFactors = FALSE)
}

# ===========================================================================
# ID conventions
# ===========================================================================
#' Classify a vector of variant IDs (internal).
#' @return list(style, prop, prefix, sep, nfield, example) where style is one of
#'   "rsid", "chr_pos", "chr_pos_alleles", "other".
#' @noRd
.id_style <- function(ids, n_sample = 5000L) {
  ids <- as.character(ids); ids <- ids[!is.na(ids) & nzchar(ids)]
  if (!length(ids)) return(list(style = "other", prop = 0, prefix = "", sep = ":",
                                nfield = 1L, example = NA_character_))
  s <- if (length(ids) > n_sample) sample(ids, n_sample) else ids

  is_rs <- grepl(.RE_RSID, s)
  # split on the first separator that actually occurs
  sep <- if (any(grepl(":", s, fixed = TRUE))) ":" else
         if (any(grepl("_", s, fixed = TRUE))) "_" else ":"
  parts <- strsplit(s, sep, fixed = TRUE)
  nf <- vapply(parts, length, integer(1))
  chr_ok <- vapply(parts, function(p) {
    if (!length(p)) return(FALSE)
    grepl("^(chr)?([0-9]{1,2}|X|Y|M|MT)$", p[1], ignore.case = TRUE) }, logical(1))
  pos_ok <- vapply(parts, function(p) length(p) >= 2 && grepl("^[0-9]+$", p[2]), logical(1))
  al_ok  <- vapply(parts, function(p) length(p) >= 4 &&
                     all(grepl("^[ACGTN]+$", toupper(p[3:4]))), logical(1))

  cand <- c(rsid            = mean(is_rs),
            chr_pos_alleles = mean(chr_ok & pos_ok & al_ok & nf >= 4),
            chr_pos         = mean(chr_ok & pos_ok & nf == 2))
  style <- names(cand)[which.max(cand)]
  prop  <- max(cand)
  if (prop < 0.5) { style <- "other"; prop <- 1 - max(cand) }
  prefix <- if (style %in% c("chr_pos","chr_pos_alleles") &&
                mean(grepl("^chr", s, ignore.case = TRUE)) > 0.5) "chr" else ""
  list(style = style, prop = prop, prefix = prefix, sep = sep,
       nfield = as.integer(stats::median(nf)), example = s[1])
}

#' Decide the alignment route from the two ID conventions (internal).
#'
#' This answers the NAMING question only. Whether the builds agree is a separate
#' question, settled by ras_weights_from_sumstats()'s build gate (which needs ID
#' overlap to run) or by external knowledge.
#' @noRd
.plan_alignment <- function(map_ids, ss_ids) {
  m <- .id_style(map_ids); e <- .id_style(ss_ids)
  route <- if (m$style == "rsid" && e$style == "rsid") "id" else
           if (m$style %in% c("chr_pos","chr_pos_alleles") && e$style == "rsid") "dict" else
           if (m$style == "rsid" && e$style != "rsid") "dict" else "manual"
  note <- switch(route,
    id   = "IDs already agree: no translation and no conversion of any kind.",
    dict = paste("IDs disagree. A dictionary on the MAP's build fixes the naming and",
                 "supplies the map's coordinates in one lookup. Position matching",
                 "alone would also work IF the builds are known to agree, but with",
                 "no ID overlap the build gate cannot verify that for you."),
    manual = "Neither side uses a recognised convention; supply a dictionary keyed on the map's own IDs.")
  why <- sprintf("map IDs look like %s (%.0f%%), sumstats IDs look like %s (%.0f%%)",
                 m$style, 100*m$prop, e$style, 100*e$prop)
  list(route = route, map_style = m, ss_style = e, why = why, note = note)
}

#' Build map-style IDs from chr/pos/allele columns (internal).
#' @noRd
.render_ids <- function(chr, pos, a3, a4, style) {
  ch <- paste0(style$prefix, .norm_chr(chr))
  if (identical(style$style, "chr_pos"))
    paste(ch, pos, sep = style$sep)
  else
    paste(ch, pos, toupper(a3), toupper(a4), sep = style$sep)
}

#' Pick the allele order (ref:alt vs alt:ref) that reproduces more of the map's
#' own IDs (internal). The order is not recoverable from an ID string alone, so
#' it is settled empirically against the map rather than assumed.
#' @noRd
.choose_allele_order <- function(dict, map_ids, style, n_probe = 20000L) {
  if (!identical(style$style, "chr_pos_alleles")) return("ref_alt")
  k <- if (nrow(dict) > n_probe) sample(nrow(dict), n_probe) else seq_len(nrow(dict))
  d <- dict[k, , drop = FALSE]
  hit_ra <- mean(.render_ids(d$CHR, d$BP, d$REF, d$ALT, style) %in% map_ids)
  hit_ar <- mean(.render_ids(d$CHR, d$BP, d$ALT, d$REF, style) %in% map_ids)
  message(sprintf("[ids] allele order probe: ref:alt %.2f%% vs alt:ref %.2f%% of dictionary rows hit the map",
                  100*hit_ra, 100*hit_ar))
  if (hit_ar > hit_ra) "alt_ref" else "ref_alt"
}

# ===========================================================================
# The dictionary
# ===========================================================================
#' Read a Variant Dictionary for Summary-Statistics Harmonisation
#'
#' Reads a table mapping rsID to (chromosome, position, REF, ALT) on ONE genome
#' build, for use by \code{\link{ras_harmonize_sumstats}}'s dictionary route and
#' by \code{\link{ras_map_ids_to_rsid}}. Sources include a dbSNP slice, the All
#' of Us Variant Annotation Table, or any delimited file with those columns.
#'
#' A plain VCF is accepted: lines beginning \code{##} are skipped and
#' \code{#CHROM} is used as the header. An \code{ALT} field carrying several
#' comma-separated alleles is expanded to one row per alternate allele, so
#' multi-allelic sites resolve by their alleles rather than by row order.
#'
#' @param src data.frame, or a path to a delimited file / VCF.
#' @param cols Optional named overrides for the auto-detected columns, e.g.
#'   \code{list(rsid = "dbsnp_rsid")}. Names are \code{rsid}, \code{chr},
#'   \code{pos}, \code{ref}, \code{alt}.
#' @param build Free-text label recorded on the result for provenance, e.g.
#'   \code{"GRCh38"}. Never inferred and never acted upon -- it is an audit
#'   note, while the build gate does the actual checking.
#'
#' @return A data frame with columns \code{RSID}, \code{CHR}, \code{BP},
#'   \code{REF}, \code{ALT}, \code{KEY} (the unordered allele-pair key), with
#'   the \code{build} label attached as an attribute.
#'
#' @seealso \code{\link{ras_harmonize_sumstats}}, \code{\link{ras_map_ids_to_rsid}}
#' @export
ras_read_variant_dictionary <- function(src, cols = list(), build = NA_character_) {
  d <- if (is.character(src)) {
    hdr <- readLines(src, n = 2000L, warn = FALSE)
    skip <- sum(grepl("^##", hdr))
    utils::read.table(src, header = TRUE, skip = skip, comment.char = "",
                      stringsAsFactors = FALSE, check.names = FALSE,
                      sep = if (grepl("\t", hdr[skip + 1L], fixed = TRUE)) "\t" else "")
  } else as.data.frame(src)
  nm <- names(d)
  pk <- function(k) cols[[k]] %||% { h <- intersect(.DSYN[[k]], nm); if (length(h)) h[1] else NA_character_ }
  cn <- lapply(c("rsid","chr","pos","ref","alt"), pk); names(cn) <- c("rsid","chr","pos","ref","alt")
  miss <- names(cn)[vapply(cn, is.na, logical(1))]
  if (length(miss)) stop("variant dictionary is missing column(s): ",
                         paste(miss, collapse = ", "),
                         "; supply them via cols=list(...)", call. = FALSE)
  out <- data.frame(RSID = as.character(d[[cn$rsid]]),
                    CHR  = .norm_chr(d[[cn$chr]]),
                    BP   = suppressWarnings(as.integer(d[[cn$pos]])),
                    REF  = toupper(as.character(d[[cn$ref]])),
                    ALT  = toupper(as.character(d[[cn$alt]])),
                    stringsAsFactors = FALSE)
  # multi-allelic ALT -> one row per alternate allele
  if (any(grepl(",", out$ALT, fixed = TRUE))) {
    sp <- strsplit(out$ALT, ",", fixed = TRUE)
    k  <- rep(seq_len(nrow(out)), lengths(sp))
    out <- out[k, , drop = FALSE]; out$ALT <- unlist(sp, use.names = FALSE)
  }
  out <- out[!is.na(out$RSID) & nzchar(out$RSID) & is.finite(out$BP), , drop = FALSE]
  out$KEY <- .allele_key(out$REF, out$ALT)
  attr(out, "build") <- build
  rownames(out) <- NULL
  message(sprintf("[dict] %d entries%s | %d distinct rsIDs",
                  nrow(out), if (is.na(build)) "" else sprintf(" (%s)", build),
                  length(unique(out$RSID))))
  out
}

# ===========================================================================
# Route B0 -- put rsIDs onto the LOCAL map
# ===========================================================================
#' @noRd
.write_plink_update_name <- function(update_name, file) {
  utils::write.table(update_name, file, row.names = FALSE, col.names = FALSE,
                     quote = FALSE, sep = "\t")
  message(sprintf("[B0] wrote %d rename pairs to %s  (plink --update-name %s)",
                  nrow(update_name), file, basename(file)))
  invisible(file)
}

#' Rename a Local Genotype Map's Variants to rsIDs (Harmonisation Route B0)
#'
#' Gives the LOCAL map rsIDs using a dictionary on the SAME genome build, so
#' that afterwards both sides are keyed on rsID, the build gate becomes
#' evaluable, and nothing about positions has to be trusted. Preferred when the
#' map is ours -- inside All of Us the dictionary is the Variant Annotation
#' Table.
#'
#' Matching is chr:pos plus the unordered allele pair, so multi-allelic sites
#' cannot be confused. Variants with no dictionary entry keep their original ID
#' (they simply will not match the summary statistics afterwards).
#'
#' @param map Data frame with columns \code{SNP}, \code{CHR}, \code{BP},
#'   \code{A1}, \code{A2}: the local genotype map, one row per genotype column,
#'   in order.
#' @param dict A dictionary from \code{\link{ras_read_variant_dictionary}}, on
#'   the map's own build.
#' @param update_name_file Optional path. When given, a two-column PLINK
#'   \code{--update-name} file (old ID, new ID; no header) is written there.
#'
#' @return A list with \code{map} (the renamed map, carrying the original IDs
#'   in \code{SNP_original}), \code{report} (counts), and \code{update_name}
#'   (the rename pairs as a data frame).
#'
#' @seealso \code{\link{ras_harmonize_sumstats}}
#' @export
ras_map_ids_to_rsid <- function(map, dict, update_name_file = NULL) {
  stopifnot(all(c("SNP","CHR","BP","A1","A2") %in% names(map)))
  key_map  <- paste(.norm_chr(map$CHR), as.integer(map$BP),
                    .allele_key(map$A1, map$A2), sep = ":")
  key_dict <- paste(dict$CHR, dict$BP, dict$KEY, sep = ":")
  dup <- duplicated(key_dict)
  i <- match(key_map, key_dict[!dup])
  new_id <- dict$RSID[!dup][i]
  hit <- !is.na(new_id)
  out <- map; out$SNP_original <- map$SNP
  out$SNP[hit] <- new_id[hit]
  upd <- data.frame(old = map$SNP[hit], new = new_id[hit], stringsAsFactors = FALSE)
  rep_ <- data.frame(n_map = nrow(map), n_renamed = sum(hit),
                     pct_renamed = round(100*mean(hit), 2),
                     n_dict_dup_keys = sum(dup), stringsAsFactors = FALSE)
  message(sprintf("[B0] renamed %d of %d local variants to rsIDs (%.1f%%)",
                  sum(hit), nrow(map), 100*mean(hit)))
  if (!is.null(update_name_file)) .write_plink_update_name(upd, update_name_file)
  list(map = out, report = rep_, update_name = upd)
}

# ===========================================================================
# Route B1 -- translate the SUMSTATS into the map's build and naming
# ===========================================================================
#' Give an external summary-statistics table the local map's coordinates and IDs
#' (internal).
#'
#' The dictionary must be on the MAP's build. rsIDs are the join key, so no
#' coordinate arithmetic happens: an rsID that dbSNP has retired or remapped
#' simply fails to resolve instead of silently landing on a neighbour.
#' @noRd
.sumstats_to_map_ids <- function(ss, dict, map, cols = list()) {
  if (is.character(ss)) ss <- .read_delim_table(ss)
  ss <- .keep_add_rows(as.data.frame(ss))
  nm <- names(ss)
  c_snp <- cols$snp %||% .pick(nm, "snp")
  c_a1  <- cols$a1  %||% .pick(nm, "a1")
  c_a2  <- cols$a2  %||% .pick(nm, "a2")
  if (is.na(c_snp)) stop("no SNP-ID column in the summary statistics", call. = FALSE)
  if (is.na(c_a1))  stop("no effect-allele column in the summary statistics", call. = FALSE)

  n0 <- nrow(ss)
  rs <- as.character(ss[[c_snp]])
  a1 <- toupper(as.character(ss[[c_a1]]))
  a2 <- if (!is.na(c_a2)) toupper(as.character(ss[[c_a2]])) else NA_character_
  if (all(is.na(a2)))
    stop("dictionary translation needs BOTH alleles: without A2 the dictionary's\n",
         "  allele pair cannot be checked, and an rsID that dbSNP has remapped\n",
         "  would be accepted silently.", call. = FALSE)
  ss_key <- .allele_key(a1, a2)

  # Join on rsID AND the unordered allele pair in one key, so that:
  #   - an rsID dbSNP has remapped to a different variant fails the allele test
  #     instead of being accepted silently, and
  #   - a multi-allelic rsID is resolved by its alleles rather than by whichever
  #     dictionary row happens to come first.
  # Vectorised: a genome-wide sumstats file is millions of rows.
  dkey <- paste(dict$RSID, dict$KEY)
  amb  <- dkey %in% dkey[duplicated(dkey)]        # same rsID+alleles twice: unusable
  d    <- dict[!amb, , drop = FALSE]
  i    <- match(paste(rs, ss_key), paste(d$RSID, d$KEY))

  cand    <- rs %in% dict$RSID                     # rsID known at all
  n_multi <- sum(is.na(i) & (paste(rs, ss_key) %in% dkey[amb]))
  n_allele_conflict <- sum(is.na(i) & cand) - n_multi

  res_chr <- d$CHR[i]; res_bp <- d$BP[i]
  res_ref <- d$REF[i]; res_alt <- d$ALT[i]
  keep <- !is.na(res_bp)

  style <- .id_style(map$SNP)
  aord  <- .choose_allele_order(dict, as.character(map$SNP), style)
  new_id <- if (style$style == "rsid") rs else if (aord == "ref_alt")
    .render_ids(res_chr, res_bp, res_ref, res_alt, style) else
    .render_ids(res_chr, res_bp, res_alt, res_ref, style)

  out <- ss[keep, , drop = FALSE]
  out[[c_snp]] <- new_id[keep]
  out$CHR <- res_chr[keep]; out$BP <- res_bp[keep]
  # keep the pre-translation coordinates for auditing
  out$BP_source <- if (!is.na(.pick(nm, "pos"))) ss[[.pick(nm, "pos")]][keep] else NA_integer_

  rep_ <- data.frame(
    n_in = n0, n_rsid_absent = sum(!cand),
    n_allele_conflict = n_allele_conflict, n_multi_mapping = n_multi,
    n_out = sum(keep), pct_out = round(100*mean(keep), 2), stringsAsFactors = FALSE)
  message(sprintf(
    "[B1] %d sumstats rows -> %d translated (%.1f%%) | rsID absent %d | allele conflict %d | multi-mapping %d",
    n0, sum(keep), 100*mean(keep), sum(!cand), n_allele_conflict, n_multi))
  list(sumstats = out, report = rep_, id_style = style, allele_order = aord)
}

# ===========================================================================
# Route B2 -- liftOver, with a mandatory round trip
# ===========================================================================
#' Default converter: shell out to the UCSC liftOver binary (internal).
#' Replaceable so the surrounding logic can be tested without the binary.
#' @noRd
.liftover_converter_ucsc <- function(chr, pos, chain, liftover_bin = "liftOver",
                                     tmpdir = tempdir()) {
  n <- length(pos)
  bed_in  <- file.path(tmpdir, "lo_in.bed")
  bed_out <- file.path(tmpdir, "lo_out.bed")
  bed_um  <- file.path(tmpdir, "lo_unmapped.bed")
  utils::write.table(
    data.frame(chr = paste0("chr", .norm_chr(chr)), start = as.integer(pos) - 1L,
               end = as.integer(pos), name = seq_len(n)),
    bed_in, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  rc <- system2(liftover_bin, c(shQuote(bed_in), shQuote(chain),
                                shQuote(bed_out), shQuote(bed_um)),
                stdout = FALSE, stderr = FALSE)
  if (!file.exists(bed_out))
    stop("liftOver produced no output (exit ", rc, "); is `", liftover_bin,
         "` on PATH and is the chain file correct?", call. = FALSE)
  ob <- utils::read.table(bed_out, header = FALSE, stringsAsFactors = FALSE)
  out <- data.frame(idx = as.integer(ob[[4]]),
                    chr_new = .norm_chr(ob[[1]]),
                    pos_new = as.integer(ob[[3]]), stringsAsFactors = FALSE)
  full <- data.frame(idx = seq_len(n), chr_new = NA_character_,
                     pos_new = NA_integer_, stringsAsFactors = FALSE)
  full[out$idx, c("chr_new","pos_new")] <- out[, c("chr_new","pos_new")]
  full
}

#' Convert positions between builds and verify by round-tripping (internal).
#'
#' `check_reverse` implements bigsnpr::snp_modifyBuild()'s guard: anything that
#' cannot return to its original chromosome AND position via from->to->from is
#' discarded. Without it, conversion-unstable positions pass silently, and those
#' are exactly the ones that can land on a different chromosome.
#' @noRd
.liftover_positions <- function(chr, pos, chain, chain_back = NULL,
                                converter = .liftover_converter_ucsc,
                                check_reverse = TRUE, cup = NULL, ...) {
  n <- length(pos); stopifnot(length(chr) == n)
  status <- rep("ok", n)
  chr0 <- .norm_chr(chr); pos0 <- as.integer(pos)

  if (!is.null(cup) && nrow(cup)) {
    cu <- data.frame(CHR = .norm_chr(cup$CHR), START = as.integer(cup$START),
                     END = as.integer(cup$END), stringsAsFactors = FALSE)
    bad <- rep(FALSE, n)
    for (cc in unique(cu$CHR)) {
      k <- which(chr0 == cc); if (!length(k)) next
      iv <- cu[cu$CHR == cc, , drop = FALSE]
      for (r in seq_len(nrow(iv)))
        bad[k[pos0[k] >= iv$START[r] & pos0[k] <= iv$END[r]]] <- TRUE
    }
    status[bad] <- "conversion_unstable"
  }

  go <- which(status == "ok")
  fw <- converter(chr0[go], pos0[go], chain, ...)
  chr_new <- rep(NA_character_, n); pos_new <- rep(NA_integer_, n)
  chr_new[go] <- fw$chr_new; pos_new[go] <- fw$pos_new
  status[go][is.na(fw$pos_new)] <- "unmapped"
  status[!is.na(chr_new) & chr_new != chr0 & status == "ok"] <- "chromosome_jump"

  if (isTRUE(check_reverse)) {
    if (is.null(chain_back))
      stop("check_reverse = TRUE needs `chain_back` (the reverse chain file).\n",
           "  Pass check_reverse = FALSE only if you accept unverified conversions.",
           call. = FALSE)
    bk <- which(status == "ok")
    if (length(bk)) {
      rv <- converter(chr_new[bk], pos_new[bk], chain_back, ...)
      back_ok <- !is.na(rv$pos_new) & rv$pos_new == pos0[bk] &
                 .norm_chr(rv$chr_new) == chr0[bk]
      back_ok[is.na(back_ok)] <- FALSE
      status[bk][!back_ok] <- "round_trip_failed"
    }
  }
  drop <- status != "ok"
  chr_new[drop] <- NA_character_; pos_new[drop] <- NA_integer_
  data.frame(chr_new = chr_new, pos_new = pos_new, status = status,
             stringsAsFactors = FALSE)
}

#' LiftOver External Summary Statistics onto the Map's Build (Route B2)
#'
#' Converts the summary statistics' coordinates to another genome build and
#' verifies every conversion by round-tripping it (\code{from -> to -> from});
#' anything that cannot return to its own chromosome and position is dropped.
#' This is the LAST-RESORT route. liftOver has real error rates (about 1.7% of
#' positions discordant GRCh37 to GRCh38, 3.1% the other way), so prefer a
#' dictionary lookup (\code{\link{ras_read_variant_dictionary}}) when one is
#' available: a lookup can only fail to find a variant, it cannot find the
#' wrong one.
#'
#' Requires the UCSC \code{liftOver} binary on \code{PATH} unless a replacement
#' \code{converter} is supplied.
#'
#' @param ss Data frame of summary statistics, or a path to one.
#' @param chain Path to the forward chain file, e.g. \code{hg19ToHg38.over.chain}.
#' @param chain_back Path to the reverse chain file. Required unless
#'   \code{check_reverse = FALSE}.
#' @param cols Optional column-name overrides, e.g. \code{list(chr = "CHROM")}.
#' @param converter Function performing the conversion; the default shells out
#'   to UCSC \code{liftOver}. Replaceable for testing.
#' @param check_reverse Logical. Verify each conversion by round-tripping it
#'   (default \code{TRUE}). Setting \code{FALSE} accepts unverified conversions.
#' @param cup Optional data frame \code{(CHR, START, END)} of conversion-unstable
#'   intervals on the SOURCE build; overlapping variants are dropped first.
#' @param ... Passed to \code{converter}.
#'
#' @return A list with \code{sumstats} (only successfully converted rows, with
#'   the pre-conversion position kept in \code{BP_source}) and \code{report}
#'   (a per-status count table).
#'
#' @seealso \code{\link{ras_harmonize_sumstats}}
#' @export
ras_liftover_sumstats <- function(ss, chain, chain_back = NULL, cols = list(),
                                  converter = .liftover_converter_ucsc,
                                  check_reverse = TRUE, cup = NULL, ...) {
  if (is.character(ss)) ss <- .read_delim_table(ss)
  ss <- .keep_add_rows(as.data.frame(ss))
  nm <- names(ss)
  c_chr <- cols$chr %||% .pick(nm, "chr"); c_pos <- cols$pos %||% .pick(nm, "pos")
  if (is.na(c_chr) || is.na(c_pos))
    stop("liftOver needs chromosome and position columns in the summary statistics",
         call. = FALSE)
  lo <- .liftover_positions(ss[[c_chr]], ss[[c_pos]], chain, chain_back,
                            converter = converter, check_reverse = check_reverse,
                            cup = cup, ...)
  keep <- lo$status == "ok"
  out <- ss[keep, , drop = FALSE]
  out$BP_source <- ss[[c_pos]][keep]
  out[[c_chr]] <- lo$chr_new[keep]; out[[c_pos]] <- lo$pos_new[keep]
  tb <- as.data.frame(table(lo$status), stringsAsFactors = FALSE)
  names(tb) <- c("status", "n")
  message(sprintf("[B2] %d rows -> %d converted (%.1f%%) | %s",
                  nrow(ss), sum(keep), 100*mean(keep),
                  paste(sprintf("%s %d", tb$status, tb$n), collapse = " | ")))
  list(sumstats = out, report = tb)
}

# ===========================================================================
# The harmonisation QC table
# ===========================================================================
#' Assemble the reproducible harmonisation QC table (internal).
#'
#' One row per stage of the pipeline, each with its own denominator, so the
#' table can be pasted into a manuscript without further arithmetic.
#' @noRd
.harmonisation_qc <- function(n_ss_rows, n_missing_beta, n_dup, n_ss_usable,
                              n_map, n_by_id, n_by_pos, n_matched,
                              status_counts, n_usable, af = NULL) {
  # `[[` on a named vector throws for an absent name, and an absent status is the
  # normal case (a clean file has no allele_mismatch rows at all), so test for
  # membership rather than indexing and hoping.
  g <- function(k) {
    if (!(k %in% names(status_counts))) return(0L)
    v <- status_counts[[k]]
    if (is.null(v) || is.na(v)) 0L else as.integer(v)
  }
  items <- list(
    c("external variants in file",        n_ss_rows,                    "external"),
    c("  dropped: missing/non-finite beta", n_missing_beta,             "external"),
    c("  dropped: duplicate ID+alleles",  n_dup,                        "external"),
    c("external variants usable",         n_ss_usable,                  "external"),
    c("  not present in target",          max(n_ss_usable - n_matched, 0L), "external"),
    c("target variants in map",           n_map,                        "target"),
    c("  matched by rsID",                n_by_id,                      "target"),
    c("  matched by chr:pos",             n_by_pos,                     "target"),
    c("matched target variants",          n_matched,                    "target"),
    c("  absent from summary statistics", g("absent_from_sumstats"),    "target"),
    c("allele-aligned variants",          g("aligned"),                 "target"),
    c("allele-flipped variants",          g("flipped"),                 "target"),
    c("  allele mismatch (dropped)",      g("allele_mismatch"),         "target"),
    c("  allele unresolved (dropped)",    g("unresolved_allele"),       "target"),
    c("palindromic variants removed",     g("strand_ambiguous"),        "target"),
    c("final usable weights",             n_usable,                     "target"))
  qc <- data.frame(
    item  = vapply(items, `[`, character(1), 1),
    n     = as.integer(vapply(items, `[`, character(1), 2)),
    denom = vapply(items, `[`, character(1), 3),
    stringsAsFactors = FALSE)
  den <- ifelse(qc$denom == "external", n_ss_rows, n_map)
  qc$pct <- round(100 * qc$n / pmax(den, 1L), 2)
  if (!is.null(af)) {
    qc <- rbind(qc, data.frame(
      item  = c("allele frequency compared", "  frequency discordant"),
      n     = c(af$n_compared, af$n_discordant),
      denom = c("target", "target"),
      pct   = round(100 * c(af$n_compared, af$n_discordant) / max(n_map, 1L), 2),
      stringsAsFactors = FALSE))
  }
  rownames(qc) <- NULL
  qc
}

# ===========================================================================
# The matcher / aligner
# ===========================================================================
#' Align External GWAS Summary Statistics to a Local Genotype Map
#'
#' Turns an external GWAS's effect estimates into a weight vector aligned to a
#' local genotype matrix, one weight per genotype column, in order. This is the
#' matcher underneath \code{\link{ras_harmonize_sumstats}}; call that instead
#' unless you have already settled the naming and build questions yourself.
#'
#' The alignment rule is deliberately conservative:
#' \itemize{
#'   \item \eqn{w_j = +\beta_j} when the summary statistics' effect allele is the
#'     allele the local dosage counts (\code{map$A1});
#'   \item \eqn{w_j = -\beta_j} when it is the other allele (\code{map$A2});
#'   \item \eqn{w_j = 0} when the variant is absent, the alleles cannot be
#'     reconciled, or (with \code{drop_ambiguous = TRUE}) the variant is
#'     strand-ambiguous.
#' }
#' Strand-ambiguous (A/T, C/G) variants are removed rather than rescued.
#' Rescuing them requires a frequency rule, which is a modelling decision this
#' layer does not make.
#'
#' @param sumstats Path to a whitespace/tab-delimited file, or a data frame.
#'   Column names are auto-detected across PLINK 1/2 and GWAS-SSF conventions;
#'   GWAS Catalog \code{hm_*} columns take precedence, being the authoritative
#'   harmonised values. A PLINK \code{TEST} column is filtered to \code{ADD}.
#' @param map Data frame with columns \code{SNP}, \code{A1}, \code{A2} (and
#'   \code{CHR}, \code{BP} when position matching or the build gate is wanted)
#'   describing the LOCAL genotype matrix, one row per column of the genotype
#'   data, in order. \code{A1} must be the allele the local dosage COUNTS.
#' @param cols Optional named list overriding auto-detected columns, e.g.
#'   \code{list(snp = "MarkerName", a1 = "Allele1", beta = "Effect")}.
#' @param drop_ambiguous Logical. Drop strand-ambiguous A/T and C/G variants
#'   (default \code{TRUE}).
#' @param on_or Logical. If the file carries \code{OR} instead of \code{BETA},
#'   take \code{log(OR)} (default \code{TRUE}).
#' @param match_by \code{"id"}, \code{"pos"}, or \code{c("id","pos")} (the
#'   default) to try IDs first and retry the leftovers on chr:pos. Position
#'   matching REQUIRES the build gate to pass, or to be switched off
#'   deliberately.
#'
#'   The two keys are not interchangeable and neither is redundant. An rsID
#'   survives a change of genome build but is NOT unique -- multi-allelic sites
#'   share one, and a real WGS file has about 1.25\% of its rows in that state.
#'   \code{chr:pos:alleles} is far sharper but is meaningless across builds, and
#'   indel anchoring conventions make it fragile. Not every variant carries an
#'   rsID at all, which is the case position matching exists for. They also
#'   stand in a specific order of dependence: the build gate can only compare
#'   positions on variants BOTH sides name, so the ID-matched subset is what
#'   licenses the position matching that then covers the variants without
#'   rsIDs. Using \code{"pos"} alone leaves nothing able to show it matched the
#'   right variants.
#' @param match_min_prop Numeric. Stop if fewer than this proportion of the
#'   smaller dataset ends up with a usable weight (default 0.2, the guard
#'   \code{bigsnpr::snp_match()} uses). A low match rate is the standard
#'   symptom of a wrong build or an ID-convention mismatch.
#' @param build_check Logical. Verify that ID-matched variants agree on
#'   position (default \code{TRUE}). Set \code{FALSE} only when the builds have
#'   been confirmed by other means; doing so while \code{match_by} includes
#'   \code{"pos"} is the one configuration that can silently produce spatially
#'   truncated weights.
#' @param target_af Optional numeric vector, length \code{nrow(map)}, giving the
#'   frequency of \code{map$A1} in the TARGET cohort. When supplied and the
#'   summary statistics carry an effect-allele frequency column, the two are
#'   compared after orientation and reported. Purely diagnostic: no variant is
#'   ever dropped on a frequency disagreement, because a real frequency
#'   difference between cohorts is expected and is not an error.
#' @param af_tol Numeric. Absolute frequency difference above which a variant is
#'   counted as discordant in that report (default 0.20).
#' @param block_size Integer. Window, in variants along the local map, for the
#'   match-rate profile returned in \code{blocks}.
#' @param pos_fallback Deprecated, kept so old calls still run. \code{TRUE} maps
#'   to \code{match_by = c("id","pos")}, \code{FALSE} to \code{match_by = "id"}.
#'
#' @return A list with
#'   \item{weights}{Numeric vector of length \code{nrow(map)}: the aligned
#'     weights, zero wherever no usable weight exists.}
#'   \item{qc}{The reproducible harmonisation QC table (see
#'     \code{\link{ras_harmonize_sumstats}}).}
#'   \item{report}{Per-status variant counts.}
#'   \item{detail}{One row per map variant: status, matching route, sign, weight.}
#'   \item{matched_by}{Counts matched by ID, by ID+alleles, and by position.}
#'   \item{build}{The build gate's verdict.}
#'   \item{blocks}{Match rate along the chromosome, in \code{block_size} windows.}
#'   \item{af}{The allele-frequency comparison, when \code{target_af} was given.}
#'
#' @seealso \code{\link{ras_harmonize_sumstats}} for the driver that chooses the
#'   naming/build route first; \code{\link{ras_scan_external}} to scan with
#'   the resulting weights.
#' @export
ras_weights_from_sumstats <- function(sumstats, map, cols = list(),
                                      drop_ambiguous = TRUE, on_or = TRUE,
                                      match_by = c("id", "pos"),
                                      match_min_prop = 0.2,
                                      build_check = TRUE,
                                      target_af = NULL,
                                      af_tol = 0.20,
                                      block_size = 500L,
                                      pos_fallback = NULL) {
  if (!is.null(pos_fallback)) {
    match_by <- if (isTRUE(pos_fallback)) c("id", "pos") else "id"
    warning("`pos_fallback` is deprecated; use match_by = ",
            deparse(match_by), call. = FALSE)
  }
  match_by <- unique(as.character(match_by))
  bad <- setdiff(match_by, c("id", "pos"))
  if (length(bad)) stop('match_by must be "id" and/or "pos", got: ',
                        paste(bad, collapse = ", "), call. = FALSE)
  if (!length(match_by)) stop("match_by is empty", call. = FALSE)

  ss <- if (is.character(sumstats)) .read_delim_table(sumstats) else as.data.frame(sumstats)
  ss <- .keep_add_rows(ss, verbose = TRUE)
  nm <- names(ss)
  cn <- list(snp  = cols$snp  %||% .pick(nm,"snp"),
             a1   = cols$a1   %||% .pick(nm,"a1"),
             a2   = cols$a2   %||% .pick(nm,"a2"),
             beta = cols$beta %||% .pick(nm,"beta"),
             or   = cols$or   %||% .pick(nm,"or"),
             chr  = cols$chr  %||% .pick(nm,"chr"),
             pos  = cols$pos  %||% .pick(nm,"pos"),
             eaf  = cols$eaf  %||% .pick(nm,"eaf"))
  if (is.na(cn$snp)) stop("could not find a SNP-ID column in the summary statistics")
  if (is.na(cn$a1))  stop("could not find an effect-allele column (A1/ALT/effect_allele)")
  if (is.na(cn$beta) && is.na(cn$or))
    stop("could not find a BETA or OR column in the summary statistics")

  b <- if (!is.na(cn$beta)) as.numeric(ss[[cn$beta]]) else {
         if (!on_or) stop("only OR present and on_or=FALSE")
         log(as.numeric(ss[[cn$or]])) }
  ssdf <- data.frame(SNP = as.character(ss[[cn$snp]]),
                     A1  = toupper(as.character(ss[[cn$a1]])),
                     A2  = if (!is.na(cn$a2)) toupper(as.character(ss[[cn$a2]])) else NA_character_,
                     BETA = b, stringsAsFactors = FALSE)
  ssdf$CHR <- if (!is.na(cn$chr)) .norm_chr(ss[[cn$chr]]) else NA_character_
  ssdf$BP  <- if (!is.na(cn$pos)) suppressWarnings(as.integer(ss[[cn$pos]])) else NA_integer_
  ssdf$EAF <- if (!is.na(cn$eaf)) suppressWarnings(as.numeric(ss[[cn$eaf]])) else NA_real_

  n_raw <- nrow(ssdf)
  ssdf <- ssdf[is.finite(ssdf$BETA), , drop = FALSE]
  n_nonfinite <- n_raw - nrow(ssdf)
  # Multi-allelic sites share one rsID -- a real WGS file has ~1.25% of its rows
  # in that state. De-duplicating on the ID alone keeps whichever row happened to
  # come first and silently discards the site's other alleles, so a local variant
  # carrying one of those alleles loses its weight for no reason. De-duplicate on
  # ID + unordered allele pair instead: only genuinely redundant rows are dropped.
  ss_has_a2 <- !all(is.na(ssdf$A2))
  ssdf$AKEY <- if (ss_has_a2) .allele_key(ssdf$A1, ssdf$A2) else NA_character_
  dup <- if (ss_has_a2) duplicated(paste(ssdf$SNP, ssdf$AKEY)) else duplicated(ssdf$SNP)
  n_dup <- sum(dup); ssdf <- ssdf[!dup, , drop = FALSE]
  n_multi <- if (ss_has_a2) sum(duplicated(ssdf$SNP)) else 0L

  map$A1 <- toupper(as.character(map$A1)); map$A2 <- toupper(as.character(map$A2))
  have_map_pos <- all(c("CHR","BP") %in% names(map))
  have_ss_pos  <- !all(is.na(ssdf$CHR)) && !all(is.na(ssdf$BP))

  idx <- rep(NA_integer_, nrow(map))
  src <- rep(NA_character_, nrow(map))

  # ---- 1. match by ID ----------------------------------------------------
  # Two tiers. Tier 1 keys on ID + unordered allele pair, so a multi-allelic
  # rsID resolves to the row whose alleles the local map actually carries rather
  # than to whichever row sorted first. Tier 2 falls back to the bare ID for
  # whatever tier 1 missed -- an allele representation the key cannot reconcile
  # still reaches the alignment step below and is scored aligned / flipped /
  # allele_mismatch exactly as before, so nothing that used to match stops.
  n_by_id <- 0L; n_by_id_allele <- 0L
  if ("id" %in% match_by) {
    if (ss_has_a2) {
      idx <- match(paste(map$SNP, .allele_key(map$A1, map$A2)),
                   paste(ssdf$SNP, ssdf$AKEY))
      n_by_id_allele <- sum(!is.na(idx))
      need <- which(is.na(idx))
      if (length(need)) idx[need] <- match(as.character(map$SNP)[need], ssdf$SNP)
    } else {
      idx <- match(as.character(map$SNP), ssdf$SNP)
    }
    n_by_id <- sum(!is.na(idx))
    src[!is.na(idx)] <- "id"
  }

  # ---- 2. build-consistency gate ----------------------------------------
  # Runs off the ID matches, which give an independent position claim from each
  # side. This is the only evidence available BEFORE trusting positions.
  bv <- list(evaluable = FALSE, n = 0L, reason = "no SNPs matched by ID")
  if (!("id" %in% match_by)) bv$reason <- "no ID matching requested"
  if (n_by_id > 0L && have_map_pos && have_ss_pos) {
    k <- which(!is.na(idx))
    bv <- .build_verdict(ssdf$BP[idx[k]], as.integer(map$BP[k]))
  } else if (n_by_id > 0L) {
    bv$reason <- "positions absent from the map and/or the summary statistics"
  }
  if (isTRUE(bv$evaluable)) {
    message(sprintf(
      "[build] %d ID-matched SNPs | median dBP %s | non-zero dBP %.1f%% -> %s",
      bv$n, format(bv$median_dbp, big.mark = ","), 100 * bv$frac_nonzero,
      if (bv$consistent) "SAME BUILD" else if (bv$piecewise)
        "PIECEWISE OFFSET (part of the chromosome agrees, part does not)"
      else "BUILD MISMATCH"))
  }
  if ("pos" %in% match_by && isTRUE(build_check)) {
    if (!isTRUE(bv$evaluable)) {
      stop(sprintf(paste0(
        "position matching requested but the genome build could not be verified (%s).\n",
        "  Positions from two different builds match the WRONG variants, and a\n",
        "  piecewise offset truncates the weights spatially, which is exactly the\n",
        "  failure a regional statistic like RAS cannot absorb.\n",
        "  Either supply IDs that overlap (so the build can be checked), or pass\n",
        "  build_check = FALSE after confirming both sides are on the same build."),
        bv$reason), call. = FALSE)
    }
    if (!bv$consistent) {
      stop(sprintf(paste0(
        "genome build mismatch: %d ID-matched SNPs have median position\n",
        "  difference %s bp (%.1f%% non-zero; quantiles %s).\n",
        "  %sRefusing to match on position -- lift the summary statistics onto the\n",
        "  map's build first (or translate rsIDs to the map's naming), then retry."),
        bv$n, format(bv$median_dbp, big.mark = ","), 100 * bv$frac_nonzero,
        paste(trimws(format(bv$quantiles, big.mark = ",")), collapse = "/"),
        if (isTRUE(bv$piecewise))
          "This is the PIECEWISE case: part of the chromosome agrees and part does not.\n  "
        else ""), call. = FALSE)
    }
  }

  # ---- 3. match by position ---------------------------------------------
  # Requires the unordered allele pair to agree as well, so multi-allelic sites
  # cannot be mismatched. NOTE the guard degrades to a bare chr:pos key when the
  # summary statistics carry no A2 column (PLINK1 .assoc.linear); that is
  # reported rather than hidden.
  n_by_pos <- 0L
  if ("pos" %in% match_by) {
    if (!have_map_pos || !have_ss_pos) {
      warning('match_by includes "pos" but positions are missing on one side; ',
              "no position matching performed", call. = FALSE)
    } else {
      use_alleles <- !all(is.na(ssdf$A2))
      need <- which(is.na(idx))
      if (length(need)) {
        if (!use_alleles)
          message("[match] summary statistics carry no A2 column: position matching ",
                  "falls back to a bare chr:pos key (no allele-pair guard)")
        ss_key <- paste(ssdf$CHR, ssdf$BP,
                        if (use_alleles) .allele_key(ssdf$A1, ssdf$A2) else "", sep = ":")
        mp_key <- paste(.norm_chr(map$CHR), as.integer(map$BP),
                        if (use_alleles) .allele_key(map$A1, map$A2) else "", sep = ":")
        hit <- match(mp_key[need], ss_key)
        idx[need] <- hit
        src[need[!is.na(hit)]] <- "pos"
        n_by_pos <- sum(!is.na(hit))
      }
    }
  }
  matched <- !is.na(idx)

  # ---- 4. allele alignment ----------------------------------------------
  w <- rep(0, nrow(map))
  status <- rep("absent_from_sumstats", nrow(map))
  sign_used <- rep(NA_real_, nrow(map))

  if (any(matched)) {
    mi <- which(matched); si <- idx[mi]
    sa1 <- ssdf$A1[si]
    la1 <- map$A1[mi];  la2 <- map$A2[mi]

    same <- sa1 == la1
    flip <- sa1 == la2
    ok   <- same | flip
    amb  <- .ambig(la1, la2)

    s <- ifelse(same, 1, ifelse(flip, -1, NA_real_))
    if (drop_ambiguous) s[which(amb)] <- NA_real_

    take <- !is.na(s)
    w[mi[take]] <- ssdf$BETA[si[take]] * s[take]
    sign_used[mi] <- s
    # NA-safe: an unresolvable allele used to fall through every ifelse() and
    # land as NA, where table() then dropped it from the report entirely.
    status[mi] <- ifelse(is.na(ok), "unresolved_allele",
                  ifelse(!ok, "allele_mismatch",
                  ifelse(!is.na(amb) & amb & drop_ambiguous, "strand_ambiguous",
                  ifelse(same, "aligned", "flipped"))))
    src[mi[!take]] <- NA_character_   # matched but unusable: no weight, no route
  }

  usable <- w != 0
  detail <- data.frame(SNP = map$SNP, status = status, matched_by = src,
                       sign = sign_used, weight = w, stringsAsFactors = FALSE)
  rep_ <- as.data.frame(table(status, useNA = "ifany"), stringsAsFactors = FALSE)
  names(rep_) <- c("status", "n_snps")
  rep_$pct <- round(100 * rep_$n_snps / nrow(map), 2)

  blocks <- .match_blocks(usable, if (have_map_pos) as.integer(map$BP) else NULL,
                          block_size = block_size)

  # ---- 5. optional allele-frequency comparison ---------------------------
  # Diagnostic only. Cohorts genuinely differ in frequency, so a disagreement is
  # information, not an error, and nothing is dropped on it. What it can catch is
  # a systematic orientation problem: if the aligned frequencies anti-correlate,
  # the effect allele was read from the wrong column.
  af <- NULL
  if (!is.null(target_af)) {
    if (length(target_af) != nrow(map))
      stop(sprintf("target_af has length %d but the map holds %d variants",
                   length(target_af), nrow(map)), call. = FALSE)
    if (all(is.na(ssdf$EAF))) {
      warning("target_af was supplied but the summary statistics carry no ",
              "effect-allele-frequency column; frequency check skipped",
              call. = FALSE)
    } else {
      ok_af <- matched & !is.na(sign_used) & is.finite(target_af)
      k <- which(ok_af)
      ext_af <- ssdf$EAF[idx[k]]
      # orient the external frequency onto map$A1
      ext_af <- ifelse(sign_used[k] < 0, 1 - ext_af, ext_af)
      good <- is.finite(ext_af)
      k <- k[good]; ext_af <- ext_af[good]
      d <- abs(ext_af - target_af[k])
      af <- list(n_compared = length(k),
                 n_discordant = sum(d > af_tol),
                 tol = af_tol,
                 cor = if (length(k) > 2L) stats::cor(ext_af, target_af[k]) else NA_real_,
                 median_abs_diff = if (length(k)) stats::median(d) else NA_real_)
      message(sprintf(
        "[freq] %d variants compared | correlation %.4f | median |dAF| %.4f | %d discordant (>%.2f, %.2f%%)",
        af$n_compared, af$cor, af$median_abs_diff, af$n_discordant, af_tol,
        100 * af$n_discordant / max(af$n_compared, 1L)))
      if (!is.na(af$cor) && af$cor < 0)
        warning("aligned allele frequencies ANTI-correlate with the target ",
                "cohort's (r = ", round(af$cor, 3), "). The effect allele may ",
                "have been read from the wrong column -- check `cols`.",
                call. = FALSE)
    }
  }

  status_counts <- stats::setNames(rep_$n_snps, rep_$status)
  qc <- .harmonisation_qc(
    n_ss_rows = n_raw, n_missing_beta = n_nonfinite, n_dup = n_dup,
    n_ss_usable = nrow(ssdf), n_map = nrow(map),
    n_by_id = n_by_id, n_by_pos = n_by_pos, n_matched = sum(matched),
    status_counts = status_counts, n_usable = sum(usable), af = af)

  message(sprintf(
"[sumstats] file rows %d | non-finite BETA dropped %d | redundant ID+allele rows dropped %d
[sumstats] multi-allelic IDs kept %d (resolved by their alleles, not by row order)
[sumstats] matched by ID %d (of which %d resolved on ID+alleles) | matched by chr:pos %d
[sumstats] local map %d SNPs -> usable weights %d (%.1f%%), zero weights %d (%.1f%%)",
    n_raw, n_nonfinite, n_dup, n_multi, n_by_id, n_by_id_allele, n_by_pos,
    nrow(map), sum(usable), 100*mean(usable), sum(!usable), 100*mean(!usable)))

  # ---- 6. coverage guards ------------------------------------------------
  # (a) overall: bigsnpr's match.min.prop, denominator = the smaller dataset.
  denom <- min(nrow(map), nrow(ssdf))
  if (denom == 0L)
    stop("the summary statistics table is empty after parsing/alignment -- nothing to match",
         call. = FALSE)
  prop  <- sum(usable) / denom
  if (prop < match_min_prop) {
    stop(sprintf(paste0(
      "only %.1f%% of the smaller dataset (%d of %d) received a usable weight,\n",
      "  below match_min_prop = %.0f%%. A low match rate is the standard symptom of\n",
      "  a wrong build or an ID-convention mismatch. Check the [build] line above,\n",
      "  or lower match_min_prop if this coverage is genuinely expected."),
      100*prop, sum(usable), denom, 100*match_min_prop), call. = FALSE)
  }
  # (b) spatial: RAS is a REGIONAL statistic, so where the misses sit matters as
  #     much as how many there are. Reported, never fatal.
  if (!is.null(blocks) && nrow(blocks) > 1) {
    wb <- which.min(blocks$rate); med_rate <- stats::median(blocks$rate)
    message(sprintf(
      "[coverage] per-%d-SNP block match rate: min %.1f%% (block %d, SNP %d-%d%s) | median %.1f%% | max %.1f%%",
      block_size, 100*blocks$rate[wb], blocks$block[wb],
      blocks$idx_start[wb], blocks$idx_end[wb],
      if (is.finite(blocks$bp_start[wb]))
        sprintf(", %.2f-%.2f Mb", blocks$bp_start[wb]/1e6, blocks$bp_end[wb]/1e6) else "",
      100*med_rate, 100*max(blocks$rate)))
    if (blocks$rate[wb] < 0.5 * med_rate)
      warning(sprintf(paste0(
        "coverage is SPATIALLY UNEVEN: block %d matches %.1f%% against a median of %.1f%%.\n",
        "  RAS cannot produce signal where the weights are zero, so an uneven map\n",
        "  biases WHERE the scan is able to detect, not just how strongly."),
        blocks$block[wb], 100*blocks$rate[wb], 100*med_rate), call. = FALSE)
  }

  list(weights = w, qc = qc, report = rep_, detail = detail,
       matched_by = c(id = n_by_id, id_allele = n_by_id_allele, pos = n_by_pos),
       build = bv, blocks = blocks, af = af)
}

# ===========================================================================
# The driver
# ===========================================================================
#' Harmonise External GWAS Summary Statistics into RAS Weights
#'
#' The entry point for summary-informed RAS. Takes a raw external
#' summary-statistics file and a local genotype map, chooses the naming/build
#' route from the data rather than from an assumption, and returns a verified
#' weight vector together with a reproducible QC table.
#'
#' The driver never guesses about the genome build. It inspects both sides' ID
#' conventions and then either matches directly, translates the summary
#' statistics through a dictionary, lifts their coordinates over, or refuses --
#' listing exactly what would unblock it. Every path ends at
#' \code{\link{ras_weights_from_sumstats}} with \code{build_check = TRUE} unless
#' the caller has explicitly asserted the builds agree or a liftOver round trip
#' has already established them.
#'
#' @param sumstats Path to the external summary statistics, or a data frame.
#' @param map Data frame with columns \code{SNP}, \code{CHR}, \code{BP},
#'   \code{A1}, \code{A2}: the local genotype map, one row per genotype column,
#'   in order. \code{A1} must be the allele the local dosage counts.
#' @param dict Optional variant dictionary ON THE MAP'S BUILD, from
#'   \code{\link{ras_read_variant_dictionary}}. Enables the dictionary route,
#'   which fixes a naming disagreement without touching any coordinate.
#' @param chain,chain_back Optional liftOver chain files, enabling the
#'   coordinate-conversion route. \code{chain_back} is required unless
#'   \code{check_reverse = FALSE}.
#' @param cup Optional data frame \code{(CHR, START, END)} of conversion-unstable
#'   intervals on the source build, passed to the liftOver route.
#' @param converter Conversion function for the liftOver route; the default
#'   shells out to the UCSC \code{liftOver} binary.
#' @param check_reverse Logical. Round-trip every liftOver conversion and drop
#'   what cannot return to its own position (default \code{TRUE}).
#' @param assume_same_build Logical. Assert, on the caller's authority, that
#'   both sides use the same genome build. Only then will positions be used
#'   without the build gate having verified them. Appropriate when the summary
#'   statistics are a GWAS Catalog harmonised (GRCh38) file and the map is
#'   GRCh38.
#' @param cols Optional column-name overrides, passed through.
#' @param ... Passed to \code{\link{ras_weights_from_sumstats}}, e.g.
#'   \code{drop_ambiguous}, \code{match_min_prop}, \code{target_af},
#'   \code{block_size}.
#'
#' @return A list with
#'   \item{weights}{Numeric vector of length \code{nrow(map)}, ready to pass to
#'     \code{\link{ras_scan_external}}.}
#'   \item{qc}{A data frame with one row per harmonisation stage
#'     (\code{item}, \code{n}, \code{denom}, \code{pct}): external variants
#'     read, dropped for a missing beta, dropped as duplicates, matched by rsID,
#'     matched by position, absent from the target, allele-aligned,
#'     allele-flipped, palindromic removed, and final usable weights. Written
#'     to be quoted directly in a manuscript.}
#'   \item{route}{The route actually taken.}
#'   \item{trail}{Human-readable audit trail of every decision.}
#'   \item{weights_result}{The full \code{\link{ras_weights_from_sumstats}}
#'     result, including \code{detail}, \code{build} and \code{blocks}.}
#'   \item{sumstats_aligned}{The summary statistics after any translation or
#'     conversion.}
#'
#' @seealso \code{\link{ras_scan_external}} to scan with these weights;
#'   \code{\link{ras_sumstats_report}} to inspect an alignment without
#'   committing to it.
#' @export
ras_harmonize_sumstats <- function(sumstats, map, dict = NULL,
                                   chain = NULL, chain_back = NULL, cup = NULL,
                                   converter = .liftover_converter_ucsc,
                                   check_reverse = TRUE,
                                   assume_same_build = FALSE,
                                   cols = list(), ...) {
  stopifnot(all(c("SNP","A1","A2") %in% names(map)))
  ss <- if (is.character(sumstats)) .read_delim_table(sumstats) else as.data.frame(sumstats)
  ss <- .keep_add_rows(ss)
  c_snp <- cols$snp %||% .pick(names(ss), "snp")
  if (is.na(c_snp)) stop("no SNP-ID column in the summary statistics", call. = FALSE)

  plan <- .plan_alignment(map$SNP, ss[[c_snp]])
  trail <- c(sprintf("plan: %s -> route '%s'", plan$why, plan$route), plan$note)
  message("[plan] ", plan$why)
  message("[plan] route = ", plan$route, "  (", plan$note, ")")

  route <- plan$route; aligned <- ss

  # --- naming ------------------------------------------------------------
  if (route == "dict") {
    if (is.null(dict)) {
      if (!is.null(chain)) {
        route <- "liftover"
      } else if (isTRUE(assume_same_build)) {
        route <- "pos_asserted"
      } else {
        stop(sprintf(paste0(
          "the two sides name variants differently (%s) and nothing was supplied to bridge them.\n",
          "  Give ONE of:\n",
          "    dict = ras_read_variant_dictionary(...)  # on the map's build -- preferred,\n",
          "                                             # a lookup cannot resolve to the wrong variant\n",
          "    chain = <hg19ToHg38.over.chain>          # only if the builds genuinely differ\n",
          "    assume_same_build = TRUE                 # only if you have confirmed both are the same build"),
          plan$why), call. = FALSE)
      }
    } else {
      tr <- .sumstats_to_map_ids(ss, dict, map, cols = cols)
      aligned <- tr$sumstats; route <- "dict"
      trail <- c(trail, sprintf("B1: %d of %d rows translated (%.1f%%); allele order %s",
                                tr$report$n_out, tr$report$n_in, tr$report$pct_out,
                                tr$allele_order))
    }
  }

  # --- coordinates -------------------------------------------------------
  if (route == "liftover") {
    lo <- ras_liftover_sumstats(aligned, chain = chain, chain_back = chain_back,
                                cols = cols, converter = converter,
                                check_reverse = check_reverse, cup = cup)
    aligned <- lo$sumstats
    trail <- c(trail, sprintf("B2: %s",
      paste(sprintf("%s %d", lo$report$status, lo$report$n), collapse = ", ")))
  }

  # --- weights, with the gate --------------------------------------------
  # Both keys, everywhere they are both available. IDs are build-invariant but
  # not unique and not universal; positions are unique but build-dependent. The
  # ID-matched subset is also what the build gate verifies positions on, so it
  # licenses the position matching that then reaches variants with no rsID.
  match_by <- switch(route,
    id           = c("id", "pos"),
    dict         = c("id", "pos"),
    liftover     = c("id", "pos"),
    pos_asserted = "pos",
    c("id", "pos"))
  build_check <- TRUE
  if (route == "pos_asserted" && isTRUE(assume_same_build)) {
    build_check <- FALSE
    trail <- c(trail, "build gate SKIPPED on the caller's assume_same_build = TRUE")
  }
  if (route == "liftover" && isTRUE(check_reverse)) {
    # The round trip already established the conversion: every surviving row
    # returned to its own source chromosome and position. Re-running the gate
    # would only ask whether IDs overlap, which after a liftOver they need not.
    build_check <- FALSE
    trail <- c(trail, "build verified by the liftOver round trip; gate not re-run")
  }

  res <- ras_weights_from_sumstats(aligned, map = map, cols = cols,
                                   match_by = match_by, build_check = build_check, ...)
  trail <- c(trail, sprintf("weights: %d usable of %d map variants (%.1f%%)",
                            sum(res$weights != 0), nrow(map),
                            100*mean(res$weights != 0)))
  message("[done] route '", route, "'")
  list(weights = res$weights, qc = res$qc, route = route, trail = trail,
       weights_result = res, sumstats_aligned = aligned)
}

# ===========================================================================
# Pre-flight inspection
# ===========================================================================
#' Inspect an External-Summary Alignment Without Committing to It
#'
#' A non-fatal report on how an external summary-statistics file lines up with a
#' local genotype map: which ID convention each side uses, the suggested route,
#' how many variants overlap by ID, what the build gate would say, and how the
#' overlap is distributed along the chromosome. Run this before
#' \code{\link{ras_harmonize_sumstats}} to see what the gates will do without
#' triggering them.
#'
#' The coverage line matters more here than in a PRS setting: RAS's statistic is
#' regional, so weights that are missing in one stretch of the chromosome remove
#' the scan's ability to detect there, rather than merely weakening it overall.
#'
#' @param ss Data frame of summary statistics, or a path to one.
#' @param map The local genotype map (\code{SNP}, and \code{CHR}/\code{BP} for
#'   the build and coverage lines).
#' @param cols Optional column-name overrides.
#' @param block_size Integer. Window, in variants, for the coverage profile.
#'
#' @return Invisibly, a list with \code{plan}, \code{n_id}, \code{build} and
#'   \code{blocks}. Called for the report it prints.
#'
#' @seealso \code{\link{ras_harmonize_sumstats}}
#' @export
ras_sumstats_report <- function(ss, map, cols = list(), block_size = 500L) {
  if (is.character(ss)) ss <- .read_delim_table(ss)
  ss <- .keep_add_rows(as.data.frame(ss))
  nm <- names(ss)
  c_snp <- cols$snp %||% .pick(nm, "snp"); c_pos <- cols$pos %||% .pick(nm, "pos")
  plan <- .plan_alignment(map$SNP, ss[[c_snp]])
  i <- match(as.character(map$SNP), as.character(ss[[c_snp]]))
  n_id <- sum(!is.na(i))
  bv <- if (n_id > 0L && !is.na(c_pos) && all(c("CHR","BP") %in% names(map)))
          .build_verdict(suppressWarnings(as.integer(ss[[c_pos]][i[!is.na(i)]])),
                         as.integer(map$BP[!is.na(i)]))
        else list(evaluable = FALSE, n = n_id, reason = "no ID overlap or no positions")
  blocks <- .match_blocks(!is.na(i), if ("BP" %in% names(map)) as.integer(map$BP) else NULL,
                          block_size = block_size)
  message(sprintf("ID convention : %s\n", plan$why))
  message(sprintf("suggested route: %s\n", plan$route))
  message(sprintf("ID overlap    : %d of %d map variants (%.1f%%)\n",
              n_id, nrow(map), 100*n_id/nrow(map)))
  if (isTRUE(bv$evaluable))
    message(sprintf("build         : median dBP %s | non-zero %.1f%% -> %s\n",
                format(bv$median_dbp, big.mark = ","), 100*bv$frac_nonzero,
                if (bv$consistent) "SAME BUILD" else if (isTRUE(bv$piecewise))
                  "PIECEWISE OFFSET" else "BUILD MISMATCH"))
  else message(sprintf("build         : not evaluable (%s)\n", bv$reason))
  if (!is.null(blocks) && nrow(blocks) > 1)
    message(sprintf("coverage      : min %.1f%% / median %.1f%% / max %.1f%% per %d-SNP block\n",
                100*min(blocks$rate), 100*stats::median(blocks$rate),
                100*max(blocks$rate), block_size))
  invisible(list(plan = plan, n_id = n_id, build = bv, blocks = blocks))
}

# ===========================================================================
# Split-free scanning
# ===========================================================================
#' RAS Scan with External Weights on All Samples
#'
#' Runs the Stage-1 forward scan with weights taken from an independent
#' external GWAS instead of a within-sample training split. Because the
#' weights are independent of the target cohort, no 50/50 split is drawn and
#' no repetitions are averaged: every sample is scanned once. Genotypes are
#' streamed from a \code{.rasbin} file in chunks. The in-memory implementation
#' is \code{\link{ras_scan_external_original}}.
#'
#' The thresholds of the changepoint detector were set on averaged profiles,
#' whose noise floor is lower than that of a single-pass profile; calibrate
#' before treating a detection from this scan as a finding (for the box-scan
#' detector, \code{\link{ras_box_calibrate}} on permuted phenotypes).
#'
#' @param geno Either the path to a \code{.rasbin} genotype file (see
#'   \code{\link{geno_to_rasbin}}, \code{\link{bed_to_rasbin}}) or an
#'   in-memory numeric genotype matrix, which is converted to a temporary
#'   \code{.rasbin} file for the run.
#' @param phenotype Numeric vector of length \eqn{n}, in the sample order of
#'   \code{geno}.
#' @param covariates Data frame with \eqn{n} rows, in the same sample order.
#' @param covariate_cols Character vector of covariate column names.
#' @param weights Numeric vector of length \eqn{N}, one weight per variant in
#'   the order of \code{geno}, typically
#'   \code{ras_harmonize_sumstats(...)$weights}. Non-finite entries count as
#'   zero.
#' @param is_continuous Logical. \code{TRUE} for a quantitative trait.
#' @inheritParams ras_scan
#' @param save_dir Character or \code{NULL}. Directory in which to save the
#'   profile as \code{ext_scan_chr<chrom>.rds}; \code{NULL} (default) saves
#'   nothing.
#'
#' @return A list with \code{x} (the SNP index of each grid point) and
#'   \code{y} (the \eqn{-\log_{10}(p)} profile).
#'
#' @seealso \code{\link{ras_harmonize_sumstats}} to build \code{weights};
#'   \code{\link{ras_scan}} for the split-based scan;
#'   \code{\link{ras_box_detect}} and \code{\link{ras_detect}} to detect
#'   regions on the profile.
#' @export
ras_scan_external <- function(geno, phenotype, covariates,
                                   covariate_cols, weights,
                                   is_continuous = TRUE,
                                   skip1 = 10, skip2 = 20,
                                   min_window_size = 5, max_window_size = 100,
                                   chunk_snps = 5000, chrom = 1,
                                   save_dir = NULL) {
  rb <- .ras_resolve_geno(geno)
  on.exit(rb$cleanup(), add = TRUE)
  rasbin_path <- rb$path
  hdr <- rasbin_header(rasbin_path)   # atomic named vector -- [[ ]], never $
  n <- as.integer(hdr[["n_samples"]])
  N <- as.integer(hdr[["n_snps"]])
  if (length(weights) != N)
    stop(sprintf("weights has length %d but the RASBIN holds %d variants",
                 length(weights), N), call. = FALSE)
  if (length(phenotype) != n || nrow(covariates) != n)
    stop(sprintf("phenotype/covariates must have %d rows to match the RASBIN", n),
         call. = FALSE)
  weights <- as.numeric(weights)
  weights[!is.finite(weights)] <- 0
  if (all(weights == 0)) stop("every weight is zero", call. = FALSE)

  scan.df <- covariates[, covariate_cols, drop = FALSE]
  scan.df$phenotype2 <- phenotype

  message(sprintf(
    "[scan] %s samples (ALL, no split) x %s variants | %s non-zero weights | chunk %d",
    format(n, big.mark = ","), format(N, big.mark = ","),
    format(sum(weights != 0), big.mark = ","), chunk_snps))

  y <- screen_forward_max_region(
    rasbin_path       = rasbin_path,
    weights           = weights,
    this.leftout      = seq_len(n),        # ALL samples: no split, no repetitions
    this.df           = scan.df,
    is_continuous     = is_continuous,
    covariate_formula = paste(covariate_cols, collapse = " + "),
    skip1             = skip1,
    skip2             = skip2,
    min_window_size   = min_window_size,
    max_window_size   = max_window_size,
    chunk_snps        = chunk_snps)

  out <- list(x = seq(1, N, by = skip1), y = as.numeric(y))
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    saveRDS(out, file.path(save_dir, sprintf("ext_scan_chr%s.rds", chrom)))
  }
  out
}

#' Summary-Informed RAS Scan on All Samples (In-Memory)
#'
#' The in-memory counterpart of \code{\link{ras_scan_external}}, for
#' cohorts small enough that a dense \eqn{n \times m} genotype matrix and its
#' PGS matrix fit in RAM. At biobank scale use
#' \code{\link{ras_scan_external}} instead: at 453,698 samples by 332,690
#' variants the dense double alone is over a petabyte.
#'
#' @param geno Numeric \eqn{n \times m} genotype dosage matrix.
#' @param phenotype Numeric vector of length \eqn{n}.
#' @param covariates Data frame with \eqn{n} rows containing
#'   \code{covariate_cols}.
#' @param covariate_cols Character vector of covariate column names.
#' @param weights Numeric vector of length \eqn{m}, aligned to the columns of
#'   \code{geno}.
#' @param is_continuous Logical. \code{TRUE} for quantitative traits.
#' @param skip1,skip2,min_window_size,max_window_size Scan geometry, as in
#'   \code{\link{ras_scan_original}}.
#' @param scan_test \code{"glm"} (per-window Wald) or \code{"score"} (Rao score
#'   test).
#' @param chrom Integer/character. Chromosome label used in output filenames.
#' @param save_dir Character or \code{NULL}. Directory to save the profile in.
#'
#' @return A list with \code{x} (grid point variant indices) and \code{y} (the
#'   profile).
#'
#' @seealso \code{\link{ras_scan_external}}, \code{\link{ras_harmonize_sumstats}}
#' @export
ras_scan_external_original <- function(geno, phenotype, covariates, covariate_cols,
                              weights, is_continuous = TRUE,
                              skip1 = 10, skip2 = 20,
                              min_window_size = 5, max_window_size = 100,
                              scan_test = c("glm","score"), chrom = 1,
                              save_dir = NULL) {
  scan_test <- match.arg(scan_test)
  stopifnot(length(weights) == ncol(geno))
  weights <- as.numeric(weights)
  weights[!is.finite(weights)] <- 0
  if (all(weights == 0)) stop("every weight is zero", call. = FALSE)
  cov.df <- covariates[, covariate_cols, drop = FALSE]

  pgs.mat <- compute_pgs_matrix(geno, seq_len(nrow(geno)), weights)
  scan.df <- cov.df; scan.df$phenotype2 <- phenotype

  y <- screen_forward_max_region_original(
    geno = geno, pgs.mat = pgs.mat, this.df = scan.df, num_signals = -1,
    is_continuous = is_continuous,
    covariate_formula = paste(covariate_cols, collapse = " + "),
    skip1 = skip1, skip2 = skip2, min_window_size = min_window_size,
    max_window_size = max_window_size, scan_test = scan_test, isPlot = FALSE)

  out <- list(x = seq(1, ncol(geno), by = skip1), y = y)
  if (!is.null(save_dir)) {
    if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
    saveRDS(out, file.path(save_dir, sprintf("ext_scan_chr%s.rds", chrom)))
  }
  out
}
