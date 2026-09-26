#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include "rasbin.h"

/* Read the .rasbin header. Returns integer(3) = c(n_samples, n_snps, dtype). */
SEXP RAS_rasbin_header(SEXP path_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    FILE *fp = fopen(path, "rb");
    if (!fp) error("rasbin_header: cannot open file '%s'", path);

    rasbin_header_t hdr;
    int rc = rasbin_read_header(fp, &hdr);
    fclose(fp);
    if (rc == -1) error("rasbin_header: file too short to contain a header: '%s'", path);
    if (rc == -2) error("rasbin_header: bad magic (not a .rasbin file): '%s'", path);

    SEXP out = PROTECT(allocVector(REALSXP, 3));
    REAL(out)[0] = (double) hdr.n_samples;
    REAL(out)[1] = (double) hdr.n_snps;
    REAL(out)[2] = (double) hdr.dtype;
    UNPROTECT(1);
    return out;
}

/* Read columns [col_start, col_end] (1-based, inclusive) from a .rasbin file.
 * Returns an n_samples x (col_end - col_start + 1) numeric matrix. */
SEXP RAS_rasbin_read_chunk(SEXP path_sexp, SEXP col_start_sexp, SEXP col_end_sexp) {
    const char *path = CHAR(STRING_ELT(path_sexp, 0));
    int64_t col_start = (int64_t) asReal(col_start_sexp);   /* 1-based, inclusive */
    int64_t col_end   = (int64_t) asReal(col_end_sexp);     /* 1-based, inclusive */

    FILE *fp = fopen(path, "rb");
    if (!fp) error("rasbin_read_chunk: cannot open file '%s'", path);

    rasbin_header_t hdr;
    int rc = rasbin_read_header(fp, &hdr);
    if (rc == -1) { fclose(fp); error("rasbin_read_chunk: file too short to contain a header: '%s'", path); }
    if (rc == -2) { fclose(fp); error("rasbin_read_chunk: bad magic (not a .rasbin file): '%s'", path); }
    if (hdr.dtype != RASBIN_DTYPE_DOUBLE && hdr.dtype != RASBIN_DTYPE_INT8) {
        fclose(fp); error("rasbin_read_chunk: unsupported dtype %d", hdr.dtype);
    }

    if (col_start < 1 || col_end < col_start || col_end > hdr.n_snps) {
        fclose(fp);
        error("rasbin_read_chunk: column range [%lld, %lld] out of bounds for n_snps=%lld",
              (long long) col_start, (long long) col_end, (long long) hdr.n_snps);
    }

    int64_t n = hdr.n_samples;
    int64_t ncols = col_end - col_start + 1;

    SEXP out = PROTECT(allocMatrix(REALSXP, (int) n, (int) ncols));
    int rc2 = rasbin_fread_columns(fp, &hdr, col_start, ncols, REAL(out));
    fclose(fp);

    if (rc2 == -1) {
        UNPROTECT(1);
        error("rasbin_read_chunk: seek failed for column range [%lld, %lld]",
              (long long) col_start, (long long) col_end);
    }
    if (rc2 == -2) {
        UNPROTECT(1);
        error("rasbin_read_chunk: short read for column range [%lld, %lld]",
              (long long) col_start, (long long) col_end);
    }

    UNPROTECT(1);
    return out;
}
