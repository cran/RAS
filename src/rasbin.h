/* .rasbin chunked genotype file format
 *
 * Header (32 bytes, fixed):
 *   bytes  0- 7  magic     "RASBIN01" (8 ASCII bytes, no NUL terminator)
 *   bytes  8-15  n_samples int64, little-endian
 *   bytes 16-23  n_snps    int64, little-endian
 *   bytes 24-27  dtype     int32 (0 = double, 1 = int8)
 *   bytes 28-31  reserved  (zero)
 *
 * Body, dtype 0 (double): n_samples * n_snps doubles, column-major. Column j
 * (0-based) occupies bytes [32 + j*n_samples*8, 32 + (j+1)*n_samples*8) so an
 * arbitrary column range is one contiguous seek+read.
 *
 * Body, dtype 1 (int8): n_samples * n_snps signed bytes, column-major, one
 * byte per genotype (values 0/1/2 = dosage, RASBIN_INT8_NA = missing). Column
 * j occupies bytes [32 + j*n_samples, 32 + (j+1)*n_samples) — same
 * fixed-stride seek pattern as dtype 0, just a 1-byte instead of 8-byte
 * element, so rasbin_fread_columns() below handles both without the caller
 * needing to know which one is on disk.
 *
 * NOTE (room for improvement, not yet done): dtype 1 spends a full byte per
 * genotype, i.e. it is ~4x larger on disk than a PLINK .bed of the same
 * dimensions (which packs 4 genotypes/byte at 2 bits each). A future dtype 2
 * could adopt that 2-bit packing for a further ~4x reduction, at the cost of
 * bit-level (non-byte-aligned) arithmetic in rasbin_fread_columns() for
 * partial-column reads. Deliberately deferred: byte-per-genotype keeps this
 * reader's fixed-stride-seek simplicity, and 4x was judged an acceptable
 * trade for that (see project discussion, 2026-07-23).
 */
#ifndef RAS_RASBIN_H
#define RAS_RASBIN_H

#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <R.h>
#include <Rinternals.h>

#define RASBIN_MAGIC "RASBIN01"
#define RASBIN_MAGIC_LEN 8
#define RASBIN_HEADER_SIZE 32
#define RASBIN_DTYPE_DOUBLE 0
#define RASBIN_DTYPE_INT8 1
#define RASBIN_INT8_NA ((int8_t) -1)

typedef struct {
    int64_t n_samples;
    int64_t n_snps;
    int32_t dtype;
} rasbin_header_t;

#ifdef _WIN32
  #define ras_fseek64(f, off, whence) _fseeki64((f), (off), (whence))
  #define ras_ftell64(f) _ftelli64((f))
#else
  #define ras_fseek64(f, off, whence) fseeko((f), (off), (whence))
  #define ras_ftell64(f) ftello((f))
#endif

/* Reads and validates the header from an already-open file positioned at 0.
 * Returns 0 on success, -1 on malformed/short header, -2 on bad magic. */
static inline int rasbin_read_header(FILE *fp, rasbin_header_t *hdr) {
    unsigned char buf[RASBIN_HEADER_SIZE];
    if (fread(buf, 1, RASBIN_HEADER_SIZE, fp) != RASBIN_HEADER_SIZE) return -1;
    if (memcmp(buf, RASBIN_MAGIC, RASBIN_MAGIC_LEN) != 0) return -2;
    memcpy(&hdr->n_samples, buf + 8,  8);
    memcpy(&hdr->n_snps,    buf + 16, 8);
    memcpy(&hdr->dtype,     buf + 24, 4);
    return 0;
}

/* Reads ncols full columns (col_start..col_start+ncols-1, 1-based, all
 * hdr->n_samples rows) from an already-open .rasbin file into `out`
 * (column-major, n_samples x ncols doubles). Transparently widens a
 * RASBIN_DTYPE_INT8 body to double (RASBIN_INT8_NA -> NA_REAL) via a small
 * heap scratch buffer, so every caller only ever deals with doubles
 * regardless of on-disk dtype. One seek + one read. Returns 0 on success,
 * -1 on seek failure, -2 on short read. */
static inline int rasbin_fread_columns(FILE *fp, const rasbin_header_t *hdr,
                                        int64_t col_start, int64_t ncols,
                                        double *out) {
    int64_t n_samples = hdr->n_samples;
    size_t elem_size = (hdr->dtype == RASBIN_DTYPE_INT8) ? sizeof(int8_t) : sizeof(double);
    int64_t offset = (int64_t) RASBIN_HEADER_SIZE +
                      (col_start - 1) * n_samples * (int64_t) elem_size;
    if (ras_fseek64(fp, offset, SEEK_SET) != 0) return -1;
    size_t nelem = (size_t) n_samples * (size_t) ncols;

    if (hdr->dtype == RASBIN_DTYPE_INT8) {
        int8_t *raw8 = (int8_t *) R_Calloc(nelem, int8_t);
        size_t nread = fread(raw8, sizeof(int8_t), nelem, fp);
        if (nread != nelem) { R_Free(raw8); return -2; }
        for (size_t i = 0; i < nelem; i++) {
            out[i] = (raw8[i] == RASBIN_INT8_NA) ? NA_REAL : (double) raw8[i];
        }
        R_Free(raw8);
    } else {
        size_t nread = fread(out, sizeof(double), nelem, fp);
        if (nread != nelem) return -2;
    }
    return 0;
}

#endif /* RAS_RASBIN_H */
