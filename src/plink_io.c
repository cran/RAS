#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdint.h>
#include "rasbin.h"

/* PLINK 1 .bed -> .rasbin (dtype 1, int8) chunked converter.
 *
 * .bed layout (SNP-major, the only mode written by any current PLINK):
 *   bytes 0-2   magic 0x6c 0x1b 0x01
 *   body        one ceil(n_samples/4)-byte block per SNP, 4 genotypes/byte,
 *               2 bits/genotype, sample i's code at bits [2*(i%4), 2*(i%4)+1]
 *               of byte i/4 within its SNP's block; unused high bits in the
 *               last byte of a block are 0.
 *   2-bit code: 00 = hom A1, 01 = missing, 10 = het, 11 = hom A2 (A1/A2 per
 *   the companion .bim). Dosage below counts copies of A1 (PLINK's own
 *   --recode A / plink2 --export A convention), so it round-trips against
 *   other PLINK-adjacent tooling if anyone compares effect signs.
 *
 * Because every SNP occupies a fixed byte span, an arbitrary SNP range is one
 * seek + one contiguous read from the .bed, exactly like rasbin_fread_columns
 * already does for .rasbin itself -- so this converter never holds more than
 * one chunk's worth of genotypes in memory, regardless of total file size.
 */

static const int8_t PLINK_CODE_TO_DOSAGE[4] = { 2, RASBIN_INT8_NA, 1, 0 };

/* RAS_bed_to_rasbin(bed_path, n_samples, n_snps, out_path, chunk_snps)
 * n_samples/n_snps are supplied by the R wrapper (parsed from .fam/.bim)
 * rather than recounted here. Returns NULL invisibly; errors on any
 * malformed input. */
SEXP RAS_bed_to_rasbin(SEXP bed_path_sexp, SEXP n_samples_sexp,
                        SEXP n_snps_sexp, SEXP out_path_sexp,
                        SEXP chunk_snps_sexp) {
    const char *bed_path = CHAR(STRING_ELT(bed_path_sexp, 0));
    const char *out_path = CHAR(STRING_ELT(out_path_sexp, 0));
    int64_t n_samples = (int64_t) asReal(n_samples_sexp);
    int64_t n_snps    = (int64_t) asReal(n_snps_sexp);
    int64_t chunk_snps = (int64_t) asReal(chunk_snps_sexp);
    if (n_samples < 1) error("RAS_bed_to_rasbin: n_samples must be >= 1");
    if (n_snps < 1) error("RAS_bed_to_rasbin: n_snps must be >= 1");
    if (chunk_snps < 1) chunk_snps = 1;

    FILE *bed_fp = fopen(bed_path, "rb");
    if (!bed_fp) error("RAS_bed_to_rasbin: cannot open '%s'", bed_path);

    unsigned char magic[3];
    if (fread(magic, 1, 3, bed_fp) != 3) {
        fclose(bed_fp); error("RAS_bed_to_rasbin: '%s' too short to contain a .bed magic header", bed_path);
    }
    if (magic[0] != 0x6c || magic[1] != 0x1b) {
        fclose(bed_fp); error("RAS_bed_to_rasbin: '%s' is not a PLINK .bed file (bad magic)", bed_path);
    }
    if (magic[2] != 0x01) {
        fclose(bed_fp);
        error("RAS_bed_to_rasbin: '%s' is in individual-major .bed mode, which is not supported "
              "(no current PLINK version writes this; re-export in SNP-major mode)", bed_path);
    }

    int64_t bytes_per_snp = (n_samples + 3) / 4;
    int64_t expected_size = 3 + n_snps * bytes_per_snp;
    if (ras_fseek64(bed_fp, 0, SEEK_END) != 0) { fclose(bed_fp); error("RAS_bed_to_rasbin: seek-to-end failed on '%s'", bed_path); }
    int64_t actual_size = ras_ftell64(bed_fp);
    if (actual_size != expected_size) {
        fclose(bed_fp);
        error("RAS_bed_to_rasbin: '%s' size (%lld bytes) does not match n_samples=%lld x n_snps=%lld "
              "(expected %lld bytes) -- check .fam/.bim counts match this .bed",
              bed_path, (long long) actual_size, (long long) n_samples, (long long) n_snps,
              (long long) expected_size);
    }

    FILE *out_fp = fopen(out_path, "wb");
    if (!out_fp) { fclose(bed_fp); error("RAS_bed_to_rasbin: cannot open '%s' for writing", out_path); }

    /* .rasbin header, dtype 1 (int8). */
    {
        unsigned char hdrbuf[RASBIN_HEADER_SIZE];
        memset(hdrbuf, 0, RASBIN_HEADER_SIZE);
        memcpy(hdrbuf, RASBIN_MAGIC, RASBIN_MAGIC_LEN);
        memcpy(hdrbuf + 8,  &n_samples, 8);
        memcpy(hdrbuf + 16, &n_snps,    8);
        int32_t dtype = RASBIN_DTYPE_INT8;
        memcpy(hdrbuf + 24, &dtype, 4);
        if (fwrite(hdrbuf, 1, RASBIN_HEADER_SIZE, out_fp) != RASBIN_HEADER_SIZE) {
            fclose(bed_fp); fclose(out_fp);
            error("RAS_bed_to_rasbin: failed writing header to '%s'", out_path);
        }
    }

    unsigned char *packed = (unsigned char *) R_Calloc((size_t) (chunk_snps * bytes_per_snp), unsigned char);
    int8_t *unpacked = (int8_t *) R_Calloc((size_t) (chunk_snps * n_samples), int8_t);

    for (int64_t cs = 1; cs <= n_snps; cs += chunk_snps) {
        int64_t ce = cs + chunk_snps - 1; if (ce > n_snps) ce = n_snps;
        int64_t ncols = ce - cs + 1;
        int64_t chunk_bed_bytes = ncols * bytes_per_snp;

        int64_t bed_offset = 3 + (cs - 1) * bytes_per_snp;
        if (ras_fseek64(bed_fp, bed_offset, SEEK_SET) != 0) {
            R_Free(packed); R_Free(unpacked); fclose(bed_fp); fclose(out_fp);
            error("RAS_bed_to_rasbin: seek failed at SNP %lld", (long long) cs);
        }
        if (fread(packed, 1, (size_t) chunk_bed_bytes, bed_fp) != (size_t) chunk_bed_bytes) {
            R_Free(packed); R_Free(unpacked); fclose(bed_fp); fclose(out_fp);
            error("RAS_bed_to_rasbin: short read from '%s' at SNP %lld", bed_path, (long long) cs);
        }

        for (int64_t k = 0; k < ncols; k++) {
            const unsigned char *col_bytes = packed + k * bytes_per_snp;
            int8_t *col_out = unpacked + k * n_samples;
            for (int64_t i = 0; i < n_samples; i++) {
                unsigned char byte = col_bytes[i >> 2];
                int shift = (int) ((i & 3) * 2);
                int code = (byte >> shift) & 0x3;
                col_out[i] = PLINK_CODE_TO_DOSAGE[code];
            }
        }

        size_t nelem = (size_t) (ncols * n_samples);
        if (fwrite(unpacked, sizeof(int8_t), nelem, out_fp) != nelem) {
            R_Free(packed); R_Free(unpacked); fclose(bed_fp); fclose(out_fp);
            error("RAS_bed_to_rasbin: failed writing chunk [%lld,%lld] to '%s'",
                  (long long) cs, (long long) ce, out_path);
        }
    }

    R_Free(packed);
    R_Free(unpacked);
    fclose(bed_fp);
    fclose(out_fp);
    return R_NilValue;
}
