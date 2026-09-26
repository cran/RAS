#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

SEXP RAS_malloc_trim(void);
SEXP RAS_malloc_trim_available(void);
SEXP RAS_rasbin_header(SEXP path_sexp);
SEXP RAS_rasbin_read_chunk(SEXP path_sexp, SEXP col_start_sexp, SEXP col_end_sexp);
SEXP RAS_fused_pgs_scan(
    SEXP path_sexp, SEXP weights_sexp, SEXP holdout_idx_sexp, SEXP mode_sexp,
    SEXP keep_mask_sexp, SEXP Zmat_sexp, SEXP Mmat_sexp, SEXP resid_sexp,
    SEXP w0_sexp, SEXP syy_sexp, SEXP df_sexp,
    SEXP skip1_sexp, SEXP skip2_sexp, SEXP min_ws_sexp, SEXP max_ws_sexp,
    SEXP chunk_snps_sexp);
SEXP RAS_chunked_gwas_continuous(SEXP path_sexp, SEXP training_idx_sexp,
                                  SEXP y_sexp, SEXP chunk_snps_sexp);
SEXP RAS_chunked_gwas_binary_clean(SEXP path_sexp, SEXP train_ok_idx_sexp,
                                    SEXP Mmat_sexp, SEXP Wmat_sexp,
                                    SEXP ryb_sexp, SEXP Syy_sexp, SEXP df_sexp,
                                    SEXP chunk_snps_sexp);
SEXP RAS_detect_fast(
    SEXP x_sexp, SEXP y_sexp, SEXP p_thresh_sexp, SEXP min_length_sexp,
    SEXP skip_sexp, SEXP window_size_sexp, SEXP slope_win_sexp,
    SEXP thresh_left_sexp, SEXP thresh_right_sexp);
SEXP RAS_davies_test(SEXP x_sexp, SEXP y_sexp);
SEXP RAS_slope_test(SEXP x_sexp, SEXP y_sexp, SEXP lower_tail_sexp);
SEXP RAS_get_break_points(SEXP x_sexp, SEXP y_sexp, SEXP t_sexp);
SEXP RAS_seg_fit_boot(SEXP x_sexp, SEXP y_sexp, SEXP n_boot_ignored);
SEXP RAS_seg_fit_single(SEXP x_sexp, SEXP y_sexp);
SEXP RAS_bed_to_rasbin(SEXP bed_path_sexp, SEXP n_samples_sexp,
                        SEXP n_snps_sexp, SEXP out_path_sexp,
                        SEXP chunk_snps_sexp);

static const R_CallMethodDef CallEntries[] = {
    {"RAS_malloc_trim",             (DL_FUNC) &RAS_malloc_trim,             0},
    {"RAS_malloc_trim_available",   (DL_FUNC) &RAS_malloc_trim_available,   0},
    {"RAS_rasbin_header",           (DL_FUNC) &RAS_rasbin_header,           1},
    {"RAS_rasbin_read_chunk",       (DL_FUNC) &RAS_rasbin_read_chunk,       3},
    {"RAS_fused_pgs_scan",          (DL_FUNC) &RAS_fused_pgs_scan,          16},
    {"RAS_chunked_gwas_continuous", (DL_FUNC) &RAS_chunked_gwas_continuous, 4},
    {"RAS_chunked_gwas_binary_clean", (DL_FUNC) &RAS_chunked_gwas_binary_clean, 8},
    {"RAS_detect_fast",             (DL_FUNC) &RAS_detect_fast,             9},
    {"RAS_davies_test",             (DL_FUNC) &RAS_davies_test,             2},
    {"RAS_slope_test",              (DL_FUNC) &RAS_slope_test,              3},
    {"RAS_get_break_points",        (DL_FUNC) &RAS_get_break_points,        3},
    {"RAS_seg_fit_boot",            (DL_FUNC) &RAS_seg_fit_boot,            3},
    {"RAS_seg_fit_single",          (DL_FUNC) &RAS_seg_fit_single,          2},
    {"RAS_bed_to_rasbin",           (DL_FUNC) &RAS_bed_to_rasbin,           5},
    {NULL, NULL, 0}
};

void R_init_RAS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
