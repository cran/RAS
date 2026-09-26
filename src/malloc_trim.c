#include <R.h>
#include <Rinternals.h>

/* malloc_trim() is a glibc extension: it exists on Linux/glibc but not on
 * other libcs such as musl (Alpine Linux), where <malloc.h> declares no such
 * function and the build fails. Compile the call only where glibc is present;
 * everywhere else release_memory() gets NA, as it already does on Windows
 * and macOS. RAS_malloc_trim_available() reports which branch was built, so
 * the unit test can demand 0/1 where malloc_trim() exists and NA elsewhere. */
#if defined(__linux__) && defined(__GLIBC__)
#include <malloc.h>
#define RAS_HAVE_MALLOC_TRIM 1
#else
#define RAS_HAVE_MALLOC_TRIM 0
#endif

SEXP RAS_malloc_trim(void) {
#if RAS_HAVE_MALLOC_TRIM
    int out = malloc_trim(0);
    return ScalarInteger(out);
#else
    return ScalarInteger(NA_INTEGER);
#endif
}

SEXP RAS_malloc_trim_available(void) {
    return ScalarLogical(RAS_HAVE_MALLOC_TRIM);
}
