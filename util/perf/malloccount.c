// LD_PRELOAD shim: counts malloc/calloc/realloc calls, prints total to stderr at exit.
// Build: gcc -O2 -shared -fPIC -o malloccount.so malloccount.c (bench.py expects it next to itself)
#include <stdio.h>
#include <stdlib.h>
extern void *__libc_malloc(size_t), *__libc_calloc(size_t, size_t), *__libc_realloc(void *, size_t);
static unsigned long long n;
void *malloc(size_t s) { __atomic_add_fetch(&n, 1, __ATOMIC_RELAXED); return __libc_malloc(s); }
void *calloc(size_t a, size_t b) { __atomic_add_fetch(&n, 1, __ATOMIC_RELAXED); return __libc_calloc(a, b); }
void *realloc(void *p, size_t s) { __atomic_add_fetch(&n, 1, __ATOMIC_RELAXED); return __libc_realloc(p, s); }
__attribute__((destructor)) static void report(void) { fprintf(stderr, "MALLOC_COUNT %llu\n", n); }
