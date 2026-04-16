#ifndef COMPAT_TBB_SCALABLE_ALLOCATOR_H
#define COMPAT_TBB_SCALABLE_ALLOCATOR_H

#include <stdlib.h>

static inline void* scalable_malloc(size_t size) {
    return malloc(size);
}

static inline void scalable_free(void* ptr) {
    free(ptr);
}

#endif
