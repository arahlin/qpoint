#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "qpoint.h"

/*
  FNV-1a hash
  Available: http://ctips.pbworks.com/w/page/7277591/FNV%20Hash
 */
#define FNV_PRIME_32 16777619
#define FNV_OFFSET_32 2166136261U
uint32_t FNV32(const char *s, size_t len) {
  uint32_t hash = FNV_OFFSET_32;
  for(size_t i = 0; i < len; i++) {
    hash = hash ^ (s[i]); // xor next byte into the bottom of the hash
    hash = hash * FNV_PRIME_32; // Multiply by prime number found to work well
  }
  return hash;
}

long get_hash(long key, size_t n) {
  return FNV32((const char *) &key, sizeof(key)) % n;
  // return SuperFastHash((const char *) &key, sizeof(key)) % n;
}

qp_pixhash_t * qp_init_pixhash(long *pix, size_t npix) {
  qp_pixhash_t *pixhash = malloc(sizeof(*pixhash));

  pixhash->count = npix;
  pixhash->buckets = malloc(npix * sizeof(qp_pix_bucket_t));
  memset(pixhash->buckets, 0, npix * sizeof(qp_pix_bucket_t));

  qp_pix_bucket_t *bucket;
  qp_pix_pair_t *pair;
  long index;
  size_t count;

  for (size_t ii = 0; ii < npix; ii++) {
    index = get_hash(pix[ii], npix);
    bucket = pixhash->buckets + index;
    count = bucket->count;

    if (count == 0) {
      bucket->pairs = malloc(sizeof(qp_pix_pair_t));
    } else {
      bucket->pairs = realloc(bucket->pairs, (count + 1) * sizeof(qp_pix_pair_t));
    }
    bucket->count++;

    pair = bucket->pairs + count;
    pair->key = pix[ii];
    pair->index = ii;
  }

  pixhash->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;

#ifdef DEBUG
  for (size_t ii = 0; ii < npix; ii++) {
    printf("idx %12d | pix %12ld | repix %12ld | collisions %2ld\n",
           ii, pix[ii], qp_repixelize(pixhash, pix[ii]),
           pixhash->buckets[ii].count);
  }
#endif

  return pixhash;
}

qp_pixhash_t * qp_copy_pixhash(qp_pixhash_t *pixhash) {
  qp_pixhash_t *new_hash = malloc(sizeof(new_hash));
  qp_pix_bucket_t *bucket, *new_bucket;

  new_hash->count = pixhash->count;
  new_hash->buckets = malloc(new_hash->count * sizeof(qp_pix_bucket_t));
  memset(new_hash->buckets, 0, new_hash->count * sizeof(qp_pix_bucket_t));

  for (size_t ii = 0; ii < new_hash->count; ii++) {
    bucket = pixhash->buckets + ii;

    if (bucket->count == 0)
      continue;

    new_bucket = new_hash->buckets + ii;
    new_bucket->count = bucket->count;
    new_bucket->pairs = malloc(new_bucket->count * sizeof(qp_pix_pair_t));
    memcpy(new_bucket->pairs, bucket->pairs,
           new_bucket->count * sizeof(qp_pix_pair_t));
  }

  new_hash->init = QP_STRUCT_INIT | QP_STRUCT_MALLOC;
  return new_hash;
}

void qp_free_pixhash(qp_pixhash_t *pixhash) {
  qp_pix_bucket_t * bucket;
  if (pixhash->init & QP_STRUCT_MALLOC) {
    for (size_t ii = 0; ii < pixhash->count; ii++) {
      bucket = pixhash->buckets + ii;
      if (bucket->count == 0)
        continue;
      free(bucket->pairs);
    }
    free(pixhash->buckets);
    free(pixhash);
  }
}

long qp_repixelize(qp_pixhash_t *pixhash, long pix) {
  long index = get_hash(pix, pixhash->count);
  qp_pix_bucket_t *bucket = pixhash->buckets + index;

  if (bucket->count == 0)
    return -1;

  for (size_t ii = 0; ii < bucket->count; ii++) {
    if (bucket->pairs[ii].key == pix)
      return bucket->pairs[ii].index;
  }

  return -1;
}
