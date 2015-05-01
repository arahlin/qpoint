#include <stdio.h>
#include <stdlib.h>
#include "qpoint.h"

int qp_get_error_code(qp_memory_t *mem) {
  return mem->error_code;
}

char * qp_get_error_string(qp_memory_t *mem) {
  return mem->error_string;
}

void qp_set_error(qp_memory_t *mem, int error_code,
                  const char *error_string) {
  mem->error_code = error_code;
  mem->error_string = (char *)error_string;
}

// set error if condition is non-zero
int qp_check_error(qp_memory_t *mem, int condition, int error_code,
                   const char *error_string) {
  if (condition) {
    qp_set_error(mem, error_code, error_string);
    return error_code;
  }
  return 0;
}
