#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define epsilon (10e-6)

typedef struct simplex_t {
  int m;
  int n;
  int *var;
  double **a;
  double *b;
  double *x;
  double *c;
  double y;
} simplex_t;

typedef struct node_t {
  int m;
  int n;
  int k;
  int h;
  double xh;
  double ak;
  double bk;
  double *min;
  double *max;
  double **a;
  double *b;
  double *x;
  double *c;
  double z;
} node_t;

typedef struct node_set_t {
  int count;
  int capacity;
  node_t **nodes;
} node_set_t;

node_set_t *init_node_set(void);
void add_node_to_set(node_t *, node_set_t *);
int size(node_set_t *);
node_t *pop_node_from_set(node_set_t *);
void free_set(node_set_t *);
void free_node(node_t *);

node_t *initial_node(int, int, double **, double *, double *);
node_t *extend(node_t *, int, int, double **, double *, double *, int, double,
               double);
bool is_integer(double *);
bool integer(node_t *);
void bound(node_t *, node_set_t *, double *, double *);
int branch(node_t *, double);
void succ(node_t *, node_set_t *, int, int, double **, double *, double *, int,
          double, double, double *, double *);
double intopt(int, int, double **, double *, double *, double *);

double simplex(int, int, double **, double *, double *, double *, double);
double xsimplex(int, int, double **, double *, double *, double *, double,
                int *, int);
void pivot(simplex_t *, int, int);
int initial(simplex_t *, int, int, double **, double *, double *, double *,
            double, int *);
void prepare(simplex_t *, int);
int select_nonbasic(simplex_t *);
int init(simplex_t *, int, int, double **, double *, double *, double *, double,
         int *);

node_set_t *init_node_set(void) {
  node_set_t *res = malloc(sizeof(node_set_t));
  res->capacity = 20;
  res->count = 0;
  res->nodes = calloc(res->capacity, sizeof(node_t *));

  for (int i = 0; i < res->capacity; i++) {
    res->nodes[i] = NULL;
  }

  return res;
}

void add_node_to_set(node_t *node, node_set_t *set) {
  int i;

  if (set->count < set->capacity) {
    for (i = 0; i < set->capacity; i++) {
      if ((set->nodes)[i] == NULL) {
        set->count++;
        set->nodes[i] = node;
        return;
      }
    }
  } else {
    set->capacity = set->capacity * 2;
    set->nodes = realloc(set->nodes, set->capacity * sizeof(node_t *));
    for (i = set->count; i < set->capacity; i++) {
      set->nodes[i] = NULL;
    }
    set->nodes[set->count] = node;
    set->count++;
  }
}

int size(node_set_t *set) { return set->count; }

node_t *pop_node_from_set(node_set_t *set) {
  node_t *node;
  for (int i = 0; i < set->capacity; i++) {
    if ((set->nodes)[i]) {
      set->count--;
      node = set->nodes[i];
      set->nodes[i] = NULL;
      break;
    }
  }

  return node;
}

void free_set(node_set_t *set) {
  free(set->nodes);
  free(set);
}

void free_node(node_t *node) {
  for (int i = 0; i < node->m + 1; i++) {
    free(node->a[i]);
  }
  free(node->a);
  free(node->b);
  free(node->c);
  free(node->x);
  free(node->min);
  free(node->max);
  free(node);
}

node_t *initial_node(int m, int n, double **a, double *b, double *c) {
  int i;

  node_t *p = malloc(sizeof(node_t));
  p->a = calloc(m + 1, sizeof(double *));
  for (i = 0; i < m + 1; i++) {
    p->a[i] = calloc(n + 1, sizeof(double));
  }
  p->b = calloc(m + 1, sizeof(double));
  p->c = calloc(n + 1, sizeof(double));
  p->x = calloc(n + 1, sizeof(double));
  p->min = calloc(n, sizeof(double));
  p->max = calloc(n, sizeof(double));
  p->m = m;
  p->n = n;

  for (i = 0; i < m; i++) {
    memcpy(p->a[i], a[i], n * sizeof(double));
  }
  memcpy(p->b, b, m * sizeof(double));
  memcpy(p->c, c, n * sizeof(double));

  for (i = 0; i < n; i++) {
    p->min[i] = -INFINITY;
    p->max[i] = +INFINITY;
  }

  return p;
}

node_t *extend(node_t *node, int m, int n, double **a, double *b, double *c,
               int k, double ak, double bk) {
  node_t *res = malloc(sizeof(node_t));
  int i, j;

  res->k = k;
  res->ak = ak;
  res->bk = bk;

  if (ak > 0 && node->max[k] < INFINITY) {
    res->m = node->m;
  } else if (ak < 0 && node->min[k] > 0) {
    res->m = node->m;
  } else {
    res->m = node->m + 1;
  }

  res->n = node->n;
  res->h = -1;

  res->a = calloc(res->m + 1, sizeof(double *));
  for (i = 0; i < res->m + 1; i++) {
    res->a[i] = calloc(res->n + 1, sizeof(double));
  }
  res->b = calloc(res->m + 1, sizeof(double));
  res->c = calloc(res->n + 1, sizeof(double));
  res->x = calloc(res->n + 1, sizeof(double));
  res->max = calloc(n, sizeof(double));
  res->min = calloc(n, sizeof(double));

  memcpy(res->max, node->max, n * sizeof(double));
  memcpy(res->min, node->min, n * sizeof(double));
  for (i = 0; i < m; i++) {
    memcpy(res->a[i], a[i], n * sizeof(double));
  }
  memcpy(res->b, b, m * sizeof(double));
  memcpy(res->c, c, n * sizeof(double));

  if (ak > 0) {
    if (res->max[k] == INFINITY || bk < res->max[k]) {
      res->max[k] = bk;
    }
  } else if (res->min[k] == -INFINITY || -bk > res->min[k]) {
    res->min[k] = -bk;
  }

  for (i = m, j = 0; j < n; j++) {
    if (res->min[j] > -INFINITY) {
      res->a[i][j] = -1;
      res->b[i] = -res->min[j];
      i++;
    }
    if (res->max[j] < INFINITY) {
      res->a[i][j] = 1;
      res->b[i] = res->max[j];
      i++;
    }
  }
  return res;
}

bool is_integer(double *xp) {
  double x = *xp;
  double r = round(x);

  if (fabs(r - x) < epsilon) {
    *xp = r;
    return true;
  } else {
    return false;
  }
}

bool integer(node_t *node) {
  int i;

  for (i = 0; i < node->n; i++) {
    if (!is_integer(&(node->x[i]))) {
      return false;
    }
  }
  return true;
}

void bound(node_t *node, node_set_t *set, double *zp, double *x) {
  if (node->z > *zp) {
    *zp = node->z;
    memcpy(x, node->x, (node->n + 1) * sizeof(double));

    for (int i = 0; i < set->capacity; i++) {
      if (!set->nodes[i] || set->nodes[i]->z >= node->z) {
        continue;
      }

      free_node(set->nodes[i]);
      set->count--;
      set->nodes[i] = NULL;
    }
  }
}

int branch(node_t *node, double z) {
  double min, max;
  int h;

  if (node->z < z) {
    return 0;
  }

  for (h = 0; h < node->n; h++) {
    if (!is_integer(&(node->x[h]))) {
      if (node->min[h] == -INFINITY) {
        min = 0;
      } else {
        min = node->min[h];
      }

      max = node->max[h];

      if (floor(node->x[h]) < min || ceil(node->x[h]) > max) {
        continue;
      }

      node->h = h;
      node->xh = node->x[h];

      return 1;
    }
  }

  return 0;
}

int initial(simplex_t *s, int m, int n, double **a, double *b, double *c,
            double *x, double y, int *var) {
  int i, j, k;
  double w;

  k = init(s, m, n, a, b, c, x, y, var);

  if (b[k] >= 0) {
    return 1;
  }

  prepare(s, k);
  n = s->n;
  s->y = xsimplex(m, n, s->a, s->b, s->c, s->x, 0, s->var, 1);

  for (i = 0; i < m + n; i++) {
    if (s->var[i] == m + n - 1) {
      if (fabs(s->x[i]) > epsilon) {
        free(s->c);
        free(s->x);
        return 0;
      } else {
        break;
      }
    }
  }

  if (i >= n) {
    for (j = k = 0; k < n; k++) {
      if (fabs(s->a[i - n][k]) > fabs(s->a[i - n][j])) {
        j = k;
      }
    }
    pivot(s, i - n, j);
    i = j;
  }

  if (i < n - 1) {
    k = s->var[i];
    s->var[i] = s->var[n - 1];
    s->var[n - 1] = k;
    for (k = 0; k < m; k++) {
      w = s->a[k][n - 1];
      s->a[k][n - 1] = s->a[k][i];
      s->a[k][i] = w;
    }
  }

  free(s->c);
  s->c = c;
  s->y = y;

  for (k = n - 1; k < n + m - 1; k++) {
    s->var[k] = s->var[k + 1];
  }

  n = s->n = s->n - 1;
  double *t = calloc(n, sizeof(double));

  bool go_to_next_k;
  for (k = 0; k < n; k++) {
    go_to_next_k = false;
    for (j = 0; j < n; j++) {
      if (k == s->var[j]) {
        t[j] = t[j] + s->c[k];
        go_to_next_k = true;
        break;
      }
    }

    if (go_to_next_k)
      continue;

    for (j = 0; j < m; j++) {
      if (s->var[n + j] == k) {
        break;
      }
    }

    s->y = s->y + s->c[k] * s->b[j];

    for (i = 0; i < n; i++) {
      t[i] = t[i] - s->c[k] * s->a[j][i];
    }
  }

  for (i = 0; i < n; i++) {
    s->c[i] = t[i];
  }

  free(t);
  free(s->x);

  return 1;
}

void prepare(simplex_t *s, int k) {
  int m = s->m;
  int n = s->n;
  int i;

  for (i = m + n; i > n; i--) {
    s->var[i] = s->var[i - 1];
  }

  s->var[n] = m + n;

  n = n + 1;
  for (i = 0; i < m; i++) {
    s->a[i][n - 1] = -1;
  }

  s->x = calloc(m + n, sizeof(double));
  s->c = calloc(n, sizeof(double));
  s->c[n - 1] = -1;
  s->n = n;
  pivot(s, k, n - 1);
}

int select_nonbasic(simplex_t *s) {
  int i;
  for (i = 0; i < s->n; i++) {
    if (s->c[i] > epsilon) {
      return i;
    }
  }
  return -1;
}

int init(simplex_t *s, int m, int n, double **a, double *b, double *c,
         double *x, double y, int *var) {
  int i, k;

  s->m = m;
  s->n = n;
  s->a = a;
  s->b = b;
  s->c = c;
  s->x = x;
  s->y = y;
  s->var = var;

  if (s->var == NULL) {
    s->var = calloc(m + n + 1, sizeof(int));
    for (i = 0; i < m + n; i++) {
      s->var[i] = i;
    }
  }

  for (k = 0, i = 1; i < m; i++) {
    if (b[i] < b[k]) {
      k = i;
    }
  }

  return k;
}

void succ(node_t *node, node_set_t *set, int m, int n, double **a, double *b,
          double *c, int k, double ak, double bk, double *zp, double *x) {
  node_t *extended = extend(node, m, n, a, b, c, k, ak, bk);

  if (extended == NULL) {
    return;
  }

  extended->z = simplex(extended->m, extended->n, extended->a, extended->b,
                        extended->c, extended->x, 0);

  if (isfinite(extended->z)) {
    if (integer(extended)) {
      bound(extended, set, zp, x);
    } else if (branch(extended, *zp)) {
      add_node_to_set(extended, set);
      return;
    }
  }

  free_node(extended);
}

void pivot(simplex_t *s, int row, int col) {
  double **a = s->a;
  double *b = s->b;
  double *c = s->c;
  int m = s->m;
  int n = s->n;
  int i, j, t;

  t = s->var[col];
  s->var[col] = s->var[n + row];
  s->var[n + row] = t;
  s->y = s->y + c[col] * b[row] / a[row][col];

  for (i = 0; i < n; i++) {
    if (i != col) {
      c[i] = c[i] - c[col] * a[row][i] / a[row][col];
    }
  }
  c[col] = -c[col] / a[row][col];

  for (i = 0; i < m; i++) {
    if (i != row) {
      b[i] = b[i] - a[i][col] * b[row] / a[row][col];
    }
  }

  for (i = 0; i < m; i++) {
    if (i != row) {
      for (j = 0; j < n; j++) {
        if (j != col) {
          a[i][j] = a[i][j] - a[i][col] * a[row][j] / a[row][col];
        }
      }
    }
  }

  for (i = 0; i < m; i++) {
    if (i != row) {
      a[i][col] = -a[i][col] / a[row][col];
    }
  }

  for (i = 0; i < n; i++) {
    if (i != col) {
      a[row][i] = a[row][i] / a[row][col];
    }
  }

  b[row] = b[row] / a[row][col];
  a[row][col] = 1 / a[row][col];
}

double xsimplex(int m, int n, double **a, double *b, double *c, double *x,
                double y, int *var, int h) {
  simplex_t *s = malloc(sizeof(simplex_t));
  int i, row, col;

  if (!initial(s, m, n, a, b, c, x, y, var)) {
    free(s->var);
    free(s);
    return NAN;
  }

  while ((col = select_nonbasic(s)) >= 0) {
    row = -1;
    for (i = 0; i < m; i++) {
      if (a[i][col] > epsilon &&
          (row < 0 || b[i] / a[i][col] < b[row] / a[row][col])) {
        row = i;
      }
    }

    if (row < 0) {
      free(s->var);
      free(s);
      return INFINITY;
    }

    pivot(s, row, col);
  }

  if (h == 0) {
    for (i = 0; i < n; i++) {
      if (s->var[i] < n) {
        x[s->var[i]] = 0;
      }
    }
    for (i = 0; i < m; i++) {
      if (s->var[n + i] < n) {
        x[s->var[n + i]] = s->b[i];
      }
    }
    free(s->var);
  } else {
    for (i = 0; i < n; i++) {
      x[i] = 0;
    }
    for (i = n; i < n + m; i++) {
      x[i] = s->b[i - n];
    }
  }

  double result = s->y;
  free(s);
  return result;
}

double simplex(int m, int n, double **a, double *b, double *c, double *x,
               double y) {
  return xsimplex(m, n, a, b, c, x, y, NULL, 0);
}

double intopt(int m, int n, double **a, double *b, double *c, double *x) {
  node_t *node = initial_node(m, n, a, b, c);
  node_set_t *h = init_node_set();
  add_node_to_set(node, h);

  double z = -INFINITY;
  node->z = simplex(node->m, node->n, node->a, node->b, node->c, node->x, 0);

  if (integer(node) || !isfinite(node->z)) {
    z = node->z;
    if (integer(node)) {
      memcpy(x, node->x, (node->n + 1) * sizeof(double));
    }
    free_node(node);
    free_set(h);
    return z;
  }

  branch(node, z);

  while (size(h) > 0) {
    node_t *node = pop_node_from_set(h);

    succ(node, h, m, n, a, b, c, node->h, 1, floor(node->xh), &z, x);
    succ(node, h, m, n, a, b, c, node->h, -1, -ceil(node->xh), &z, x);
    free_node(node);
  }

  free_set(h);

  if (z == -INFINITY) {
    return NAN;
  } else {
    return z;
  }
}

int main(void) {
  int m;
  int n;
  double **a;
  double *b;
  double *c;
  double *x;
  int i, j;

  scanf("%d %d", &m, &n);

  a = calloc(m + n, sizeof(double *));
  b = calloc(m + n, sizeof(double));
  c = calloc(n + 1, sizeof(double));
  x = calloc(n + m + 1, sizeof(double));

  for (i = 0; i < n; i++) {
    scanf("%lf", &c[i]);
  }
  for (i = 0; i < m; i++) {
    a[i] = calloc(n + 1, sizeof(double));
    for (j = 0; j < n; j++) {
      scanf("%lf", &a[i][j]);
    }
  }
  for (i = 0; i < m; i++) {
    scanf("%lf", &b[i]);
  }

  double z = intopt(m, n, a, b, c, x);
  printf("%lf\n", z);

  free(b);
  for (i = 0; i < m; i++) {
    free(a[i]);
  }
  free(a);
  free(c);
  free(x);
}
