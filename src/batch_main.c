/* Standalone JSON batch runner for the GeneScript C runtime. */
#define _POSIX_C_SOURCE 200809L
#include <ctype.h>
#include <errno.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double gs_calc_gc_content(const char *seq);
long gs_calc_length(const char *seq);
char *gs_calc_reverse(const char *seq);
char *gs_calc_complement(const char *seq);
char *gs_calc_reverse_complement(const char *seq);
char *gs_calc_translate(const char *seq);

#define MAX_SEQUENCES 100
#define MAX_REQUEST_BYTES (1024u * 1024u)
#define MAX_SEQUENCE_BASES 100000
#define MAX_WORKERS 4
#define MAX_OPERATIONS 6

typedef struct { char *data; size_t len; } String;
typedef enum { OP_LENGTH, OP_GC, OP_REVERSE, OP_COMPLEMENT, OP_RC, OP_TRANSLATE } Op;
typedef struct {
    int valid;
    char error[128];
    long length;
    double gc;
    char *strings[MAX_OPERATIONS];
} Result;
typedef struct {
    String *sequences;
    Result *results;
    Op *operations;
    size_t nseq, nops, next;
    pthread_mutex_t mutex;
    int mutex_ready;
} Work;

static char input[MAX_REQUEST_BYTES + 1];
static size_t input_len, pos;
static const char *parse_error;

static void skip_ws(void) {
    while (pos < input_len && (input[pos] == ' ' || input[pos] == '\t' ||
           input[pos] == '\r' || input[pos] == '\n')) pos++;
}

static int hex4(unsigned *out) {
    unsigned v = 0;
    for (int i = 0; i < 4; i++) {
        if (pos >= input_len) return 0;
        unsigned char c = (unsigned char)input[pos++];
        v <<= 4;
        if (c >= '0' && c <= '9') v |= c - '0';
        else if (c >= 'a' && c <= 'f') v |= c - 'a' + 10;
        else if (c >= 'A' && c <= 'F') v |= c - 'A' + 10;
        else return 0;
    }
    *out = v;
    return 1;
}

static int append_utf8(char *dst, size_t cap, size_t *n, unsigned cp) {
    unsigned char bytes[4]; size_t count;
    if (cp <= 0x7f) { bytes[0] = (unsigned char)cp; count = 1; }
    else if (cp <= 0x7ff) {
        bytes[0] = 0xc0 | (cp >> 6); bytes[1] = 0x80 | (cp & 63); count = 2;
    } else if (cp <= 0xffff) {
        bytes[0] = 0xe0 | (cp >> 12); bytes[1] = 0x80 | ((cp >> 6) & 63);
        bytes[2] = 0x80 | (cp & 63); count = 3;
    } else if (cp <= 0x10ffff) {
        bytes[0] = 0xf0 | (cp >> 18); bytes[1] = 0x80 | ((cp >> 12) & 63);
        bytes[2] = 0x80 | ((cp >> 6) & 63); bytes[3] = 0x80 | (cp & 63); count = 4;
    } else return 0;
    if (*n + count >= cap) return 0;
    memcpy(dst + *n, bytes, count); *n += count;
    return 1;
}

/* Decode JSON strings, including escapes and surrogate pairs. */
static int parse_string(String *out) {
    if (pos >= input_len || input[pos++] != '"') return 0;
    size_t cap = input_len - pos + 1, n = 0;
    char *s = (char *)malloc(cap);
    if (!s) { parse_error = "out of memory"; return 0; }
    while (pos < input_len) {
        unsigned char c = (unsigned char)input[pos++];
        if (c == '"') { s[n] = '\0'; out->data = s; out->len = n; return 1; }
        if (c < 0x20) { free(s); return 0; }
        if (c != '\\') {
            if (n + 1 >= cap) { free(s); return 0; }
            s[n++] = (char)c;
            continue;
        }
        if (pos >= input_len) { free(s); return 0; }
        c = (unsigned char)input[pos++];
        switch (c) {
            case '"': case '\\': case '/': if (n + 1 >= cap) { free(s); return 0; } s[n++] = (char)c; break;
            case 'b': if (n + 1 >= cap) { free(s); return 0; } s[n++] = '\b'; break;
            case 'f': if (n + 1 >= cap) { free(s); return 0; } s[n++] = '\f'; break;
            case 'n': if (n + 1 >= cap) { free(s); return 0; } s[n++] = '\n'; break;
            case 'r': if (n + 1 >= cap) { free(s); return 0; } s[n++] = '\r'; break;
            case 't': if (n + 1 >= cap) { free(s); return 0; } s[n++] = '\t'; break;
            case 'u': {
                unsigned cp;
                if (!hex4(&cp)) { free(s); return 0; }
                if (cp >= 0xd800 && cp <= 0xdbff) {
                    unsigned low;
                    if (pos + 2 > input_len || input[pos++] != '\\' || input[pos++] != 'u' || !hex4(&low) || low < 0xdc00 || low > 0xdfff) { free(s); return 0; }
                    cp = 0x10000 + ((cp - 0xd800) << 10) + (low - 0xdc00);
                } else if (cp >= 0xdc00 && cp <= 0xdfff) { free(s); return 0; }
                if (!append_utf8(s, cap, &n, cp)) { free(s); return 0; }
                break;
            }
            default: free(s); return 0;
        }
    }
    free(s); return 0;
}

static int expect(char c) { skip_ws(); if (pos >= input_len || input[pos] != c) return 0; pos++; return 1; }

static int op_from_string(const String *s, Op *op) {
#define MATCH(text, value) (s->len == sizeof(text)-1 && memcmp(s->data, text, sizeof(text)-1) == 0 && ((*op = value), 1))
    return MATCH("LENGTH", OP_LENGTH) || MATCH("GC_CONTENT", OP_GC) ||
           MATCH("REVERSE", OP_REVERSE) || MATCH("COMPLEMENT", OP_COMPLEMENT) ||
           MATCH("REVERSE_COMPLEMENT", OP_RC) || MATCH("TRANSLATE", OP_TRANSLATE);
#undef MATCH
}

static int parse_string_array(String **items, size_t *count, size_t limit) {
    if (!expect('[')) return 0;
    skip_ws();
    if (expect(']')) { *items = NULL; *count = 0; return 1; }
    String *a = (String *)calloc(limit, sizeof(*a));
    if (!a) { parse_error = "out of memory"; return 0; }
    size_t n = 0;
    for (;;) {
        if (n == limit || n >= limit) { parse_error = "array limit exceeded"; goto fail; }
        skip_ws();
        if (!parse_string(&a[n])) { if (!parse_error) parse_error = "malformed JSON string or truncated request"; goto fail; }
        n++;
        skip_ws();
        if (expect(']')) break;
        if (!expect(',')) { parse_error = "expected ',' or ']'"; goto fail; }
    }
    *items = a; *count = n; return 1;
fail:
    for (size_t i = 0; i < n; i++) free(a[i].data);
    free(a); return 0;
}

static int parse_request(Op **ops, size_t *nops, String **seqs, size_t *nseq) {
    int have_ops = 0, have_seqs = 0;
    String *keys = NULL;
    if (!expect('{')) { parse_error = "request must be a JSON object"; return 0; }
    skip_ws();
    while (!expect('}')) {
        String key = {0};
        if (!parse_string(&key) || !expect(':')) { free(key.data); parse_error = "malformed JSON object"; goto fail; }
        if (key.len == 10 && memcmp(key.data, "operations", 10) == 0) {
            free(key.data);
            if (have_ops) { parse_error = "duplicate operations field"; goto fail; }
            String *values = NULL; size_t n = 0;
            if (!parse_string_array(&values, &n, MAX_OPERATIONS)) goto fail;
            Op *parsed = (Op *)calloc(n ? n : 1, sizeof(*parsed));
            if (!parsed) { for (size_t i=0;i<n;i++) free(values[i].data); free(values); parse_error="out of memory"; goto fail; }
            for (size_t i = 0; i < n; i++) {
                if (!op_from_string(&values[i], &parsed[i])) { parse_error = "unsupported operation"; }
                for (size_t j = 0; j < i; j++) if (parsed[j] == parsed[i] && !parse_error) parse_error = "duplicate operation";
                free(values[i].data);
            }
            free(values);
            if (parse_error) { free(parsed); goto fail; }
            *ops = parsed; *nops = n; have_ops = 1;
        } else if (key.len == 9 && memcmp(key.data, "sequences", 9) == 0) {
            free(key.data);
            if (have_seqs) { parse_error = "duplicate sequences field"; goto fail; }
            if (!parse_string_array(&keys, nseq, MAX_SEQUENCES)) goto fail;
            *seqs = keys; have_seqs = 1;
        } else {
            free(key.data); parse_error = "unknown request field"; goto fail;
        }
        skip_ws();
        if (expect('}')) break;
        if (!expect(',')) { parse_error = "expected ',' or '}'"; goto fail; }
    }
    skip_ws();
    if (pos != input_len) { parse_error = "trailing data after JSON request"; goto fail; }
    if (!have_ops || !have_seqs) { parse_error = "request requires operations and sequences"; goto fail; }
    if (*nops == 0 || *nseq == 0) { parse_error = "operations and sequences must not be empty"; goto fail; }
    return 1;
fail:
    if (keys) { for (size_t i=0;i<*nseq;i++) free(keys[i].data); free(keys); }
    if (have_seqs) *seqs = NULL;
    if (have_ops) { free(*ops); *ops = NULL; }
    return 0;
}

static void json_string(const char *s, size_t len) {
    putchar('"');
    for (size_t i = 0; i < len; i++) {
        unsigned char c = (unsigned char)s[i];
        switch (c) {
            case '"': fputs("\\\"", stdout); break;
            case '\\': fputs("\\\\", stdout); break;
            case '\b': fputs("\\b", stdout); break;
            case '\f': fputs("\\f", stdout); break;
            case '\n': fputs("\\n", stdout); break;
            case '\r': fputs("\\r", stdout); break;
            case '\t': fputs("\\t", stdout); break;
            default: if (c < 0x20) printf("\\u%04x", c); else putchar(c);
        }
    }
    putchar('"');
}

static void free_results(Result *r, size_t n, size_t nops) {
    if (!r) return;
    for (size_t i = 0; i < n; i++) for (size_t j = 0; j < nops; j++) free(r[i].strings[j]);
    free(r);
}

static void do_one(Work *w, size_t i) {
    String *s = &w->sequences[i]; Result *r = &w->results[i];
    if (s->len > MAX_SEQUENCE_BASES) { snprintf(r->error, sizeof r->error, "sequence exceeds %d bases", MAX_SEQUENCE_BASES); return; }
    for (size_t k = 0; k < s->len; k++) {
        int c = toupper((unsigned char)s->data[k]);
        if (c != 'A' && c != 'C' && c != 'G' && c != 'T') {
            snprintf(r->error, sizeof r->error, "invalid nucleotide at position %zu", k + 1); return;
        }
    }
    char *seq = (char *)malloc(s->len + 1);
    if (!seq) { snprintf(r->error, sizeof r->error, "out of memory"); return; }
    for (size_t k=0; k<s->len; k++) seq[k] = (char)toupper((unsigned char)s->data[k]);
    seq[s->len] = '\0';
    r->valid = 1;
    for (size_t j = 0; j < w->nops; j++) {
        switch (w->operations[j]) {
            case OP_LENGTH: r->length = gs_calc_length(seq); break;
            case OP_GC: r->gc = gs_calc_gc_content(seq); break;
            case OP_REVERSE: r->strings[j] = gs_calc_reverse(seq); break;
            case OP_COMPLEMENT: r->strings[j] = gs_calc_complement(seq); break;
            case OP_RC: r->strings[j] = gs_calc_reverse_complement(seq); break;
            case OP_TRANSLATE: r->strings[j] = gs_calc_translate(seq); break;
        }
        if ((w->operations[j] == OP_REVERSE || w->operations[j] == OP_COMPLEMENT ||
             w->operations[j] == OP_RC || w->operations[j] == OP_TRANSLATE) && !r->strings[j]) {
            r->valid = 0; snprintf(r->error, sizeof r->error, "out of memory");
            for (size_t k=0;k<w->nops;k++) { free(r->strings[k]); r->strings[k]=NULL; }
            break;
        }
    }
    free(seq);
}

static void *worker(void *arg) {
    Work *w = (Work *)arg;
    for (;;) {
        pthread_mutex_lock(&w->mutex); size_t i = w->next++; pthread_mutex_unlock(&w->mutex);
        if (i >= w->nseq) break;
        do_one(w, i);
    }
    return NULL;
}

static const char *op_name(Op op) {
    static const char *names[] = {"LENGTH", "GC_CONTENT", "REVERSE", "COMPLEMENT", "REVERSE_COMPLEMENT", "TRANSLATE"};
    return names[op];
}

static void emit_error(const char *message) {
    fputs("{\"error\":", stdout); json_string(message, strlen(message)); fputs("}\n", stdout);
    fprintf(stderr, "genescript-batch: %s\n", message);
}

static int valid_utf8(const unsigned char *s, size_t n) {
    for (size_t i = 0; i < n;) {
        unsigned c = s[i++];
        if (c < 0x80) continue;
        unsigned need, cp, min;
        if (c >= 0xc2 && c <= 0xdf) { need=1; cp=c&0x1f; min=0x80; }
        else if (c >= 0xe0 && c <= 0xef) { need=2; cp=c&0x0f; min=0x800; }
        else if (c >= 0xf0 && c <= 0xf4) { need=3; cp=c&0x07; min=0x10000; }
        else return 0;
        if (i + need > n) return 0;
        for (unsigned j=0;j<need;j++) {
            unsigned d=s[i++]; if ((d&0xc0)!=0x80) return 0;
            cp=(cp<<6)|(d&0x3f);
        }
        if (cp < min || cp > 0x10ffff || (cp >= 0xd800 && cp <= 0xdfff)) return 0;
    }
    return 1;
}

int main(void) {
    for (;;) {
        int c = getchar();
        if (c == EOF) break;
        if (input_len == MAX_REQUEST_BYTES) { emit_error("request exceeds 1 MiB limit"); return 2; }
        input[input_len++] = (char)c;
    }
    input[input_len] = '\0';
    if (!valid_utf8((const unsigned char *)input, input_len)) { emit_error("request is not valid UTF-8"); return 2; }
    Op *ops = NULL; String *seqs = NULL; size_t nops = 0, nseq = 0;
    if (!parse_request(&ops, &nops, &seqs, &nseq)) {
        emit_error(parse_error ? parse_error : "malformed or truncated JSON request"); return 2;
    }
    Result *results = (Result *)calloc(nseq, sizeof(*results));
    if (!results) { emit_error("out of memory"); goto fatal; }
    Work w = { .sequences=seqs, .results=results, .operations=ops, .nseq=nseq, .nops=nops };
    pthread_t threads[MAX_WORKERS]; size_t started = 0;
    if (pthread_mutex_init(&w.mutex, NULL) == 0) {
        w.mutex_ready = 1;
        size_t count = nseq > MAX_WORKERS ? MAX_WORKERS - 1 : nseq - 1;
        for (; started < count; started++) {
            int rc = pthread_create(&threads[started], NULL, worker, &w);
            if (rc != 0) { fprintf(stderr, "genescript-batch: worker creation failed (%s); completing with available workers\n", strerror(rc)); break; }
        }
        worker(&w); /* collector participates; also guarantees completion after partial thread creation */
        for (size_t i=0;i<started;i++) {
            int rc = pthread_join(threads[i], NULL);
            if (rc != 0) {
                fprintf(stderr, "genescript-batch: pthread_join failed (%s); terminating without freeing worker state\n", strerror(rc));
                fflush(stderr);
                /* A failed join does not prove the worker has stopped. Exiting
                   the process is the only safe way to avoid freeing memory it
                   may still access. */
                _Exit(1);
            }
        }
        pthread_mutex_destroy(&w.mutex);
    } else {
        fprintf(stderr, "genescript-batch: mutex initialization failed; using serial fallback\n");
        for (size_t i=0;i<nseq;i++) do_one(&w, i);
    }
    fputs("{\"results\":[", stdout);
    for (size_t i=0;i<nseq;i++) {
        if (i) putchar(',');
        fputs("{\"sequence\":", stdout); json_string(seqs[i].data, seqs[i].len);
        if (!results[i].valid) { fputs(",\"error\":", stdout); json_string(results[i].error, strlen(results[i].error)); }
        else {
            for (size_t j=0;j<nops;j++) {
                printf(",\"%s\":", op_name(ops[j]));
                if (ops[j] == OP_LENGTH) printf("%ld", results[i].length);
                else if (ops[j] == OP_GC) printf("%.2f", results[i].gc);
                else json_string(results[i].strings[j], strlen(results[i].strings[j]));
            }
        }
        putchar('}');
    }
    fputs("]}\n", stdout);
    free_results(results, nseq, nops);
    for (size_t i=0;i<nseq;i++) free(seqs[i].data);
    free(seqs); free(ops);
    return ferror(stdout) ? 1 : 0;
fatal:
    if (results) free(results);
    for (size_t i=0;i<nseq;i++) free(seqs[i].data);
    free(seqs); free(ops); return 2;
}
