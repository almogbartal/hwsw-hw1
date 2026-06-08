/*
 * Naive delta-based Huffman compressor for 1024x1024 grayscale raw images.
 *
 * Usage:
 *   ./huffimg compress   <input_dir> <output_dir>   (reads *.raw, writes *.huf)
 *   ./huffimg decompress <input_dir> <output_dir>   (reads *.huf, writes *.raw)
 *
 * Each .raw file must be exactly 1024*1024 = 1048576 bytes (one byte per pixel).
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <dirent.h>
#include <sys/stat.h>

#define IMG_WIDTH  1024
#define IMG_HEIGHT 1024
#define IMG_SIZE   (IMG_WIDTH * IMG_HEIGHT)
#define ALPHABET   256
#define MAX_CODE_LEN 256

typedef struct Node {
    uint64_t freq;
    int symbol;             /* 0..255 for leaves, -1 for internal nodes */
    struct Node *left;
    struct Node *right;
} Node;

typedef struct {
    char bits[MAX_CODE_LEN + 1];  /* '0'/'1' string, null terminated */
    int  length;
} Code;

/* ---------- bit I/O ---------- */

typedef struct {
    FILE *fp;
    uint8_t byte;
    int nbits;        /* number of bits currently buffered (0..7) */
} BitWriter;

static void bw_write_bit(BitWriter *bw, int bit) {
    bw->byte = (uint8_t)((bw->byte << 1) | (bit & 1));
    bw->nbits++;
    if (bw->nbits == 8) {
        fputc(bw->byte, bw->fp);
        bw->byte = 0;
        bw->nbits = 0;
    }
}

static void bw_flush(BitWriter *bw) {
    if (bw->nbits > 0) {
        bw->byte = (uint8_t)(bw->byte << (8 - bw->nbits));
        fputc(bw->byte, bw->fp);
        bw->byte = 0;
        bw->nbits = 0;
    }
}

typedef struct {
    FILE *fp;
    uint8_t byte;
    int nbits;        /* bits remaining in the current byte */
} BitReader;

static int br_read_bit(BitReader *br) {
    if (br->nbits == 0) {
        int c = fgetc(br->fp);
        if (c == EOF) return -1;
        br->byte = (uint8_t)c;
        br->nbits = 8;
    }
    int bit = (br->byte >> 7) & 1;
    br->byte = (uint8_t)(br->byte << 1);
    br->nbits--;
    return bit;
}

/* ---------- Huffman tree ---------- */

static Node *new_leaf(int sym, uint64_t freq) {
    Node *n = (Node *)malloc(sizeof(Node));
    n->freq = freq;
    n->symbol = sym;
    n->left = n->right = NULL;
    return n;
}

static Node *new_internal(Node *l, Node *r) {
    Node *n = (Node *)malloc(sizeof(Node));
    n->freq = (l ? l->freq : 0) + (r ? r->freq : 0);
    n->symbol = -1;
    n->left = l;
    n->right = r;
    return n;
}

static void free_tree(Node *n) {
    if (!n) return;
    free_tree(n->left);
    free_tree(n->right);
    free(n);
}

/* Naive O(n^2) Huffman: repeatedly scan for the two smallest-freq nodes. */
static Node *build_tree(const uint64_t freqs[ALPHABET]) {
    Node *pool[ALPHABET];
    int n = 0;
    for (int i = 0; i < ALPHABET; i++) {
        if (freqs[i] > 0) pool[n++] = new_leaf(i, freqs[i]);
    }
    if (n == 0) return NULL;
    if (n == 1) {
        /* Only one symbol used: wrap in an internal node so its code is "0". */
        return new_internal(pool[0], NULL);
    }

    while (n > 1) {
        int i1 = 0, i2 = 1;
        if (pool[i2]->freq < pool[i1]->freq) { int t = i1; i1 = i2; i2 = t; }
        for (int i = 2; i < n; i++) {
            if (pool[i]->freq < pool[i1]->freq) {
                i2 = i1;
                i1 = i;
            } else if (pool[i]->freq < pool[i2]->freq) {
                i2 = i;
            }
        }
        Node *parent = new_internal(pool[i1], pool[i2]);

        int big = i1 > i2 ? i1 : i2;
        int small = i1 < i2 ? i1 : i2;
        pool[small] = parent;
        pool[big]   = pool[n - 1];
        n--;
    }
    return pool[0];
}

static void gen_codes(Node *node, Code codes[ALPHABET], char *buf, int depth) {
    if (!node) return;
    if (node->symbol >= 0) {
        if (depth == 0) {
            /* Shouldn't happen because build_tree wraps the single-symbol case. */
            codes[node->symbol].bits[0] = '0';
            codes[node->symbol].bits[1] = '\0';
            codes[node->symbol].length = 1;
        } else {
            buf[depth] = '\0';
            memcpy(codes[node->symbol].bits, buf, (size_t)depth + 1);
            codes[node->symbol].length = depth;
        }
        return;
    }
    if (node->left) {
        buf[depth] = '0';
        gen_codes(node->left, codes, buf, depth + 1);
    }
    if (node->right) {
        buf[depth] = '1';
        gen_codes(node->right, codes, buf, depth + 1);
    }
}

/* ---------- compress / decompress ---------- */

/*
 * File format (little-endian, machine native — fine for homework use):
 *   "HUF1"                              4 bytes magic
 *   uint32_t width
 *   uint32_t height
 *   uint64_t total_bits                 number of valid bits in the payload
 *   uint64_t freqs[256]                 frequency table
 *   bitstream payload                   ceil(total_bits / 8) bytes
 */

static int compress_file(const char *in_path, const char *out_path) {
    FILE *in = fopen(in_path, "rb");
    if (!in) { perror(in_path); return -1; }

    uint8_t *pixels = (uint8_t *)malloc(IMG_SIZE);
    size_t r = fread(pixels, 1, IMG_SIZE, in);
    fclose(in);
    if (r != IMG_SIZE) {
        fprintf(stderr, "%s: expected %d bytes, got %zu\n", in_path, IMG_SIZE, r);
        free(pixels);
        return -1;
    }

    uint8_t *deltas = (uint8_t *)malloc(IMG_SIZE);
    deltas[0] = pixels[0];
    for (int i = 1; i < IMG_SIZE; i++) {
        deltas[i] = (uint8_t)(pixels[i] - pixels[i - 1]);
    }
    free(pixels);

    uint64_t freqs[ALPHABET] = {0};
    for (int i = 0; i < IMG_SIZE; i++) freqs[deltas[i]]++;

    Node *root = build_tree(freqs);
    Code codes[ALPHABET];
    memset(codes, 0, sizeof(codes));
    char buf[MAX_CODE_LEN + 1];
    gen_codes(root, codes, buf, 0);

    uint64_t total_bits = 0;
    for (int i = 0; i < ALPHABET; i++) {
        total_bits += freqs[i] * (uint64_t)codes[i].length;
    }

    FILE *out = fopen(out_path, "wb");
    if (!out) { perror(out_path); free(deltas); free_tree(root); return -1; }

    fwrite("HUF1", 1, 4, out);
    uint32_t w = IMG_WIDTH, h = IMG_HEIGHT;
    fwrite(&w, sizeof(uint32_t), 1, out);
    fwrite(&h, sizeof(uint32_t), 1, out);
    fwrite(&total_bits, sizeof(uint64_t), 1, out);
    fwrite(freqs, sizeof(uint64_t), ALPHABET, out);

    BitWriter bw = { out, 0, 0 };
    for (int i = 0; i < IMG_SIZE; i++) {
        Code *c = &codes[deltas[i]];
        for (int j = 0; j < c->length; j++) {
            bw_write_bit(&bw, c->bits[j] - '0');
        }
    }
    bw_flush(&bw);

    fclose(out);
    free(deltas);
    free_tree(root);
    return 0;
}

static int decompress_file(const char *in_path, const char *out_path) {
    FILE *in = fopen(in_path, "rb");
    if (!in) { perror(in_path); return -1; }

    char magic[4];
    if (fread(magic, 1, 4, in) != 4 || memcmp(magic, "HUF1", 4) != 0) {
        fprintf(stderr, "%s: bad magic\n", in_path);
        fclose(in);
        return -1;
    }
    uint32_t w, h;
    uint64_t total_bits;
    uint64_t freqs[ALPHABET];
    fread(&w, sizeof(uint32_t), 1, in);
    fread(&h, sizeof(uint32_t), 1, in);
    fread(&total_bits, sizeof(uint64_t), 1, in);
    fread(freqs, sizeof(uint64_t), ALPHABET, in);

    Node *root = build_tree(freqs);
    if (!root) {
        fprintf(stderr, "%s: empty frequency table\n", in_path);
        fclose(in);
        return -1;
    }

    uint64_t total_pixels = (uint64_t)w * (uint64_t)h;
    uint8_t *deltas = (uint8_t *)malloc(total_pixels);

    BitReader br = { in, 0, 0 };
    for (uint64_t i = 0; i < total_pixels; i++) {
        Node *node = root;
        while (node->symbol < 0) {
            int bit = br_read_bit(&br);
            if (bit < 0) {
                fprintf(stderr, "%s: unexpected EOF in bitstream\n", in_path);
                free(deltas); free_tree(root); fclose(in);
                return -1;
            }
            node = bit ? node->right : node->left;
            if (!node) {
                fprintf(stderr, "%s: invalid bit sequence\n", in_path);
                free(deltas); free_tree(root); fclose(in);
                return -1;
            }
        }
        deltas[i] = (uint8_t)node->symbol;
    }

    fclose(in);
    free_tree(root);

    uint8_t *pixels = (uint8_t *)malloc(total_pixels);
    pixels[0] = deltas[0];
    for (uint64_t i = 1; i < total_pixels; i++) {
        pixels[i] = (uint8_t)(pixels[i - 1] + deltas[i]);
    }
    free(deltas);

    FILE *out = fopen(out_path, "wb");
    if (!out) { perror(out_path); free(pixels); return -1; }
    fwrite(pixels, 1, total_pixels, out);
    fclose(out);
    free(pixels);
    return 0;
}

/* ---------- directory walker ---------- */

static int has_suffix(const char *s, const char *suf) {
    size_t ls = strlen(s), lf = strlen(suf);
    return ls >= lf && strcmp(s + ls - lf, suf) == 0;
}

static int process_dir(const char *in_dir, const char *out_dir, int compress_mode) {
    DIR *dir = opendir(in_dir);
    if (!dir) { perror(in_dir); return -1; }

    if (mkdir(out_dir, 0755) != 0) {
        /* ignore EEXIST */
    }

    const char *in_ext  = compress_mode ? ".raw" : ".huf";
    const char *out_ext = compress_mode ? ".huf" : ".raw";

    struct dirent *entry;
    int ok = 0, fail = 0;
    while ((entry = readdir(dir)) != NULL) {
        if (entry->d_name[0] == '.') continue;
        if (!has_suffix(entry->d_name, in_ext)) continue;

        char in_path[2048];
        char out_path[2048];
        char base[1024];

        size_t name_len = strlen(entry->d_name);
        size_t ext_len  = strlen(in_ext);
        size_t base_len = name_len - ext_len;
        if (base_len >= sizeof(base)) continue;
        memcpy(base, entry->d_name, base_len);
        base[base_len] = '\0';

        snprintf(in_path,  sizeof(in_path),  "%s/%s",   in_dir,  entry->d_name);
        snprintf(out_path, sizeof(out_path), "%s/%s%s", out_dir, base, out_ext);

        printf("%s %s -> %s\n",
               compress_mode ? "compress  " : "decompress",
               in_path, out_path);

        int rc = compress_mode
            ? compress_file(in_path, out_path)
            : decompress_file(in_path, out_path);
        if (rc == 0) ok++; else fail++;
    }
    closedir(dir);
    printf("Done: %d ok, %d failed\n", ok, fail);
    return fail == 0 ? 0 : 1;
}

int main(int argc, char **argv) {
    if (argc != 4) {
        fprintf(stderr,
                "Usage:\n"
                "  %s compress   <input_dir> <output_dir>\n"
                "  %s decompress <input_dir> <output_dir>\n",
                argv[0], argv[0]);
        return 1;
    }
    int compress_mode;
    if      (strcmp(argv[1], "compress")   == 0) compress_mode = 1;
    else if (strcmp(argv[1], "decompress") == 0) compress_mode = 0;
    else {
        fprintf(stderr, "First arg must be 'compress' or 'decompress'\n");
        return 1;
    }
    return process_dir(argv[2], argv[3], compress_mode);
}
