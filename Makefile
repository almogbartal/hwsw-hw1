CC      := cc
CFLAGS  := -O0 -Wall -Wextra
LDLIBS  := -lm

BINS := firefly_sync_unoptimized firefly_sync_optimized

all: $(BINS)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(BINS)

.PHONY: all clean
