CC      := cc
CFLAGS  := -g -O0 -Wall -Wextra -fno-omit-frame-pointer
LDLIBS  := -lm

BINS := firefly_sync_unoptimized firefly_sync_optimized

all: $(BINS)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(BINS)

.PHONY: all clean
