CC      := cc
CFLAGS  := -g -O0 -Wall -Wextra -fno-omit-frame-pointer
LDLIBS  := -lm

BINS := optimized unoptimized unoptimized_orig 

all: $(BINS)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(BINS)

.PHONY: all clean
