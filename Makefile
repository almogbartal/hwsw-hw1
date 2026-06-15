CC      := cc
CFLAGS  := -g -O0 -Wall -Wextra -fno-omit-frame-pointer
LDLIBS  := -lm

BINS := unoptimized opt_1 opt_1_2 opt_3 opt_4 opt_5 opt_6 opt_7 opt_1_2_3 unoptimized_orig

all: $(BINS)

%: %.c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(BINS)

.PHONY: all clean
