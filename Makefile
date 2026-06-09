CC      := cc
CFLAGS  := -O0 -Wall -Wextra
LDLIBS  := -lm

BIN := firefly_kuramoto

all: $(BIN)

$(BIN): $(BIN).c
	$(CC) $(CFLAGS) -o $@ $< $(LDLIBS)

clean:
	rm -f $(BIN)

.PHONY: all clean
