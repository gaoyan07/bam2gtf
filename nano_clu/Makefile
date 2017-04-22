CC      = gcc
CFLAGS  = -c -Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
DFLAGS  = -g -Wall
LDFLAGS = -lz -lm -lpthread

BIN_DIR =
SRC_DIR = 
SOURCE  = $(wildcard $(SRC_DIR)*.c)
OBJS    = $(SOURCE:.c=.o)

PROG    = $(BIN_DIR)rest
DEBUG   = $(BIN_DIR)gdb_rest
LIT     = $(BIN_DIR)lit_rest
DMARCRO = -D __DEBUG__
LMARCRO = -D __LIT__

all: $(SOURCE) $(PROG)

$(PROG): $(OBJS)
	$(CC) $(OBJS) -o $@ $(LDFLAGS)

.c.o:
	$(CC) $(CFLAGS) $< -o $@ 

$(DEBUG):
	$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LDFLAGS)
$(LIT):
	$(CC) $(DFLAGS) $(SOURCE) $(LMARCRO) $(DMARCRO) -o $@ $(LDFLAGS)

clean:
	rm -f *.o $(PROG) $(DEBUG)

