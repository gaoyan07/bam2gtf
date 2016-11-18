CC      =	gcc
CFLAGS  =	-Wall -O3 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
DFLAGS  =	-g -Wall  
HTSLIB_DIR = ./htslib
HTSLIB  =   $(HTSLIB_DIR)/libhts.a
LIB     =	$(HTSLIB) -lm -lz -lpthread
INCLUDE = -I ./htslib


BIN_DIR =	.
SRC_DIR =   .

HTS_ALL =   hts_all
SOURCE  =	$(wildcard ${SRC_DIR}/*.c) 
OBJS    =	$(SOURCE:.c=.o)

BIN     =	$(BIN_DIR)/gtools

GDB_DEBUG   =   $(BIN_DIR)/gdb_gtools
NOR_DEBUG   =   $(BIN_DIR)/debug_gtools
DMARCRO =	-D __DEBUG__

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

#all:       $(SOURCE) $(BIN) 
all:       $(HTS_ALL) $(BIN) 
gdb_gtools: $(SOURCE) $(GDB_DEBUG) 
debug_gtools: $(SOURCE) $(NOR_DEBUG)


_end_flag = '$'
$(HTS_ALL):
	cd $(HTSLIB_DIR); make; cd ../
$(BIN): $(OBJS)
		$(CC) $(OBJS) -o $@ $(LIB)

$(GDB_DEBUG):
		$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LIB)
$(NOR_DEBUG):
		$(CC) $(CFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LIB)

clean:
		rm -f $(SRC_DIR)/*.o $(BIN)

clean_debug:
		rm -f $(SRC_DIR)/*.o $(GDB_DEBUG) $(NOR_DEBUG)
