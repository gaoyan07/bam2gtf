CC      =	gcc
CFLAGS  =	-Wall -O2 -Wno-unused-variable -Wno-unused-but-set-variable -Wno-unused-function
DFLAGS  =	-g -Wall  
HTSLIB_DIR = ./htslib
HTSLIB  =   $(HTSLIB_DIR)/libhts.a
LIB     =	$(HTSLIB) -lm -lz -lpthread
COMP_LIB=	-lz
INCLUDE = -I $(HTSLIB_DIR)


BIN_DIR =	.
SRC_DIR =   .

HTS_ALL =   hts_all
SOURCE  =	$(wildcard ${SRC_DIR}/*.c) 
COMP_SOURCE = compare_gtf.c utils.c gtf.c
OBJS    =	$(SOURCE:.c=.o)

BIN     =	$(BIN_DIR)/gtools

RMATS       =   $(BIN_DIR)/rmats_gtools
RGDB_DEBUG   =   $(BIN_DIR)/rgdb_gtools
GDB_DEBUG   =   $(BIN_DIR)/gdb_gtools
NOR_DEBUG   =   $(BIN_DIR)/debug_gtools
DMARCRO 	=	-D __DEBUG__
RDMARCRO    =   -D __DEBUG__ -D _RMATS_
RMARCRO     =   -D _RMATS_
COMP_GTF	= 	$(BIN_DIR)/comp-gtf
GDB_COMP	=   $(BIN_DIR)/gdb_comp-gtf
COMP_D		=	-D COMP_MAIN

.c.o:
		$(CC) -c $(CFLAGS) $(INCLUDE) $< -o $@

#all:       $(SOURCE) $(BIN) 
all:       $(HTS_ALL) $(BIN) 
gdb_gtools: $(SOURCE) $(GDB_DEBUG) 
rmats_gtools: $(SOURCE) $(RMATS)
rgdb_gtools: $(SOURCE) $(RGDB_DEBUG) 
debug_gtools: $(SOURCE) $(NOR_DEBUG)
comp-gtf:	$(COMP_GTF)
gdb_comp-gtf: $(GDB_COMP)


$(HTS_ALL):
	cd $(HTSLIB_DIR); make;
$(BIN): $(OBJS)
		$(CC) $(OBJS) -o $@ $(LIB)

$(RMATS):
		$(CC) $(CFLAGS) $(SOURCE) $(RMARCRO) -o $@ $(LIB)
$(GDB_DEBUG):
		$(CC) $(DFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LIB)
$(RGDB_DEBUG):
		$(CC) $(DFLAGS) $(SOURCE) $(RDMARCRO) -o $@ $(LIB)
$(NOR_DEBUG):
		$(CC) $(CFLAGS) $(SOURCE) $(DMARCRO) -o $@ $(LIB)
$(COMP_GTF):
		$(CC) $(CFLAGS) $(COMP_SOURCE) $(COMP_D) -o $@ $(COMP_LIB)
$(GDB_COMP):
		$(CC) $(DFLAGS) $(COMP_SOURCE) $(COMP_D) -o $@ $(COMP_LIB)


clean:
		rm -f $(SRC_DIR)/*.o $(BIN) $(RMATS)

clean_debug:
		rm -f $(SRC_DIR)/*.o $(GDB_DEBUG) $(RGDB_DEBUG) $(NOR_DEBUG) $(GDB_COMP)
