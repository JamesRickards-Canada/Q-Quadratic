#Call "./configure" to configure the library, then "make" to build the project. "make clean" clears up the .o files.

#Name of the output library
TARGET = qquadratic

#The object files to build
OBJS = qq_base.o qq_bqf.o qq_bqf_int.o qq_geometry.o qq_quat.o qq_quat_int.o qq_visual.o

#Nothing after here should be modified, unless you know what you are doing.

#PARI_LIB is folder where libpari.so/libpari.dylib is found, PARI_INCLUDE is where the pari.h header file is found, and PARI_CFG is the location of pari.cfg.
PARI_LOC = $(TARGET).cfg
PARI_CFG = $(shell grep "CFG=" "qquadratic.cfg" -s | cut -d"'" -f2)
ifeq ($(PARI_CFG), )
	PARI_CFG = /usr/local/lib/pari/pari.cfg
endif
PARI_LIB = $(shell grep "libdir=" "$(PARI_CFG)" -s | cut -d"'" -f2)
PARI_INCLUDE = $(shell grep "includedir=" "$(PARI_CFG)" -s | cut -d"'" -f2)/pari

#Naming the library file to include the version of pari/gp.
VER = $(shell grep "pari_release=" "$(PARI_CFG)" -s | cut -d"'" -f2 | tr . - | cut -d"-" -f1,2,3)
DYN = lib$(TARGET)-$(VER).so

#Compiling options
CC = cc
CFLAGS = -O3 -Wall -fno-strict-aliasing -fPIC -march=native
RM = rm -f

#System check as -shared option for linker fails on MacOS
ifneq ($(shell uname -s), Darwin)
    OS_FLAG = -Wl,-shared
else
    OS_FLAG =
endif

#Recipes
all: $(DYN)

#Build the shared library object
$(DYN): $(OBJS)
	$(CC) -o $@ -shared	$(CFLAGS) $(OS_FLAG) $(OBJS) -lc -lm -L$(PARI_LIB) -lpari

#Make the object files
%.o: %.c
	$(CC) -c $(CFLAGS) -I. -I$(PARI_INCLUDE) $<

#Clear all .o files
clean:
	$(RM) *.o $(ALL)
