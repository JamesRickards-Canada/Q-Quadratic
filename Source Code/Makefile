PARI_LIB = /usr/local/lib
PARI_INCLUDE = /usr/local/include

OBJS = qq_bqf.o qq_geometry.o qq_quat_int.o qq_base.o qq_bqf_int.o qq_quat.o qq_visual.o

TARGET = qquadratic
SHELL  = /bin/sh

CFLAGS     = -O3 -Wall -fno-strict-aliasing
EXTRACFLAGS=

CC         = cc
CPPFLAGS   = -I. -I$(PARI_INCLUDE)
LD         = cc
LDFLAGS    = -O3 -Wall -fno-strict-aliasing    -Wl,--export-dynamic 
MODLD      = cc
MODLDFLAGS = -shared  $(CFLAGS) $(DLCFLAGS) -Wl,-shared 
EXTRAMODLDFLAGS = -lc -lm -L$(PARI_LIB) -lpari
EXTRALIBS  =

DLCFLAGS   = -fPIC

DYN = lib$(TARGET).so
ALL = $(DYN)

all: $(DYN)

$(DYN): $(OBJS)
	$(MODLD) -o $@ $(MODLDFLAGS) $(EXTRACFLAGS) $(OBJS) $(EXTRAMODLDFLAGS)
	cp $@ ..

%.o: %.c
	$(CC) -c $(CFLAGS) $(EXTRACFLAGS) $(CPPFLAGS) $(DLCFLAGS) $<
clean:
	-$(RM) *.o $(ALL)
