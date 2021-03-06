#Commands
CC=gcc
MPICC=mpicc
RM=rm -f
MAKEDEPEND=makedepend

#flags
CFLAGS=-Wall -g
LDFLAGS=-lm
OPTI_FLAG =-Ofast

#Files
LBM_SOURCES=main.c lbm_phys.c lbm_init.c lbm_struct.c lbm_comm.c lbm_config.c
LBM_HEADERS=$(wildcards:*.h)
LBM_OBJECTS=$(LBM_SOURCES:.c=.o)

TARGET=lbm display

all: $(TARGET)

%.o: %.c
	$(MPICC) $(CFLAGS) $(OPTI_FLAG) -c -o $@ $<

lbm: $(LBM_OBJECTS)
	$(MPICC) $(CFLAGS) $(OPTI_FLAG) -o $@ $^ $(LDFLAGS)

display: display.c
	$(CC) $(CFLAGS) $(OPTI_FLAG) -o $@ display.c

clean:
	$(RM) $(LBM_OBJECTS)
	$(RM) $(TARGET)

depend:
	$(MAKEDEPEND) -Y. $(LBM_SOURCES) display.c

run:
	mpirun -np 4 ./lbm

diff:
	diff -y resultat.raw ../Reference/resultat.raw
	
gif:
	./gen_animate_gif.sh ./resultat.raw ./output.gif

opengif:
	xdg-open ./output.gif

edit:
	subl *.c
	subl *.h

.PHONY: clean all depend

# DO NOT DELETE

main.o: lbm_config.h lbm_struct.h lbm_phys.h lbm_comm.h lbm_init.h
lbm_phys.o: lbm_config.h lbm_struct.h lbm_phys.h lbm_comm.h
lbm_init.o: lbm_phys.h lbm_struct.h lbm_config.h lbm_comm.h lbm_init.h
lbm_struct.o: lbm_struct.h lbm_config.h
lbm_comm.o: lbm_comm.h lbm_struct.h lbm_config.h
lbm_config.o: lbm_config.h
display.o: lbm_struct.h lbm_config.h
