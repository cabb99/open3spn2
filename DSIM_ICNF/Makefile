CMPLR      =  gcc

CFLGS      =  -O3 -Wall -finline-functions -funroll-all-loops 

INCLD      =  -I/usr/local/include

LIBSD      =  -L/usr/local/lib -lm

EXECU      =  icnf.exe

SRCES      =  genr_boxd.c \
              genr_dnac.c \
              genr_rtmp_adna.c \
              genr_rtmp_bdna.c \
              genr_topo.c \
              genr_trns.c \
              prog_icnf.c \
              genr_ions.c \
              read_ions.c \
              read_sequ.c \
              wrte_lammps.c \
              wrte_rpsf.c \
              wrte_rxyz.c

OBJCS      =  ${SRCES:.c=.o}

all: ${EXECU}

%.o: %.c
	${CMPLR} ${CFLGS} ${DFLGS} ${INCLD} -c $< -o $@

${EXECU}: ${OBJCS}
	${CMPLR} ${CFLGS} ${DFLGS} ${INCLD} -o $@ ${OBJCS} ${LIBSD}

clean-objs:
	rm -f ${OBJCS}

clean:
	rm -f ${OBJCS} ${EXECU}
