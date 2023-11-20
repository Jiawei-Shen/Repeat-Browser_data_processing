include ../../../inc/common.mk

L += $(MYSQLLIBS) -lm -lz
MYLIBDIR = ../../../lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkhgap.a ${MYLIBDIR}/jkweb.a

MYF = NearByGenes

O = $(MYF).o

all: ${O} $(MYLIBS)
	${CC} ${COPT} -o $(MYF) $O ${MYLIBS} $L
