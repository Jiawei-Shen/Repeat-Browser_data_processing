include ../../../inc/common.mk

L += -lm -lz
MYLIBDIR = ../../../lib/${MACHTYPE}
MYLIBS =  ${MYLIBDIR}/jkhgap.a ${MYLIBDIR}/jkweb.a

MYF = TEstat

O = $(MYF).o

all: ${O} $(MYLIBS)
	${CC} ${COPT} -o $(MYF) $O ${MYLIBS} $L
