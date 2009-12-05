
#
# Makefile para compilar programas en C++
#

PROYECTO  = bignum
FUENTES   = *.cc

C_FLAGS   = -Wall -O2 -g $(INCDIR)
LIBRERIAS = -lm -lstdc++
COMPILAR  = g++ $(C_FLAGS) -c $< -o $@
ENLAZAR   = g++ -o $(PROYECTO) $(OBJS) $(LIBRERIAS)

SRCS = $(wildcard $(FUENTES))
OBJS = $(SRCS:.cc=.o)

all: $(PROYECTO)

$(PROYECTO): .depend $(OBJS)
	@echo "Enlazando $(OBJS) -> $@"
	@$(ENLAZAR)

%.o: %.cc .depend
	@echo "Compilando $< -> $@"
	@$(COMPILAR)

run: $(PROYECTO)
clean:
	@rm -f $(OBJS)
	@rm -f .depend

.PHONY: all run clean

.DEFAULT:
	@g++ -M *.cc > .depend

sinclude .depend
