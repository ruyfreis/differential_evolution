#Compilador, opcoes de compilacao.
CC        = g++
#CFLAGS    = -g -Wall 	 #use para debug
CFLAGS   = -O3    #use para compilar com otimização
#CFLAGS   = -g -Wall  #use para fazer valgrind
LDFLAGS   = -lm 

#Fontes 
SRCDIR = ./
SRCS   = $(shell ls $(SRCDIR)/*.cpp 2>/dev/null)
DEPS   = $(shell ls $(SRCDIR)/*.h 2>/dev/null)
OBJS   = $(SRCS:.c=.o)

#Nome do Executavel
exec1  = lib_de

#Regras de compilacao
all : $(exec1)

$(exec1): $(OBJS) $(DEPS)
	$(CC) $(CFLAGS) -o $(exec1) $(OBJS) $(LDFLAGS)

clean : 
	rm -f $(exec1) $(SRCDIR)/*.o
	
