CC = FCCpx
CFLAGS = -Nclang -Ofast
LDFLAGS = -lstdc++

ifeq ($(fugaku),no)
CC = g++
CFLAGS = -Wall -Wextra
CFLAGS += -O3 -march=native
endif

CFLAGS += -fopenmp
LDFLAGS += -fopenmp

ifeq ($(seq),std)
CFLAGS += -D_SEQ_STD
endif

ifeq ($(seq),pdq)
CFLAGS += -D_SEQ_PDQ
endif

ifeq ($(merge),heap)
CFLAGS += -D_MERGE_WITH_HEAP
endif

ifeq ($(merge),sort)
CFLAGS += -D_MERGE_BY_SORT
endif

OBJS = main.o
PROGRAM = a.out

.PHONY: all clean

all: $(PROGRAM)

clean:
	rm -f $(PROGRAM) *.o

$(PROGRAM): $(OBJS)
	$(CC) $(OBJS) $(LDFLAGS) -o $(PROGRAM)

%.o: %.cpp
	$(CC) $(CFLAGS) -c $<
