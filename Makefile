# Makefile
# Automatic EQ

CC      = gcc
CFLAGS  = -g -std=c99 -Wall -lncurses -lportaudio -lsndfile -lfftw3 -lm
EXE  = AutoEQ
SRCS = surgEQnew.c inputlib.c

#OBJS = $(SRCS:.c=.o) /usr/local/lib/libportaudio.dylib /usr/local/lib/libsndfile.dylib 

$(EXE):	$(OBJS) $(HDRS) $(SRCS)
	$(CC) $(CFLAGS) $(SRCS) -o $@ 

clean:
	rm -f *~ core $(EXE) *.o