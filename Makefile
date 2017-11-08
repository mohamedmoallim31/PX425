########################################
# Makefile for PX425 2017 Assignment 2 #
########################################

# C compiler and flags
CC=gcc -O3
CFLAGS=-Wall -Wextra -I/opt/local/include

# Command to use for linking and executable
LD=gcc
EXE=wave_eqn
LDFLAGS=-L/opt/local/lib -lpng -lm

all: wave_eqn

# Default build target
wave_eqn : makePNG.o wave_eqn.o mt19937ar.o
	$(CC) $(CFLAGS) -o $(EXE) makePNG.o wave_eqn.o  mt19937ar.o $(LDFLAGS)

# Solution
#sol : makePNG.o wave_eqn_sol.o mt19937ar.o
#	$(CC) $(CFLAGS) -o $(EXE)_sol makePNG.o wave_eqn_sol.o  mt19937ar.o $(LDFLAGS)

# Extreme
#xtreme : makePNG.o wave_eqn_xtreme.o mt19937ar.o
#	$(CC) $(CFLAGS) -o $(EXE)_xtreme makePNG.o wave_eqn_xtreme.o  mt19937ar.o $(LDFLAGS)


movie:
	yes | ffmpeg -f image2 -i snapshot%08d.png -pix_fmt yuv420p -vcodec libx264  -level 30 -maxrate 10000000 -bufsize 10000000 -b 1200k -f mp4 movie.mp4

# Purge build files and executable
clean :
	rm -rf *.o *.mod ./$(EXE)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

