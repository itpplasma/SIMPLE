FFLAGS = -O
OBJS  =  nleq1.o wnorm.o linalg_nleq1.o zibconst.o zibsec.o zibmon.o

nleq1: main_nleq1.o $(OBJS)
	$(FC) $(DFLAGS) main_nleq1.o $(OBJS) -o nleq1

nleq1-easy: main_nleq1_easy.o $(OBJS) nleq1e.o
	$(FC) $(DFLAGS) main_nleq1_easy.o $(OBJS) nleq1e.o -o nleq1-easy

clean:
	rm -f $(OBJS) nleq1 nleq1-easy main_nleq1.o main_nleq1_easy.o nleq1e.o
