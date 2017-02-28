CC = g++
CCFLAGS = -O2
interpolation:
	$(CC) $(CCFLAGS) -ointerpolation.exe utils/rational.o utils/vect.o interpolation.cpp
gauss:
	$(CC) $(CCFLAGS) -ogauss.exe utils/rational.o utils/matrix.o utils/vect.o gauss.cpp
clean:
	$(RM) *.o run.exe



