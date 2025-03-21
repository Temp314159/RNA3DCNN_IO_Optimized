CC=g++ -O3

LAPACK=lapack_LINUX.a blas_LINUX.a libf2c.a

v1=ExtractList.o pdb.o pdbWrite.o VECTOR.o calc_crmsd_VECTOR.o calc_crmsd_RNA.o reshapeResidue.o calc_score.o $(LAPACK)
v2=Voxelizate.o pdb.o pdbWrite.o VECTOR.o calc_crmsd_VECTOR.o calc_crmsd_RNA.o reshapeResidue.o calc_score.o $(LAPACK)

all: o1
o1:$(v1)
	$(CC) -no-pie $(v1) -o ExtractList
o2:$(v2)
	$(CC) -no-pie $(v2) -o Voxelizate
	
clean:
	rm -f *.o *.out *.h.gch
c:clean

%.o:%.cpp
	$(CC) -c -fPIC $< -o $@
