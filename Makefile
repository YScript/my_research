CC = ifort

BIN = a.out 

OBJ = main.o globalPara.o random.o init_coor.o caculate_e.o change_part.o attchange_snake.o replica.o

MODSRCS = globalPara.mod random.mod init_coor.mod caculate_e.mod change_part.mod attchange_snake.mod replica.mod

.PHONY: all clean

$(BIN) : $(OBJ)
	$(CC) -o $(BIN) $(OBJ)

globalPara.mod : globalPara.o globalPara.f90
	$(CC) -c globalPara.f90 -o globalPara.mod

globalPara.o : globalPara.f90
	$(CC) -c globalPara.f90 

random.mod : globalPara.mod random.o random.f90
	$(CC) -c random.f90 

random.o : globalPara.mod random.f90
	$(CC) -c random.f90 

caculate_e.mod : globalPara.mod caculate_e.o caculate_e.f90
	$(CC) -c caculate_e.f90

caculate_e.o : globalPara.mod caculate_e.f90
	$(CC) -c caculate_e.f90

init_coor.mod : globalPara.mod random.mod caculate_e.mod init_coor.o init_coor.f90
	$(CC) -c init_coor.f90 

init_coor.o : globalPara.mod random.mod caculate_e.mod init_coor.f90
	$(CC) -c init_coor.f90

change_part.mod : globalPara.mod change_part.o change_part.f90
	$(CC) -c change_part.f90

change_part.o : globalPara.mod change_part.f90
	$(CC) -c change_part.f90

attchange_snake.mod : globalPara.mod random.mod caculate_e.mod change_part.mod attchange_snake.o attchange_snake.f90
	$(CC) -c attchange_snake.f90

attchange_snake.o : globalPara.mod random.mod caculate_e.mod change_part.mod attchange_snake.f90
	$(CC) -c attchange_snake.f90

replica.mod : globalPara.mod random.mod caculate_e.mod change_part.mod attchange_snake.mod replica.o replica.f90
	$(CC) -c replica.f90

replica.o : globalPara.mod random.mod caculate_e.mod change_part.mod attchange_snake.mod replica.f90
	$(CC) -c replica.f90

main.o : globalPara.mod random.mod caculate_e.mod init_coor.mod replica.mod main.f90
	$(CC) -c main.f90

clean :
	rm $(BIN) $(OBJ) $(MODSRCS)
