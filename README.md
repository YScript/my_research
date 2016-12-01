# my_research
the project is about codes in my Master's research ;
the project is done by myself with fortran language, f90 files ,90 standard.
it is a Process Oriented Program;

contains:
  main.f90            --the main program;
  globalPara.f90      --global parameter module;
  random.f90          --a random generator function module;
  caculate_e.f90      --caculate system configuration energy module;
  init_coor.f90       --system configuration coordinates and energy and atom character initialization module;
  change_part.f90     --atom's coordinates and character changing module; 
  attchange_snake.f90 --choose atoms to attempt change its properties，find which system-config is allowed to be accepted,by snake moving;
  replica.f90         --replication exchanging module,contains choose atom,attempt change and accept possibility;
  
##  this project is mainly about self-assembly of grafted Homopolymer PS and diblock Copolymer PS-b-PDMAEMA solution mixed system;
    the single solvent is methanol with atom type 3,firstly.It would be double solvents THF and methanol later ;
