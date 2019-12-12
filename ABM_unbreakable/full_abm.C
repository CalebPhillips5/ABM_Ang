// Running the main function as a function to test for vegf_prod values.
#include "general_libraries.h"
#include "main_abm_oxy.h"
#include "cell_abm.h"
#include "mpi.h" 

int main(int argc, char** argv){
  
  MPI_Init(&argc, &argv);
  
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);
  
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  int name_len;
  MPI_Get_processor_name(processor_name, &name_len);
  
  MPI_Comm proc_comm;
  MPI_Comm_split(MPI_COMM_WORLD, world_rank, world_rank, &proc_comm);
  
  int proc_rank, proc_size;
  MPI_Comm_rank(proc_comm, &proc_rank);
  MPI_Comm_size(proc_comm, &proc_size);

  // You need to pass proc_comm to each libmesh code in the following way:
  //void main_code(MPI_Comm lib_comm,....
  //And then use it on libmeshinit:
  //  LibMeshInit init (argc, argv, lib_comm);

  // initialize parameters
  double nutrient_diffusion,nutrient_production,nutrient_consumption,vegf_diffusion,vegf_production,vegf_consumption;
  ifstream myfile;
  std::ifstream read;
  std::string line;
  // open and read file
  read.open("files-All_8000.dat"); // FILE NAME GOES HERE
  if(!read.is_open())
    std::cout << "Error opening.\n";
  // assemble matrix from file
  std::vector< std::vector<double> > PARAMS;
  while(std::getline(read,line)){
    std::istringstream(line) >> nutrient_diffusion >> nutrient_production >> nutrient_consumption >> vegf_diffusion >> vegf_production >> vegf_consumption;
    std::vector<double> read_par(6,0.0);
    read_par[0] = 1000.00*(0.5+nutrient_diffusion);
    read_par[1] = 200.00*(0.5+nutrient_production);
    read_par[2] = 0.8000*(0.5+nutrient_consumption);
    read_par[3] = 8000.0*(0.5+vegf_diffusion);
    read_par[4] = 300.00*(0.5+vegf_production);
    read_par[5] = 10.0000*(0.5+vegf_consumption); //2.0
    PARAMS.push_back(read_par);
  }
  read.close();
  const unsigned int size = PARAMS.size();
  cout << "Number of set of parameters = " << size << endl;
  int start = 0;
  for(int i=0;i<world_rank;i++)
    start+=(size+i)/world_size;
  int size_local = (size+world_rank)/world_size;
  int end = start+size_local-1;
  printf("WORLD RANK/SIZE: %d/%d \t ROW RANK/SIZE: %d/%d \t MATRIX SIZE: %d\t START/END: %d/%d\n",	world_rank, world_size, proc_rank, proc_size,(size+world_rank)/world_size,start,end);
  /*
  for(int r=start;r<=end;r++)
    for(int c=0;c<6;c++)
      printf("WORLD RANK %d \t VAR/VALUE: %d/%lf\n",	world_rank, c,PARAMS[r][c]);
  */
  //  HOW YOUR CODE WILL LOOK LIKE!
  
  for(int j=start;j<=end;j++){
    main_code(argc,argv,proc_comm,PARAMS[j],j);
    cout << "Finished = " << j << endl;
  }
  
  MPI_Finalize();
}
//	   main_code(argc,argv,200,5,0.2,1000,200,4,1);
