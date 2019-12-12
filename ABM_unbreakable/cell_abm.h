#ifndef CELL_ABM_H_
#define CELL_ABM_H_

class Cell{
 public:
  //========== Cell position ==========
  double x,y;
  //========== Cell properties ==========
  double N_radius,C_radius,A_radius,uptake,cal;
  double v[2],F[2],prev_F[2];
  int time,state,prev_state,stalks,flag1,flag2;
  int closer,anasto,activate;
  void set(double X,double Y,double RN,double R,double RA,double Uptake, int Time,int State);
  void print();
};

void update_states(list<Cell> &Cells_local,list<Cell*> &Cells_tip,list<Cell> &Cells_vessel_start,list<Cell> &Cell_vessel_end,int time,Ran& ran,EquationSystems& equation_systems,const int outside_cells);
void potential_vegf(double& psi_x,double& psi_y,double r_x,double r_y,double vegf_x,double vegf_y,double R_N,double R,double M,int m);
void compute_forces(list<Cell> &Cells_local,list<Cell*> &Cells_tip,list<Cell> &Cells_vessel_start,list<Cell> &Cells_vessel_end,EquationSystems& equation_systems,double height,int& outside_cells,int time);
inline double distance(Cell cell_A,double pos_x,double pos_y){return sqrt(pow(cell_A.x-pos_x,2)+pow(cell_A.y-pos_y,2));}
inline double point_distance(const Point& A,const Point& B){return sqrt(pow(A(0)-B(0),2)+pow(A(1)-B(1),2));}
void potential_rep(double& psi_x,double& psi_y,double r_x,double r_y,double R_N,double R,double M,int m);
void divide_tumor(list<Cell> &Cells_local,std::list<Cell>::iterator &it,Ran& ran,double height);
void print_order(list<Cell> &Cells_local,double domain_diameter,string s,int file_number,int t);
void save_cells(list<Cell> &Cells_local,double domain_diameter,string s,int file_number,int t);
void divide_endo(list<Cell> &Cells_local,list<Cell*> &Cells_tip,std::list<Cell>::iterator &it);
void init_cond_cells(list<Cell> &Cells_local,const Parameters& all_parameters,Ran& ran);
void potential_adh(double& phi_x,double& phi_y,double r_x,double r_y,double R_A,int n);
void restart_function(list<Cell> &Cells_local,list<Cell*> &Cells_tip,string name_s);
void normal(double& N_x,double& N_y,double x_cel,double y_cel,double height);
void restart_save(list<Cell> &Cells_local,string name_s);
double euclid_norm(double x,double y);

#endif
