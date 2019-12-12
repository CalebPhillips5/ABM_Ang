#include "general_libraries.h"
#include "cell_abm.h"
#include "main_abm_oxy.h"
#include "chrono"

//###########################################################################
//                    0  - Dead Cells
//                    1  - Tumor Cells
//                    2  - Proliferative Tumor Cells
//                    3  - Hypoxic Tumor Cells
//                    4  - Dying Tumor Cells
//                    5  - G1 Tumor Cells
//                    6  - Normoxic Cells
//                    7  - Endothelial Cells
//                    8  - Tip Cells
//                    9  - Stalk Cells
// 		      10 - Activated Tip Cells
// 		      11 - Growing Stalk Cells
//###########################################################################

vector<Cell> Cells_global;

Number exact_value_ox(const Point& ,const Parameters& ,const std::string& ,const std::string& ){return 0.6;}
double triangle_area(const Point& e1,const Point& e2,const Point& e3){return 0.5*fabs((e1(0)-e3(0))*(e2(1)-e1(1))-(e1(0)-e2(0))*(e3(1)-e1(1)));}
double scalar_prod(const Point& e1,const Point& e2){return e1(0)*e2(0)+e1(1)*e2(1);}

void main_code(int argc, char** argv,MPI_Comm lib_comm,vector<double> Parameters,int file_number){
  // The parameters follow this order
  //double nut_d,double nut_coeff,double con_t,double vegf_diff,double vegf_prod,double vegf_cons
  //========== Initialize the library ==========
  LibMeshInit init (argc, argv,lib_comm);
  //========== Declare cells list ==========****
  list <Cell> Cells_local;
  list <Cell> Cells_vessel_start;
  list <Cell> Cells_vessel_end;
  list <Cell*> Tip_local;
  vector<Cell> Cells_tip;
  //========== Read argumments ==========
  GetPot input_file("options.in");
  const bool read_solution         = input_file("read_solution",0);
  const bool verbose               = input_file("verbose",0);
  const bool unbreakable           = input_file("unbreakable",0);
  const bool deactivation          = input_file("deactivation",0);
  const bool anastomosis           = input_file("anastomosis",0);
  const int max_tip_cells          = input_file("max_tip_cells",100);
  const bool dirichlet             = input_file("dirichlet",0);
  const unsigned int n_elements    = input_file("n_elements",20);
  const unsigned int n_timesteps   = input_file("n_timesteps",100);
  const unsigned int initial_tum   = input_file("initial_tum",20);
  int init_timestep                = input_file("init_timestep",0);
  const unsigned int print_inter   = input_file("print_inter",100);
  const unsigned int ic_type       = input_file("ic_type",1);
  const unsigned int rand_seed     = input_file("rand_seed",7);
  const unsigned int max_outside   = input_file("max_outside",1500);
  const std::string mesh_name      = input_file("mesh_name", "test.msh");
  const unsigned int read_mesh     = input_file("read_mesh",0);
  const double initial_con         = input_file("initial_con",20.0);
  const double con_b               = input_file("con_b",1.0);
  const double con_n               = input_file("con_n",1.0);
  const double end_r               = input_file("end_r",0.1);
  const double vegf_ths            = input_file("vegf_ths",0.1);
  const double domain_diameter     = input_file("domain_diameter",1.0);
  const double time_step           = input_file("time_step",0.05);
  const double prol_intens         = input_file("prol_intens",1.0);
  const double nucleus_radius      = input_file("nucleus_radius",5.295);
  const double cell_radius         = input_file("cell_radius",9.953);
  const double action_prop         = input_file("action_prop",1.214);
  const double c_ccr               = input_file("c_ccr",10.0);
  const double c_eea               = input_file("c_eea",0.588836);
  //int tip_cell_c                   = input_file("tip_cell_c",0);
  Ran ran(rand_seed);
  int outside_cells = 0;
  int total_tumor = initial_tum;
  char buffer[30];
  ofstream out_data;
  std::stringstream cell_holder;
  std::string extension = ".txt";
  cell_holder << "cells_dataSA" << file_number << extension;
  std::string file = cell_holder.str();
  out_data.open(file);
  //========== Create a uniform mesh ==========
  SerialMesh mesh(init.comm());
  if(!read_solution){
    if(read_mesh){
      mesh.read(mesh_name);
    }
    else{
      MeshTools::Generation::build_cube (mesh,
					 n_elements,
					 n_elements,
					 0,
					 0., domain_diameter,
					 0., domain_diameter,
					 0., 0.,
					 TRI3);
      mesh.write("mesh.msh");
    }
  }
  else{
    mesh.read(mesh_name);
  }
  if(verbose){mesh.print_info();}
  double H_MAX = 0;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  for ( ; el != end_el ; ++el){
    const Elem* elem = *el;
    if(H_MAX<elem->hmax())
      H_MAX = elem->hmax();
  }
  //========== Create an equation systems object and set a simulation-specific parameter ==========
  //0 - nut_d, 1 - nut_coeff, 2 - con_t, 3 - vegf_diff, 4 - vegf_prod, 5 - vegf_cons
  
  EquationSystems equation_systems(mesh);
  equation_systems.parameters.set<Real> ("c_ccr")           = c_ccr;
  equation_systems.parameters.set<Real> ("c_eea")           = c_eea;
  equation_systems.parameters.set<Real> ("vegf_prod")       = Parameters[4];
  equation_systems.parameters.set<Real> ("vegf_cons")       = Parameters[5];
  equation_systems.parameters.set<Real> ("nut_coeff")       = Parameters[1];
//  equation_systems.parameters.set<Real> ("nut_coeff")       = 0;
  equation_systems.parameters.set<Real> ("vegf_diff")       = Parameters[3];
  equation_systems.parameters.set<Real> ("vegf_ths")        = vegf_ths;
  equation_systems.parameters.set<Real> ("tip_distance")    = 8;
  equation_systems.parameters.set<Real> ("repulsion_distance")    = 1.1;
  equation_systems.parameters.set<Real> ("nut_d")           = Parameters[0];
  equation_systems.parameters.set<Real> ("con_b")           = con_b;
  equation_systems.parameters.set<Real> ("con_t")           = Parameters[2];
  equation_systems.parameters.set<Real> ("con_n")           = con_n;
  equation_systems.parameters.set<Real> ("end_r")           = end_r;
  equation_systems.parameters.set<Real> ("time_step_size")  = time_step;
  equation_systems.parameters.set<int>  ("initial_tumor")   = total_tumor;
  equation_systems.parameters.set<bool> ("deactivation")    = deactivation;
  equation_systems.parameters.set<bool> ("unbreakable")     = unbreakable;
  equation_systems.parameters.set<bool> ("anastomosis")     = anastomosis;
  equation_systems.parameters.set<int>  ("max_tip_cells")   = max_tip_cells;
  equation_systems.parameters.set<int>  ("ic_type")         = ic_type;
  equation_systems.parameters.set<Real> ("initial_con")     = initial_con;
  equation_systems.parameters.set<Real> ("nucleus_radius")  = nucleus_radius;
  equation_systems.parameters.set<Real> ("cell_radius")     = cell_radius;
  equation_systems.parameters.set<Real> ("action_radius")   = action_prop*cell_radius;
  equation_systems.parameters.set<Real> ("apop_time")       = 8.6;
  equation_systems.parameters.set<Real> ("g1_time")         = 9.0;
  equation_systems.parameters.set<Real> ("lysing_time")     = 6.0;
  equation_systems.parameters.set<Real> ("calc_time")       = 360.0;
  equation_systems.parameters.set<Real> ("cellc_time")      = 18.0;
  equation_systems.parameters.set<Real> ("hypoxic_thrs")    = 0.3;
  equation_systems.parameters.set<Real> ("prol_intes")      = prol_intens/1.921166667;
  equation_systems.parameters.set<Real> ("apop_intes")      = 1.0/786.61;
  equation_systems.parameters.set<Real> ("delta_tt")        = 1.0;
  equation_systems.parameters.set<Real> ("f_NS")            = 1.0;
  equation_systems.parameters.set<Real> ("lambda_cell")     = 0.1;
  equation_systems.parameters.set<Real> ("domain_diameter") = domain_diameter;
  equation_systems.parameters.set<Real> ("max_out_cells")   = max_outside;
  equation_systems.parameters.set<Real> ("h_max_mesh")      = H_MAX;
  equation_systems.parameters.set<list <Cell> *>("l_cells") = &Cells_local;
  //========== Declare the system ==========
  TransientLinearImplicitSystem &vegf = equation_systems.add_system<TransientLinearImplicitSystem> ("VEGF");
  TransientLinearImplicitSystem &nutr = equation_systems.add_system<TransientLinearImplicitSystem> ("Nutrient");
  if(!read_solution){
    //========== Timestep always 0 ==============
    init_timestep = 0;
    //========== Declare the variables ==========
    unsigned int v_var = vegf.add_variable("vegf_var", FIRST);
    unsigned int nut_var = nutr.add_variable("nut_var", FIRST);
    //nutr.add_variable("nut_var", FIRST);
    //========== Matrix assembly and initial condition functions ==========
    vegf.attach_assemble_function(assemble_vegf);
    nutr.attach_assemble_function(assemble_ox);
    nutr.attach_init_function(initial_condition_ox);
    //========== Initialize the data structures for the equation system ==========
    std::set<boundary_id_type> boundary_ids_one;
    if(read_mesh){
      boundary_ids_one.insert(0);
    }
    else{
      boundary_ids_one.insert(0);
      boundary_ids_one.insert(1);
      boundary_ids_one.insert(2);
      boundary_ids_one.insert(3);
    }
    std::vector<unsigned int> variables;
    variables.push_back(v_var);
    ConstFunction<Number> bc_value(0.0);
    if(dirichlet){
      DirichletBoundary dirichlet_bc_one(boundary_ids_one,variables,&bc_value);
      vegf.get_dof_map().add_dirichlet_boundary(dirichlet_bc_one);
    }
    std::vector<unsigned int> variablesn;
    variablesn.push_back(nut_var);
    ConstFunction<Number> bcnut_value(1.0);
    //DirichletBoundary dirichlet_bc_two(boundary_ids_one,variablesn,&bcnut_value);
    //nutr.get_dof_map().add_dirichlet_boundary(dirichlet_bc_two);
    equation_systems.init();
    init_cond_cells(Cells_local,equation_systems.parameters,ran);
  }
  else{
    std::stringstream ss;
    std::string ext = ".e";
    ss << "saved_solution" << file_number << "_" << setfill('0') << setw(5) << init_timestep << ext;
    std::string results = ss.str();
    //equation_systems.read(results, libMeshEnums::READ);
    equation_systems.read(results);
    vegf.update();
    nutr.update();
    //========== Initialize the data structures for the equation system ==========
    unsigned int v_var = vegf.variable_number("vegf_var");
    unsigned int nut_var = nutr.variable_number("nut_var");
    std::set<boundary_id_type> boundary_ids_one;
    if(read_mesh){
      boundary_ids_one.insert(0);
    }
    else{
      boundary_ids_one.insert(0);
      boundary_ids_one.insert(1);
      boundary_ids_one.insert(2);
      boundary_ids_one.insert(3);
    }
    std::vector<unsigned int> variables;
    variables.push_back(v_var);
    ConstFunction<Number> bc_value(0.0);
    if(dirichlet){
      DirichletBoundary dirichlet_bc_one(boundary_ids_one,variables,&bc_value);
      vegf.get_dof_map().add_dirichlet_boundary(dirichlet_bc_one);
    }
    std::vector<unsigned int> variablesn;
    variablesn.push_back(nut_var);
    ConstFunction<Number> bcnut_value(1.0);
    //DirichletBoundary dirichlet_bc_two(boundary_ids_one,variablesn,&bcnut_value);
    //nutr.get_dof_map().add_dirichlet_boundary(dirichlet_bc_two);
    vegf.attach_assemble_function(assemble_vegf);
    nutr.attach_assemble_function(assemble_ox);
    equation_systems.reinit();
    std::stringstream abm;
    std::string extension = ".txt";
    abm << "results" << file_number << "_" << setfill('0') << setw(5) << init_timestep << extension;
    std::string Caleb = abm.str();
    restart_function(Cells_local,Tip_local,Caleb);
  }
  if(verbose){
    equation_systems.print_info();
    sprintf(buffer, "exo%d.e",file_number);
    ExodusII_IO(mesh).write_equation_systems(buffer,equation_systems);
    save_cells(Cells_local,domain_diameter,"saida",file_number,0);
    print_order(Cells_local,domain_diameter,"output",file_number,init_timestep);
    save_grad_vegf(equation_systems,"grad_vegf",file_number,0);
    out_data << 0.0 << " " << Cells_global.size() << " " << total_tumor << " " << outside_cells << " " << initial_con << endl;
  }
  //========== Loop over time ==========
  unsigned int t_step = init_timestep;
  do{
    //========== Increase time_step counter ==========
    t_step++;
    vegf.time = t_step*time_step;
    nutr.time = t_step*time_step;
    //========== Copy solution of previous time step ==========
    *vegf.old_local_solution = *vegf.current_local_solution;
    *nutr.old_local_solution = *nutr.current_local_solution;


    //========== Solve ABM system ==========
    update_states(Cells_local,Tip_local,Cells_vessel_start,Cells_vessel_end,t_step,ran,equation_systems,outside_cells);

    if(verbose && t_step%print_inter==0){
      print_order(Cells_local,domain_diameter,"output",file_number,-t_step);
    }
    compute_forces(Cells_local,Tip_local,Cells_vessel_start,Cells_vessel_end,equation_systems,domain_diameter,outside_cells,t_step);
    /////----- Solve systems -----/////
    vegf.solve();
    nutr.solve();
    if(verbose && t_step%print_inter==0){
      save_cells(Cells_local,domain_diameter,"saida",file_number,t_step);
      print_order(Cells_local,domain_diameter,"output",file_number,t_step);
      save_grad_vegf(equation_systems,"grad_vegf",file_number,t_step);
      //-- Save system to restart later --//
      std::stringstream abm;
      std::string extension = ".txt";
      abm << "results" << file_number << "_" << setfill('0') << setw(5) << t_step << extension;
      std::string Caleb = abm.str();
      restart_save(Cells_local,Caleb);
      std::stringstream ss;
      std::string ext = ".e";
      ss << "saved_solution" << file_number << "_" << setfill('0') << setw(5) << t_step << ext;
      std::string results = ss.str();
      libmesh_assert_equal_to(libMesh::processor_id(), 0);
      //equation_systems.write(results, libMeshEnums::WRITE);
      //========== Write output to paraview ==========
      ExodusII_IO exo(mesh);
      exo.append(true);
      exo.write_timestep(buffer, equation_systems, t_step+1, vegf.time);
    }
    double confluence = 0.;
    int Quiescent = 0, Proliferative = 0, Hypoxic = 0, Necrotic = 0, Stalk = 0, Endothelial = 0;
    std::list<Cell>::iterator it;
    for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
      confluence += std::pow((*it).C_radius,2)/std::pow(0.5*domain_diameter,2);
      if((*it).state == 1) Quiescent += 1;
      if((*it).state == 2) Proliferative += 1;
      if((*it).state == 3) Hypoxic += 1;
      if((*it).state == 0) Necrotic += 1;
      if((*it).state == 9 || (*it).state == 11) Stalk += 1;
      if((*it).state == 7) Endothelial += 1;
    }
    int total_t = Quiescent + Proliferative + Hypoxic + Necrotic;
    if(verbose){
      cout << "==================================================" << endl;
      cout << "Confluence          = " << confluence << endl;
      cout << "Number of cells     = " << Cells_local.size() << endl;
      cout << "Tumor cells         = " << total_t << endl;
      cout << "Necrotic cells      = " << Necrotic << endl; // 5
      cout << "Quiescent cells     = " << Quiescent << endl; // 2
      cout << "Proliferative cells = " << Proliferative << endl; // 3 
      cout << "Hypoxic cells       = " << Hypoxic << endl; // 4
      cout << "Endothelial cells   = " << Endothelial << endl; // 7
      cout << "Tip cells           = " << Tip_local.size() << endl; // 1
      cout << "Stalk cells         = " << Stalk << endl; // 6
      cout << "Outside cells       = " << outside_cells << endl;
      cout << "Time                = " << t_step << endl;
      cout << "==================================================" << endl;
    }
    out_data << t_step << " " << Necrotic << " " << Quiescent << " " << Proliferative << " " << Hypoxic << " " << Endothelial << " " << Tip_local.size() << " " << Stalk << endl;
  }while(t_step<n_timesteps);
  out_data.close();
  Cells_global.clear();
}

void save_grad_vegf(EquationSystems& es,string s,int file_number,int t){
  const char *c = s.c_str();
  char n[100],name[200];
  sprintf(n,"%d_%05d.m",file_number,t);
  strcpy(name,c);
  strcat(name,n);
  stringstream ss;
  string name_s;
  ss << name;
  ss >> name_s;
  ofstream out_file;
  out_file.open (name_s);
  //=======*** Saving the data =====***==//
  const MeshBase& mesh = es.get_mesh();
  out_file << "grad = zeros(" << mesh.n_elem() << "," << 4 << ");" <<endl;
  out_file << "grad = [";
  const unsigned int dim = mesh.mesh_dimension();
  TransientLinearImplicitSystem & system = es.get_system<TransientLinearImplicitSystem>("VEGF");
  const unsigned int v_nut = system.variable_number("vegf_var");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_nut;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  for ( ; el != end_el ; ++el){
    const Elem* elem = *el;
    dof_map.dof_indices(elem,dof_indices);
    dof_map.dof_indices(elem,dof_indices_nut,v_nut);
    fe->reinit(elem);
    unsigned int qp=0;
    Gradient gradient_vegf;
    for(unsigned int l=0; l<dof_indices.size(); l++){
      gradient_vegf.add_scaled(dphi[l][qp], system.old_solution(dof_indices[l]));
    }
    out_file << elem->centroid()(0) << " " << elem->centroid()(1) << " " << gradient_vegf(0) << " " << gradient_vegf(1) << endl;
  }
  out_file << "];";
  out_file.close();
}


void assemble_vegf(EquationSystems& es,const std::string& libmesh_dbg_var(system_name)){
  libmesh_assert_equal_to (system_name, "VEGF");
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  TransientLinearImplicitSystem & system = es.get_system<TransientLinearImplicitSystem>("VEGF");
  const unsigned int v_nut = system.variable_number("vegf_var");
  const Real time_size = es.parameters.get<Real>("time_step_size");
  const Real vegf_diff = es.parameters.get<Real>("vegf_diff");
  const Real vegf_prod = es.parameters.get<Real>("vegf_prod");
  const Real vegf_cons = es.parameters.get<Real>("vegf_cons");
  const Real ac_radius = es.parameters.get<Real>("action_radius");
  const Real h_max_msh = es.parameters.get<Real>("h_max_mesh");
  const Real height    = es.parameters.get<Real>("domain_diameter");
  list <Cell> *Cells_local = es.parameters.get<list <Cell> *>("l_cells");
  const double h_bin   = ac_radius+h_max_msh;
  const int number_bins = ceil(height/h_bin);
  const int total_bins  = number_bins*number_bins;
  //========== Generate bins ==========
  vector< list < Cell * > > Cell_Bins(total_bins);
  std::list<Cell>::iterator it;
  for(it = Cells_local->begin(); it != Cells_local->end(); ++it){
    int ix = floor((*it).x/h_bin);
    int jy = floor((*it).y/h_bin);
    int xy = ix+jy*number_bins;
    if((*it).x<0 || (*it).y<0 || (*it).x>height || (*it).y> height || jy>=number_bins || ix>=number_bins || xy >=total_bins){
      cout << "Error" << endl;
      cout << "Cell = ( " << (*it).x << " , " << (*it).y << " ) = ( " << ix << " , " << jy << " ) = " << xy << endl;
      cout << "State = " << (*it).state << endl;
      cout << "total_bins  = " << total_bins << endl;
      cout << "number_bins = " << number_bins << endl;
      cout << "h_bin       = " << h_bin << endl;
      getchar();
    }
    Cell_Bins[xy].push_back(&*it);
  }
  //========== Continue stnd ==========
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_nut;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  for ( ; el != end_el ; ++el){
    const Elem* elem = *el;
    //========== Computed volume fraction ==========
    Point e_center = elem->centroid();
    int ix = floor(e_center(0)/h_bin);
    int jy = floor(e_center(1)/h_bin);
    double p_t,p_n,p_h,p_e,p_s;
    double elem_volume = elem->volume();
    p_t = p_n = p_h = p_e = p_s = 0.;
    for(int xx = -1; xx<=1; xx++){
      if(ix>0 && ix+xx<number_bins){
	for(int yy = -1; yy<=1; yy++){
	  if(jy>0 && jy+yy<number_bins){
	    int bin_xy = (ix+xx)+(jy+yy)*number_bins;
	    std::list<Cell*>::iterator cell_ab;
	    for(cell_ab = Cell_Bins[bin_xy].begin(); cell_ab != Cell_Bins[bin_xy].end(); ++cell_ab){
	      if( (*(*cell_ab)).state == 0 || (*(*cell_ab)).state == 4) continue;
	      else{
		Point p( (*(*cell_ab)).x,(*(*cell_ab)).y,0.);
		double area_computed = 0.;
		if(elem->type()==3){
		  area_computed = area_inside_tri(p,elem->point(0),elem->point(1),elem->point(2),(*(*cell_ab)).C_radius,elem_volume);
		}
		else{
		  area_computed = area_inside_quad(p,elem->point(0),elem->point(1),elem->point(2),elem->point(3),(*(*cell_ab)).C_radius);
		}
		if(area_computed > 0){
		  if((*(*cell_ab)).state == 6) p_n+=(area_computed/elem_volume);
		  else if((*(*cell_ab)).state == 3) p_h+=(area_computed/elem_volume);
		  else if((*(*cell_ab)).prev_state == 7) p_e+=(area_computed/elem_volume);
		  else if((*(*cell_ab)).prev_state == 13) p_e+=(area_computed/elem_volume);
		  else if((*(*cell_ab)).prev_state == 14) p_e+=(area_computed/elem_volume);
		  else p_t+=(area_computed/elem_volume);
		}
	      }
	    }
	  }
	}
      }
    }
    //========== Continue stnd ==========    
    dof_map.dof_indices(elem,dof_indices);
    dof_map.dof_indices(elem,dof_indices_nut,v_nut);
    fe->reinit(elem);
    Ke.resize(dof_indices.size(),dof_indices.size());
    Fe.resize(dof_indices.size());
    for (unsigned int qp=0; qp<qrule.n_points(); qp++){
      Number nut_old = 0.0;
      for(unsigned int l=0; l<phi.size(); l++){
	nut_old += phi[l][qp]*system.old_solution(dof_indices_nut[l]);
      }
      for (unsigned int i=0; i<phi.size(); i++){
	Fe(i) += JxW[qp]*(nut_old
			  +time_size*p_h*vegf_prod
			  )*phi[i][qp];
	for (unsigned int j=0; j<phi.size(); j++){
	  Ke(i,j) += JxW[qp]*time_size*vegf_diff*dphi[i][qp]*dphi[j][qp];
	  Ke(i,j) += JxW[qp]*(1.0+time_size*(p_h*vegf_prod+(p_e+p_s)*vegf_cons))*phi[i][qp]*phi[j][qp];
	}
      }
    }
    dof_map.heterogenously_constrain_element_matrix_and_vector(Ke,Fe,dof_indices);
    system.matrix->add_matrix (Ke,dof_indices);
    system.rhs->add_vector    (Fe,dof_indices);
  }
}

double area_inside_quad(const Point& p,const Point& e1,const Point& e2,const Point& e3,const Point& e4,double radius){
  double elem_volumeA = triangle_area(e1,e2,e3);
  double elem_volumeB = triangle_area(e3,e4,e1);
  double area_div = area_inside_tri(p,e1,e2,e3,radius,elem_volumeA)+area_inside_tri(p,e3,e4,e1,radius,elem_volumeB);
  return area_div;
}

double area_inside_tri(const Point& p,const Point& e1,const Point& e2,const Point& e3,double radius,double elem_volume){
  //========== Check if the center of the cell is inside the element ==========
  int cs_inside = 0;
  double test_area = triangle_area(p,e1,e2)+triangle_area(p,e1,e3)+triangle_area(p,e2,e3);
  if(fabs(elem_volume-test_area)<=1.e-10) cs_inside = 1;
  //========== Define some parameters ========================================*
  double chord[3];
  double angle;
  double area_cs[3];
  double area_t[3];
  std::vector<unsigned int> vert_inside;
  std::vector<unsigned int> vert_out;
  std::vector<unsigned int> edge_inside;
  std::vector<unsigned int> edge_out;
  double area_circ = M_PI*pow(radius,2);
  Point v_e[3];
  v_e[0].add(e1);v_e[1].add(e2);v_e[2].add(e3);
  double t_e[3][2];  
  int part_edge[3][2];
  for(int i=0;i<3;i++){
    for(int j=0;j<2;j++){
      part_edge[i][j]=0;
      t_e[i][j]=0.;
    }
  }
  Point d_v[4];
  Point d_a_e[3];
  //========== Compute the vertices inside the cell =========================**
  for(unsigned int i=0;i<3;i++){
    if(point_distance(p,v_e[i])<=radius){
      vert_inside.push_back(i);
    }
    else{
      vert_out.push_back(i);
    }
  }
  const unsigned int nv = vert_inside.size();
  //========== Compute the edges cutted by the cell =========================**
  for(unsigned int i=0;i<3;i++){
    double a_eq,b_eq,c_eq;
    int push_this = 0;
    int vi = (i+1)%3;
    d_a_e[i](0)=v_e[i](0)-v_e[vi](0);
    d_a_e[i](1)=v_e[i](1)-v_e[vi](1);
    a_eq = scalar_prod(d_a_e[i],d_a_e[i]);
    b_eq = 2*(scalar_prod(d_a_e[i],v_e[i])-scalar_prod(d_a_e[i],p));
    c_eq = scalar_prod(v_e[i],v_e[i])+scalar_prod(p,p)-2*scalar_prod(v_e[i],p)-pow(radius,2);
    if(pow(b_eq,2)>=4.*a_eq*c_eq){
      t_e[i][0]=(-b_eq-sqrt(pow(b_eq,2)-4.*a_eq*c_eq))/(2.*a_eq);
      t_e[i][1]=(-b_eq+sqrt(pow(b_eq,2)-4.*a_eq*c_eq))/(2.*a_eq);
      for(unsigned int j=0;j<2;j++){
	d_v[j].zero();
	d_v[j].add(v_e[i]);
	d_v[j].add_scaled(d_a_e[i],t_e[i][j]);
	if(point_distance(d_v[j],v_e[i]) < point_distance(v_e[i],v_e[vi]) && point_distance(d_v[j],v_e[vi]) < point_distance(v_e[i],v_e[vi])){part_edge[i][j]=1;push_this=1;}
      }
      if(push_this){edge_inside.push_back(i);}
    }
  }
  const unsigned int ne = edge_inside.size();
  //###########################################################################
  //                    0 vertices and 0 edges
  //###########################################################################
  if(nv==0 && ne==0){
    if(!cs_inside){
      return 0.;
    }
    else{
      return area_circ;
    }
  }
  //###########################################################################
  //                    0 vertices and 1 edges
  //###########################################################################
  if(nv==0 && ne==1){
    for(unsigned i=0; i<2; i++){
      d_v[i].zero();
      d_v[i].add(v_e[edge_inside[0]]);
      d_v[i].add_scaled(d_a_e[edge_inside[0]],t_e[edge_inside[0]][i]);
    }
    chord[0] = point_distance(d_v[0],d_v[1]);
    angle = 2.*asin(chord[0]/(2.*radius));
    area_cs[0] = 0.5*pow(radius,2)*(angle-sin(angle));
    if(!cs_inside){
      return area_cs[0];
    }
    else{
      return area_circ-area_cs[0];
    }
  }
  //###########################################################################
  //                    0 vertices and 2 edges
  //###########################################################################
  if(nv==0 && ne==2){
    if(!cs_inside){
      for(unsigned j=0; j<ne; j++){
	for(unsigned i=0; i<2; i++){
	  d_v[i].zero();
	  d_v[i].add(v_e[edge_inside[j]]);
	  d_v[i].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][i]);
	}
	chord[j] = point_distance(d_v[0],d_v[1]);
	angle = 2.*asin(chord[j]/(2.*radius));
	area_cs[j] = 0.5*pow(radius,2)*(angle-sin(angle));
      }
      if(chord[0]>chord[1]){return area_cs[0]-area_cs[1];}
      else{return area_cs[1]-area_cs[0];}
    }
    else{
      for(unsigned j=0; j<ne; j++){
	for(unsigned i=0; i<2; i++){
	  d_v[i].zero();
	  d_v[i].add(v_e[edge_inside[j]]);
	  d_v[i].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][i]);
	}
	chord[j] = point_distance(d_v[0],d_v[1]);
	angle = 2.*asin(chord[j]/(2.*radius));
	area_cs[j] = 0.5*pow(radius,2)*(angle-sin(angle));
      }
      return area_circ-area_cs[0]-area_cs[1];
    }
  }
  //###########################################################################
  //                    0 vertices and 3 edges
  //###########################################################################
  if(nv==0 && ne==3){
    if(!cs_inside){
      for(unsigned j=0; j<ne; j++){
	for(unsigned i=0; i<2; i++){
	  d_v[i].zero();
	  d_v[i].add(v_e[edge_inside[j]]);
	  d_v[i].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][i]);
	}
	chord[j] = point_distance(d_v[0],d_v[1]);
	angle = 2.*asin(chord[j]/(2.*radius));
	area_cs[j] = 0.5*pow(radius,2)*(angle-sin(angle));
      }
      if(chord[0]>chord[1] && chord[0]>chord[2]){return area_cs[0]-area_cs[1]-area_cs[2];}
      else if(chord[1]>chord[0] && chord[1]>chord[2]){return area_cs[1]-area_cs[0]-area_cs[2];}
      else{return area_cs[2]-area_cs[0]-area_cs[1];}
    }
    else{
      for(unsigned j=0; j<ne; j++){
	for(unsigned i=0; i<2; i++){
	  d_v[i].zero();
	  d_v[i].add(v_e[edge_inside[j]]);
	  d_v[i].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][i]);
	}
	chord[j] = point_distance(d_v[0],d_v[1]);
	angle = 2.*asin(chord[j]/(2.*radius));
	area_cs[j] = 0.5*pow(radius,2)*(angle-sin(angle));
      }
      return area_circ-area_cs[0]-area_cs[1]-area_cs[2];
    }
  }
  //###########################################################################
  //                    1 vertices and 2 edges
  //###########################################################################
  if(nv==1 && ne==2){
    for(unsigned j=0; j<ne; j++){
      d_v[j].zero();
      d_v[j].add(v_e[edge_inside[j]]);
      if(part_edge[edge_inside[j]][0]){d_v[j].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][0]);}
      else{                            d_v[j].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][1]);}
    }
    chord[0] = point_distance(d_v[0],d_v[1]);
    angle = 2.*asin(chord[0]/(2.*radius));
    area_cs[0] = 0.5*pow(radius,2)*(angle-sin(angle));
    area_t[0] = triangle_area(v_e[vert_inside[0]],d_v[0],d_v[1]);
    if(!cs_inside){return area_t[0]+area_cs[0];}
    else{
      double m1=(p(1)-v_e[vert_inside[0]](1))/(p(0)-v_e[vert_inside[0]](0));
      double m2=(d_v[0](1)-d_v[1](1))/(d_v[0](0)-d_v[1](0));
      Point Mid;
      Mid(0)=(p(1)-d_v[0](1)-m1*p(1)+m2*d_v[0](0))/(m2-m1);
      Mid(1)=p(1)+m1*(Mid(0)-p(0));
      double distAB = point_distance(p,v_e[vert_inside[0]]);
      double distAM = point_distance(p,Mid);
      double distBM = point_distance(v_e[vert_inside[0]],Mid);
      if(distAM<distAB && distBM<distAB){
	return area_circ-area_cs[0]+area_t[0];
      }
      else
	return area_t[0]+area_cs[0];
    }
  }
  //###########################################################################
  //                    1 vertices and 3 edges
  //###########################################################################
  if(nv==1 && ne==3){
    unsigned int count_j=0;
    for(unsigned j=0; j<ne; j++){
      if(part_edge[edge_inside[j]][0]){
	d_v[count_j].zero();
	d_v[count_j].add(v_e[edge_inside[j]]);
	d_v[count_j].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][0]);
	count_j++;
      }
      if(part_edge[edge_inside[j]][1]){
	d_v[count_j].zero();
	d_v[count_j].add(v_e[edge_inside[j]]);
	d_v[count_j].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][1]);
	count_j++;
      }
    }
    for(unsigned j=0; j<2; j++){
      double distA,distB;
      int posA,posB;
      distA = distB = 1e5;
      posA = posB = 0;
      for(int i=0;i<4;i++){
	double dist_t = point_distance(v_e[vert_out[j]],d_v[i]);
	if(dist_t<distA){
	  distB=distA;
	  posB=posA;
	  distA = dist_t;
	  posA = i;
	}
	else if(dist_t<distB){
	  distB = dist_t;
	  posB = i;
	}
      }
      chord[j] = point_distance(d_v[posA],d_v[posB]);
      angle = 2.*asin(chord[j]/(2.*radius));
      area_cs[j] = 0.5*pow(radius,2)*(angle-sin(angle));
      area_t[j] = triangle_area(v_e[vert_out[j]],d_v[posA],d_v[posB]);
    }
    if(!cs_inside){return elem_volume-area_t[0]+area_cs[0]-area_t[1]+area_cs[1];}
    else{return elem_volume-area_t[0]+area_cs[0]-area_t[1]+area_cs[1];}
  }
  //###########################################################################
  //                    2 vertices and 2 edges
  //###########################################################################
  if(nv==2 && ne==2){
    for(unsigned j=0; j<ne; j++){
      d_v[j].zero();
      d_v[j].add(v_e[edge_inside[j]]);
      if(part_edge[edge_inside[j]][0]){d_v[j].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][0]);}
      else{                            d_v[j].add_scaled(d_a_e[edge_inside[j]],t_e[edge_inside[j]][1]);}
    }
    chord[0] = point_distance(d_v[0],d_v[1]);
    angle = 2.*asin(chord[0]/(2.*radius));
    area_cs[0] = 0.5*pow(radius,2)*(angle-sin(angle));
    area_t[0] = triangle_area(v_e[vert_out[0]],d_v[0],d_v[1]);
    if(!cs_inside){return elem_volume-area_t[0]+area_cs[0];}
    else{return elem_volume-area_t[0]+area_cs[0];}
  }
  //###########################################################################
  //                    3 vertices and 0 edges
  //###########################################################################
  if(nv==3 && ne==0){
    if(!cs_inside){return elem_volume;}
    else{return elem_volume;}
  }
  return 0.;
}

void initial_condition_ox(EquationSystems& es,const std::string& libmesh_dbg_var(system_name)){
  libmesh_assert_equal_to(system_name,"Nutrient");
  TransientLinearImplicitSystem & sys = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  es.parameters.set<Real>("time") = sys.time = 0;
  sys.project_solution(exact_value_ox,NULL,es.parameters);
}

void assemble_ox(EquationSystems& es,const std::string& libmesh_dbg_var(system_name)){
  libmesh_assert_equal_to (system_name, "Diffusion");
  const MeshBase& mesh = es.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  TransientLinearImplicitSystem & system = es.get_system<TransientLinearImplicitSystem>("Nutrient");
  const unsigned int v_nut = system.variable_number("nut_var");
  const Real time_size = es.parameters.get<Real>("time_step_size");
  const Real nut_d = es.parameters.get<Real>("nut_d");
  const Real con_b = es.parameters.get<Real>("con_b");
  const Real con_t = es.parameters.get<Real>("con_t");
  const Real con_n = es.parameters.get<Real>("con_n");
  const Real end_r = es.parameters.get<Real>("end_r");
  const int nut_coeff = es.parameters.get<Real>("nut_coeff");
  const Real ac_radius = es.parameters.get<Real>("action_radius");
  const Real h_max_msh = es.parameters.get<Real>("h_max_mesh");
  const Real height    = es.parameters.get<Real>("domain_diameter");
  list <Cell> *Cells_local = es.parameters.get<list <Cell> *>("l_cells");
  const double h_bin   = ac_radius+h_max_msh;
  const int number_bins = ceil(height/h_bin);
  const int total_bins  = number_bins*number_bins;
  //========== Generate bins ==========
  vector< list < Cell * > > Cell_Bins(total_bins);
  std::list<Cell>::iterator it;
  for(it = Cells_local->begin(); it != Cells_local->end(); ++it){
    int ix = floor((*it).x/h_bin);
    int jy = floor((*it).y/h_bin);
    int xy = ix+jy*number_bins;
    if((*it).x<0 || (*it).y<0 || (*it).x>height || (*it).y> height || jy>=number_bins || ix>=number_bins || xy >=total_bins){
      cout << "Error" << endl;
      cout << "Cell = ( " << (*it).x << " , " << (*it).y << " ) = ( " << ix << " , " << jy << " ) = " << xy << endl;
      cout << "State = " << (*it).state << endl;
      cout << "total_bins  = " << total_bins << endl;
      cout << "number_bins = " << number_bins << endl;
      cout << "h_bin       = " << h_bin << endl;
      getchar();
    }
    Cell_Bins[xy].push_back(&*it);
  }
  //========== Continue stnd ==========    
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);
  const std::vector<Real>& JxW = fe->get_JxW();
  const std::vector<std::vector<Real> >& phi = fe->get_phi();
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  DenseMatrix<Number> Ke;
  DenseVector<Number> Fe;
  std::vector<dof_id_type> dof_indices;
  std::vector<dof_id_type> dof_indices_nut;
  MeshBase::const_element_iterator       el     = mesh.active_local_elements_begin();
  const MeshBase::const_element_iterator end_el = mesh.active_local_elements_end();
  for ( ; el != end_el ; ++el){
    const Elem* elem = *el;
    //========== Computed volume fraction ==========
    Point e_center = elem->centroid();
    int ix = floor(e_center(0)/h_bin);
    int jy = floor(e_center(1)/h_bin);
    double p_t,p_n,p_e,p_tip,p_s,p_gs;
    double elem_volume = elem->volume();
    p_t = p_n = p_e = p_tip = p_s = p_gs = 0.;
    for(int xx = -1; xx<=1; xx++){
      if(ix>0 && ix+xx<number_bins){
	for(int yy = -1; yy<=1; yy++){
	  if(jy>0 && jy+yy<number_bins){
	    int bin_xy = (ix+xx)+(jy+yy)*number_bins;
	    std::list<Cell*>::iterator cell_ab;
	    for(cell_ab = Cell_Bins[bin_xy].begin(); cell_ab != Cell_Bins[bin_xy].end(); ++cell_ab){
	      if( (*(*cell_ab)).state == 0 || (*(*cell_ab)).state == 4) continue;
	      else{
		Point p( (*(*cell_ab)).x,(*(*cell_ab)).y,0.);
		double area_computed = 0.;
		if(elem->type()==3){
		  area_computed = area_inside_tri(p,elem->point(0),elem->point(1),elem->point(2),(*(*cell_ab)).C_radius,elem_volume);
		}
		else{
		  area_computed = area_inside_quad(p,elem->point(0),elem->point(1),elem->point(2),elem->point(3),(*(*cell_ab)).C_radius);
		}
		if(area_computed > 0){
		  if((*(*cell_ab)).state == 6) p_n+=(area_computed/elem_volume);
//		  else if((*(*cell_ab)).state >= 7 && (*(*cell_ab)).state != 13) p_e+=(area_computed/elem_volume);
		  else if((*(*cell_ab)).state >= 7 && (*(*cell_ab)).activate == 2) p_s+=(area_computed/elem_volume);
		  else if((*(*cell_ab)).state >= 7 && (*(*cell_ab)).prev_state == 7) p_e+=(area_computed/elem_volume);
                 else if( (*(*cell_ab)).state >= 0 && (*(*cell_ab)).state <= 6) p_t+=(area_computed/elem_volume);
		}
	      }
	    }
	  }
	}
      }
    }
    double p_b;
    if((p_t+p_n+p_e)>1.0){p_b = 0.;}
    else{p_b=1.0-(p_t+p_n+p_e);}
    //========== Continue stnd ==========  
    dof_map.dof_indices(elem,dof_indices);
    dof_map.dof_indices(elem,dof_indices_nut,v_nut);
    fe->reinit(elem);
    Ke.resize(dof_indices.size(),dof_indices.size());
    Fe.resize(dof_indices.size());
    for (unsigned int qp=0; qp<qrule.n_points(); qp++){
      Number nut_old = 0.0;
      for(unsigned int l=0; l<phi.size(); l++){
	nut_old += phi[l][qp]*system.old_solution(dof_indices_nut[l]);
      }
      for (unsigned int i=0; i<phi.size(); i++){
	Fe(i) += JxW[qp]*(nut_old
			  +time_size*nut_coeff*(p_e+p_tip+p_s+p_gs)*end_r
			  )*phi[i][qp];
	for (unsigned int j=0; j<phi.size(); j++){
	  Ke(i,j) += JxW[qp]*time_size*nut_d*dphi[i][qp]*dphi[j][qp];
	  Ke(i,j) += JxW[qp]*(1.0+time_size*(p_b*con_b+
					     p_n*con_n+
					     p_t*con_t+
					     nut_coeff*(p_e+p_tip+p_s+p_gs)*end_r))*phi[i][qp]*phi[j][qp];
	}
      }
    }
    dof_map.heterogenously_constrain_element_matrix_and_vector(Ke,Fe,dof_indices);
    system.matrix->add_matrix (Ke,dof_indices);
    system.rhs->add_vector    (Fe,dof_indices);
  }
}
