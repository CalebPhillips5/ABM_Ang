#include "general_libraries.h"
#include "cell_abm.h"
#include "chrono"
#include "libmesh/fe_interface.h"
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
//		      12 - Waiting Stalk Cells
//		      13 - Deactivated Endothelial Cells
//		      14 - Activated through Anastomosis Endothelial Cells
//###########################################################################

void divide_tumor(list<Cell> &Cells_local,std::list<Cell>::iterator &it,Ran& ran,double height){
  double x = (*it).x,y = (*it).y;
  double rN = (*it).N_radius/sqrt(2.0);
  double nx,ny;
  unsigned int try_set_cell=1;
  while(try_set_cell){
    nx = rN*(2.0*ran.doub()-1.0);
    ny = sqrt(rN*rN - nx*nx);
    Point pos_n;
    pos_n(0) = x-nx;
    pos_n(1) = y-ny;
    Point pos_p;
    pos_p(0) = x+nx;
    pos_p(1) = y+ny;
    Point center;
    center(0) = 0.5*height;
    center(1) = 0.5*height;
    if(point_distance(pos_p,center)>0.5*height || point_distance(pos_n,center)>0.5*height){
      try_set_cell = 1;
    }
    else
      try_set_cell = 0;
  }
  (*it).set(x+nx,y+ny,rN,(*it).C_radius/sqrt(2.0),(*it).A_radius/sqrt(2.0),(*it).uptake,(*it).time,5);
  Cell a;
  a.set(x-nx,y-ny,rN,(*it).C_radius,(*it).A_radius,(*it).uptake,(*it).time,5);
  a.activate = 0;
  a.anasto = 0;
  Cells_local.insert (it,a);
  //Cells_local.push_back(a);
}

void divide_endo(list<Cell> &Cells_local,list<Cell*> &Cells_tip,std::list<Cell>::iterator &it){
  double x = (*it).x,y = (*it).y;
  double rN = (*it).N_radius/sqrt(2.0);
  double nx,ny;
  int tip_cell_location = 0;
  double minimum = 3.0*(*it).A_radius;
  std::list<Cell*>::iterator it_p;
  int counter = -1;
  for(it_p = Cells_tip.begin(); it_p != Cells_tip.end(); ++it_p){
    counter++;
    double d = sqrt(pow((*it).x-(**it_p).x,2)+pow((*it).y-(**it_p).y,2));
    if(d < minimum){
      minimum = d;
      tip_cell_location = counter;
    }
  }
  it_p = Cells_tip.begin();
  advance(it_p,tip_cell_location);
  nx = ((**it_p).x-(*it).x)/4;
  ny = ((**it_p).y-(*it).y)/4;
  Cell a;
  cout << "Prevstate = " << (*it).prev_state << endl;
  if((*it).prev_state==13){
    a.set(x-nx,y-ny,(*it).N_radius/sqrt(2.0),(*it).C_radius/sqrt(2.0),(*it).A_radius/sqrt(2.0),600,(*it).time,11);
    a.prev_state = (*it).prev_state;
    a.anasto = (*it).anasto;
    a.activate = (*it).activate;
    (*it).set(x+nx,y+ny,(*it).N_radius/sqrt(2.0),(*it).C_radius/sqrt(2.0),(*it).A_radius/sqrt(2.0),600,(*it).time,11);
    (*it).prev_state = a.prev_state;
    Cells_local.insert (it,a);
  }
  else{
    a.set(x+nx,y+ny,(*it).N_radius/sqrt(2.0),(*it).C_radius/sqrt(2.0),(*it).A_radius/sqrt(2.0),600,(*it).time,11);
    a.prev_state = (*it).prev_state;
    a.anasto = (*it).anasto;
    a.activate = (*it).activate;
    (*it).set(x-nx,y-ny,(*it).N_radius/sqrt(2.0),(*it).C_radius/sqrt(2.0),(*it).A_radius/sqrt(2.0),600,(*it).time,11);
    (*it).prev_state = a.prev_state;
    Cells_local.insert (it,a);
  }
  //Cells_local.push_back(a);
}

void update_states(list<Cell> &Cells_local,list<Cell*> &Cells_tip,list<Cell> &Cells_vessel_start,list<Cell> &Cells_vessel_end,int time,Ran& ran,EquationSystems& equation_systems,const int outside_cells){
  const Real N_radius       = equation_systems.parameters.get<Real>("nucleus_radius");
  const Real C_radius       = equation_systems.parameters.get<Real>("cell_radius");
  const Real A_radius       = equation_systems.parameters.get<Real>("action_radius");
  const Real tau_A          = equation_systems.parameters.get<Real>("apop_time");
  const Real tau_P          = equation_systems.parameters.get<Real>("cellc_time");
  const bool deactivation   = equation_systems.parameters.get<bool>("deactivation");
  const bool anastomosis    = equation_systems.parameters.get<bool>("anastomosis");
  const Real tau_E          = tau_P + 20.0;
  const Real tau_g1         = equation_systems.parameters.get<Real>("g1_time");
  const Real tau_c          = equation_systems.parameters.get<Real>("calc_time");
  const Real delta_tt       = equation_systems.parameters.get<Real>("delta_tt");
  const double max_c_radius = equation_systems.parameters.get<Real> ("cell_radius");
  const Real alfa_p_barra   = equation_systems.parameters.get<Real>("prol_intes");
  const Real alfa_a         = equation_systems.parameters.get<Real>("apop_intes");
  const Real tau_NL         = equation_systems.parameters.get<Real>("lysing_time");
  const Real f_NS           = equation_systems.parameters.get<Real>("f_NS");
  const Real sigma_h        = equation_systems.parameters.get<Real>("hypoxic_thrs");
  const Real sigma_n        = 0.25*sigma_h;
  const Real N_out_max      = equation_systems.parameters.get<Real>("max_out_cells");
  const double height       = equation_systems.parameters.get<Real>("domain_diameter");
  const double vegf_ths     = equation_systems.parameters.get<Real>("vegf_ths");
  const double tip_distance = equation_systems.parameters.get<Real>("tip_distance");
  const int max_tip_cells   = equation_systems.parameters.get<int>("max_tip_cells");
  vector<int> remove_cells;
  std::vector<Number> soln;  
  equation_systems.build_solution_vector(soln);
  const MeshBase& mesh = equation_systems.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  //===== VEGF =====//
  TransientLinearImplicitSystem & system = equation_systems.get_system<TransientLinearImplicitSystem>("VEGF");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);
  //===== Nutrient =====//
  TransientLinearImplicitSystem & nut_sys = equation_systems.get_system<TransientLinearImplicitSystem>("Nutrient");
  const DofMap& sys_dof_map = nut_sys.get_dof_map();
  FEType sys_fe_type = sys_dof_map.variable_type(0);
  UniquePtr<FEBase> sys_fe (FEBase::build(dim, sys_fe_type));
  QGauss sys_qrule (dim, sys_fe_type.default_quadrature_order());
  sys_fe->attach_quadrature_rule (&sys_qrule);
  std::vector<dof_id_type> dof_indices;
  //========== Loop over cells ==========
  std::list<Cell>::iterator it;
  int i_count = -1;
  for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
    i_count++;  
    double sigma_cel = 0.6;
    double nut_cel = 0.6;
    double tau = time - (*it).time;
//cout << "which cell we are on: " << i_count << endl;
    //========== Compute vegf concentration ==========
    if((*it).state == 7){
      const Point p((*it).x,(*it).y,0.);
      double dpo=sqrt(pow(p(0)-500.0,2)+pow(p(1)-500.0,2));
      if(dpo<500.0){
	// Get point locator structure
	UniquePtr<PointLocatorBase> point_locator = mesh.sub_point_locator();
	// Get element where point p is
	const Elem* top_elem = (*point_locator)(p);
	// Get the location (on the reference element) of the point p located in physical space
	Point p_master = FEInterface::inverse_map(dim, fe_type, top_elem, p);
	dof_map.dof_indices (top_elem, dof_indices);
	fe->reinit (top_elem);
	const unsigned int dof_size = dof_indices.size();
	std::vector<Real> point_phi(dof_size);
	// Get shape function
	for (unsigned int i=0; i != dof_size; i++)
	  point_phi[i] = FEInterface::shape(dim, fe_type, top_elem, i, p_master);
	// Compute solution at point p
	Number vegf_p = 0.0;
	for(unsigned int l=0; l<dof_size; l++){
	  vegf_p += point_phi[l]*system.old_solution(dof_indices[l]);
	}
	sigma_cel = vegf_p;
      }
      else{
	sigma_cel = 0;
      } 
      //cout << "size of active_local_eles = " << count << endl;
// =============================================================================
    //========== Endothelial Cell Growth ==========
    if((*it).N_radius < N_radius/2){
      if(tau < tau_E){ // I believe tau_E is 38. 
	if(1+ (10.0+tau_g1+(tau -tau_E-10))/(10.0+tau_g1) > 0){
        //(*it).N_radius = ((N_radius/2)/sqrt(2.0))*sqrt(1+ (10.0+tau_g1+(tau -tau_E))/(10.0+tau_g1));//make endo cells grow fast, it may make up the distance
	(*it).N_radius = ((N_radius/2)/sqrt(2.0))*sqrt(1+ (10.0+tau_g1+(tau -tau_E-10))/(10.0+tau_g1)); // -10 from tau_E
        (*it).C_radius = ((C_radius/2)/sqrt(2.0))*sqrt(1+ (10.0+tau_g1+(tau -tau_E-10))/(10.0+tau_g1));
        (*it).A_radius = (*it).C_radius*1.214;
	}
      }
      else{
        (*it).N_radius = N_radius/2;
        (*it).C_radius = C_radius/2;
        (*it).A_radius = A_radius/2;
      }
      continue;
    }
//========== Endo-Tip Transition =========
      double minimum = 100.;
      double d = 0.;
      if(Cells_tip.size() < max_tip_cells )//Just one tip cell TRYING TO STUDY THIS
        if(sigma_cel >= vegf_ths){
          if(sqrt(pow((*it).x-0.5*height,2) + pow((*it).y-0.5*height,2)) < 0.5*height-60){ // if not on boundary
            if(Cells_tip.size() > 0){
              std::list<Cell*>::iterator it_p;
              int counter = -1;
              for(it_p = Cells_tip.begin(); it_p != Cells_tip.end(); ++it_p){
                counter++;
                d = sqrt(pow((*it).x-(**it_p).x,2)+pow((*it).y-(**it_p).y,2));
                if(d < minimum) minimum = d;
              }
              if(minimum > tip_distance*A_radius){
                //========== if no branching then uncomment this if
                //if(Cells_local[i].prev_state ==7){
                (*it).state = 10;
		(*it).anasto = 0;
		(*it).activate = 0;
		(*it).time = 0;
                Cells_tip.push_back(&(*it));
                //}
              } // end if distance
            } // end if Cell_tip.size()>0
            else{
              (*it).state = 10;
	      (*it).anasto = 0;
	      (*it).activate = 0;
		(*it).time = 0;
              Cells_tip.push_back(&(*it));
            } // end else
          } // end if not on boundary
        }
    } // end if state == 7


    //========== Compute nutrient concentration ==========
    if((*it).state == 1 || (*it).state == 2 || (*it).state == 3){
auto start2 = std::chrono::high_resolution_clock::now();
      Point p((*it).x,(*it).y,0.);
      // Get point locator structure
      UniquePtr<PointLocatorBase> point_locator = mesh.sub_point_locator();
      // Get element where point p is
      const Elem* top_elem = (*point_locator)(p);
      // Get the location (on the reference element) of the point p located in physical space
      Point p_master = FEInterface::inverse_map(dim, sys_fe_type, top_elem, p);
      sys_dof_map.dof_indices (top_elem, dof_indices);
      sys_fe->reinit (top_elem);
      const unsigned int dof_size = dof_indices.size();
      std::vector<Real> point_phi(dof_size);
      // Get shape function
      for (unsigned int i=0; i != dof_size; i++)
	point_phi[i] = FEInterface::shape(dim, sys_fe_type, top_elem, i, p_master);
      // Compute solution at point p
      Number vegf_p = 0.0;
      for(unsigned int l=0; l<dof_size; l++){
	vegf_p += point_phi[l]*nut_sys.old_solution(dof_indices[l]);
      }
      nut_cel = vegf_p;
    //========== Skip healthy cells ==========
//    if((*it).state == 6) continue;
    //========== Transition to hypoxic ==========
    if(((*it).state == 1 || (*it).state == 2) && nut_cel < sigma_h){
      (*it).time = time;
      (*it).prev_state = (*it).state;
      (*it).state = 3;
    }
    //========== Transition from hypoxic ==========
    else if((*it).state == 3){
      if(nut_cel > sigma_h){ // become proliferative
	(*it).state = 1;
      }
      else if(nut_cel < sigma_n){ // become necrotic
	(*it).state = 0;
	(*it).time = time;
	(*it).prev_state = (*it).state;
      }
    }
    //========== Tumor cells ==========
    if((*it).state == 1){
    double alfa_p = alfa_p_barra*(nut_cel-sigma_h)/(1 - sigma_h)*(1-outside_cells/N_out_max);
      double prob_mit = (1.0-exp(-delta_tt*alfa_p))/2.0;
      double prob_apo = ((1.0-exp(-delta_tt*alfa_a))/2.0)+prob_mit;
      double rand_num = ran.doub();
      //========== Mitosis ==========
      if(prob_mit >= rand_num){
        (*it).time = time;
        (*it).prev_state = (*it).state;
        (*it).state = 2;
        continue;
      }
      //========== Apoptosis ==========
      if(prob_apo >= rand_num){
        (*it).time = time;
        (*it).prev_state = (*it).state;
        (*it).state = 4;
        continue;
      }
      continue;
    }
    //========== Proliferative cells ==========
    if((*it).state == 2){
      if(tau >= (tau_P-tau_g1) ){
        divide_tumor(Cells_local,it,ran,height);
      }
      continue;
    }

   } // end if ========Compute nutrient concentration ======= 


    //========== Increase calcification ==========   changed else if to if
    if((*it).state == 0 && (*it).cal != 1){
      if(tau < tau_NL){
	(*it).C_radius = C_radius * sqrt(1+(f_NS*tau)/tau_NL);
	(*it).N_radius = N_radius * sqrt(1+(f_NS*tau)/tau_NL);
	//(*it).A_radius = 1.214*(*it).C_radius;
	(*it).cal = tau/tau_c;
	continue;
      }
      if(tau == tau_NL){
	(*it).C_radius = C_radius * pow(1.0+f_NS, 1.0/3.0);
	(*it).N_radius = N_radius * pow(1.0+f_NS, 1.0/3.0);
	//(*it).A_radius = 1.214*(*it).C_radius; 
	(*it).cal = tau/tau_c;
	continue;
      }
      if(tau > tau_NL && tau <= tau_c){
	if((tau - tau_NL)< 10){
	  (*it).C_radius = (*it).C_radius - (C_radius-N_radius)/10.0;
	  (*it).N_radius = (*it).N_radius - (C_radius-N_radius)/10.0;
	}
	else{
	  (*it).C_radius = N_radius;
	  (*it).N_radius = 0.0;
	}
	//(*it).A_radius = 1.214*(*it).C_radius;  
	(*it).cal = tau/tau_c;
	continue;
      }
      if (tau >= tau_c) (*it).cal = 1.0;
      continue;
    }
    
	//========== Tip Cell Time ==========
    if((*it).state == 8){
      (*it).time += 1;
    }
    //========== Proliferative stalk cells ==========
    if((*it).state == 9){
      if(tau >= (tau_E-tau_g1)){
	divide_endo(Cells_local,Cells_tip,it);
//	(*it).time = time; //UNCOMMENTED THIS
      }
      continue;
    }
//
    if((*it).state == 12){
     cout << "LOOKING AT STATE 12 tau " << tau << endl;  
	if(tau >= 38){ // 46+tau_E-tau_g1)){
	   (*it).state = 9;
	   (*it).time = time;
	}
    }

    //============= Tumor Cell============
    if((*it).state == 5){
      if(tau < tau_P){
	(*it).N_radius = (N_radius/sqrt(2.0))*sqrt(1+ (tau_g1+(tau -tau_P))/tau_g1);
	(*it).C_radius = (C_radius/sqrt(2.0))*sqrt(1+ (tau_g1+(tau -tau_P))/tau_g1);
	(*it).A_radius = (*it).C_radius*1.214;
      }
      else{
	(*it).prev_state = 2;
	(*it).N_radius = N_radius;
	(*it).C_radius = C_radius;
	(*it).A_radius = A_radius;
	(*it).state = 1;
	(*it).time = time;
      }
      continue;
    }
    //========== Stalk Cell Growth ==========
    if((*it).state == 11){
      if(tau < tau_E){
	if(1+ (10.0+tau_g1+(tau -tau_E-10))/(10.0+tau_g1) > 0){
	(*it).N_radius = ((N_radius/2)/sqrt(2.0))*sqrt(1+ (10.0+tau_g1+(tau -tau_E-10))/(10.0+tau_g1)); // -10 from tau_E
	(*it).C_radius = ((C_radius/2)/sqrt(2.0))*sqrt(1+ (10.0+tau_g1+(tau -tau_E-10))/(10.0+tau_g1));
//	cout << "stalk cell time " << (*it).time << endl;
//	cout << "time " << time << endl;
//	cout << "tau " << tau << endl;
	}
      }
      else{
	(*it).N_radius = N_radius/2;
	(*it).C_radius = C_radius/2;
	(*it).A_radius = A_radius/2;
	(*it).state = 9;
	(*it).time = time;
cout << "Cell was full grown!" << endl;
      }
      continue;
    }
    //========== Update apoptotic size ==========
    if((*it).state == 4){
      if(tau >= tau_A){
	remove_cells.push_back(i_count);
      }
      else{ 
	if((tau_A - tau) == 1)  (*it).cal = (*it).C_radius/tau_A;
	(*it).C_radius = (*it).C_radius - (*it).cal;
	(*it).N_radius = (*it).N_radius - 0.5*(*it).cal;
      }
      continue;
    }
  }
  for(unsigned int i = 0; i < remove_cells.size(); i++){
    it = Cells_local.begin();
    advance(it,remove_cells[i]-i);
    Cells_local.erase(it);
  }

 if(deactivation == 0){ 
  
  // ERNESTO HERE IS WHERE I TRY TO IMPLEMENT THE DEACTIVATION.
  //  Two flags, first to signify disconnect one way,
  //  Second to signify disconnect the other.
  // If both are 1, deactivate. 
  // Reset all flags for the timestep.
  // I ran into some errors and gave up for the night. (I could use the sleep).

//HERE IS THE BEGINNING OF TRYING TO IMPLEMENT PREV_STATE 15 ITERATORS
//I think the best way is to just iterate from prev_state 15 -> prev_state 15 until the end. 
Cells_vessel_start.clear();

  std::vector<int> count;
  std::list<Cell>::iterator iter_15;
int i = -1;
for (iter_15 = Cells_local.begin(); iter_15 != Cells_local.end(); iter_15++){
i++;
	if( (*iter_15).prev_state == 15){
	count.push_back(i);
      Cells_vessel_start.push_back(*iter_15);
}
}

  std::list<Cell>::iterator iter_bgn;
  std::list<Cell>::iterator iter_bgn2;
  std::list<Cell>::iterator iter_end;
  std::list<Cell>::iterator iter_end2;
  std::list<Cell>::iterator iter_kill;
if(Cells_vessel_start.size() > 1){
for(int i=0; i < Cells_vessel_start.size()-1; i++){
	iter_bgn2 = Cells_local.begin();
	iter_end2 = Cells_local.begin();
//	advance(iter_bgn2,count[i]);
//	advance(iter_end2,count[i+1]);
	advance(iter_bgn2,count[i]-1);
	advance(iter_end2,count[i+1]+1);
//	cout << "count 1, count 2 : " << count[i] << " " << count[i+1] << endl;
  bool one_break = 0;
  bool deactivate_cells = 0;
        int broken;
int counter = -1;

  for(iter_bgn = iter_bgn2; iter_bgn != iter_end2; ++iter_bgn){
    counter++;
    if( ((*iter_bgn).state==7 && (*iter_bgn).prev_state!=7) || ((*iter_bgn).state>=8 && (*iter_bgn).state!=13) ){
      std::list<Cell>::iterator previous_cell;
      previous_cell = iter_bgn;
      previous_cell--;
      double dist = sqrt(pow((*previous_cell).x-(*iter_bgn).x,2)+pow((*previous_cell).y-(*iter_bgn).y,2));
      if(dist > 3.0*(*iter_bgn).A_radius){
        one_break=1;
 //       cout << "broken1 = " << counter << endl;
        break;
      }
    }
  }// for iter_bgn
  if(one_break){
        int counted = -1;
    for(iter_end = iter_end2; iter_end != iter_bgn; --iter_end){
        counted++;
      if( ((*iter_end).state==7 && (*iter_end).prev_state!=7) || ((*iter_end).state>=8 && (*iter_end).state!=13) ){
        std::list<Cell>::iterator previous_cell;
        previous_cell = iter_end;
        previous_cell--;
        double dist = sqrt(pow((*previous_cell).x-(*iter_end).x,2)+pow((*previous_cell).y-(*iter_end).y,2));
        if(dist > 1.55*max_c_radius){
          deactivate_cells=1;
          broken = counted;
  //        cout << "broken2 = " << counted << endl;
	  break;
        }
      }
    }
  } //if one break
  if(deactivate_cells){
    std::list<Cell>::iterator iter_kill2;
        // THIS IS WHERE CALEB TRIES HIS MAGIC

//		cout << "in deactivate cells if" << endl;
        for(iter_kill2 = iter_bgn; iter_kill2 != iter_end; iter_kill2++){
//		cout << "deactivating cells" << endl;
                (*iter_kill2).state=13;
                (*iter_kill2).activate=0;
        }
   }//if deactivate

}
 }
} // END IF DEACTIVATION

              std::list<Cell*>::iterator it_p;
              for(it_p = Cells_tip.begin(); it_p != Cells_tip.end(); ++it_p){
                 if( (**it_p).state != 8 && (**it_p).state != 10 ){
                        Cells_tip.erase(it_p);
                 } else{
	            continue;
                 }
              }

   if(anastomosis == 0){
    Cells_vessel_start.clear();
    int k = -1;
    std::list<Cell>::iterator iter_activate;
    std::vector<int> count;
	// this gives us the location of all activated tip cells
	for(iter_activate = Cells_local.begin(); iter_activate != Cells_local.end(); iter_activate++){
        k++;
		if( (*iter_activate).activate == 2){
			(*iter_activate).activate = 0;
		}
	  if( (*iter_activate).activate == 1){
             count.push_back(k); 
	     Cells_vessel_start.push_back(*iter_activate); // this gives us a list of all anasto activated guys.
	  } // end if
        } // end for
  std::list<Cell>::iterator iter_bgn;
  std::list<Cell>::iterator iter_end;
  std::list<Cell>::iterator iter_bgn2;
  std::list<Cell>::iterator iter_end2;
  // Here is the idea. We want to check all the activated anastomosis cells. The ones that go together will 
  // be right next to each other. So we could just check the activated cells that are closest together,
  // and then iterate between them to activate. Then we don't have to keep track of which cells are paired up,
  // it will be obvious. Now to code that.... 
if(Cells_vessel_start.size() > 1){
  for(int i = 0; i < Cells_vessel_start.size()-1; i++){
     iter_bgn2 = Cells_local.begin(); 
     iter_end2 = Cells_local.begin(); 
		  // advance to the ith activated guy and the i+1 activated.
     advance(iter_bgn2,count[i]);
		  // you can't just check the next one, you have to check all of them.
    for(int j = 1; j < Cells_vessel_start.size()-i; j++){
	iter_end2 = Cells_local.begin();
		  // 
        advance(iter_end2,count[i+j]);
        double dist = sqrt( pow((*iter_bgn2).x-(*iter_end2).x,2) + pow( (*iter_bgn2).y - (*iter_end2).y,2 ) );
       if (dist < 3*(*iter_bgn2).A_radius ){
		  // activate the cells between them. But leave them at activate == 1
int count_deactivated = 0;
         for( iter_bgn = iter_bgn2; iter_bgn != iter_end2; ++iter_bgn) {
                 if ( (*iter_bgn).activate == 0 ){
			(*iter_bgn).activate = 2;
 //                       cout << "hey man we activated cells!" << endl;
		 }
                 if ( (*iter_bgn).state == 13 ){
			count_deactivated += 1;
		 }
		  
         } // end for iter_bgn 
// if theres deactivated guys in there, deactivate all of those cells. just throw activate back to 0. That's all.
      if(count_deactivated > 0){
// we still have the right i and j. so that's money.
	iter_bgn2 = Cells_local.begin();
	iter_end2 = Cells_local.begin();
        advance(iter_bgn2,count[i]);
        advance(iter_end2,count[i+j]);
         for( iter_bgn = iter_bgn2; iter_bgn != iter_end2; ++iter_bgn) {
		if( (*iter_bgn).activate == 2 )
	             (*iter_bgn).activate = 0;	
	 }
	} // end if count_deactivated
       }  // end if
    } // end for j
  } // end for i
} // if cell.size() > 1

   }// end if anastomosis


}

void normal(double& N_x,double& N_y,double x_cel,double y_cel,double height){
  if (x_cel < height/2 && x_cel != 0){
    //coeficientes da reta que passa pela célula e o centro da meia circunferência	
    double a = (height/2 - y_cel)/(height/2 - x_cel);
    double b = y_cel - a*x_cel;
    //calcula o ponto que intersepta a meia circunferência com a reta acima (origem da normal)
    //Obs: sinal da raiz é negativo pois o interesse é no menor x (a esquerda do ponto central)
    double p_x = (-2*a*b+a*height+height - sqrt(-4*a*b*height+2*a*height*height-4*b*b+4*b*height))/(2*a*a+2);
    double p_y = a*p_x + b;
    //calcular vetor normal vezes distância
    N_x = x_cel - p_x;
    N_y = y_cel - p_y;
  }
  if (x_cel > height/2){
    //coeficientes da reta que passa pela célula e o centro da meia circunferência	
    double a = (height/2 - y_cel)/(height/2 - x_cel);
    double b = y_cel - a*x_cel;
    //calcula o ponto que intersepta a meia circunferência com a reta acima (origem da normal)
    //Obs: sinal da raiz é positivo pois o interesse é no maior x (a direita do ponto central)
    double p_x = (-2*a*b+a*height+height + sqrt(-4*a*b*height+2*a*height*height-4*b*b+4*b*height))/(2*a*a+2);
    double p_y = a*p_x + b;
    //calcular vetor normal vezes distância
    N_x = x_cel - p_x;
    N_y = y_cel - p_y;
  }
  if (x_cel == height/2){
    //calcular vetor normal vezes distância
    N_x = 0;
    if (y_cel > height/2)  N_y = y_cel-height;
    else   N_y = y_cel;
  }
}

void compute_forces(list<Cell> &Cells_local,list<Cell*> &Cells_tip,list<Cell> &Cells_vessel_start,list<Cell> &Cells_vessel_end,EquationSystems& equation_systems,double height,int& outside_cells,int time){
  //========== Define parameters ==========
  const double c_ccr              = equation_systems.parameters.get<Real> ("c_ccr");
  const double c_eea              = equation_systems.parameters.get<Real> ("c_eea");
  const bool unbreakable          = equation_systems.parameters.get<bool> ("unbreakable");
  const double repulsion_distance = equation_systems.parameters.get<Real> ("repulsion_distance");
  const double max_c_radius       = equation_systems.parameters.get<Real> ("cell_radius");
  const double h_bin              = 3.0*max_c_radius;
  const unsigned int number_bins  = ceil(height/h_bin);
  const unsigned int total_bins   = number_bins*number_bins;
  const double c_cvg              = 2.0*c_ccr;
  //========== Build vegf solution ==========
  std::vector<Number> soln;
  equation_systems.build_solution_vector(soln);
  const MeshBase& mesh = equation_systems.get_mesh();
  const unsigned int dim = mesh.mesh_dimension();
  TransientLinearImplicitSystem & system = equation_systems.get_system<TransientLinearImplicitSystem>("VEGF");
  const DofMap& dof_map = system.get_dof_map();
  FEType fe_type = dof_map.variable_type(0);
  UniquePtr<FEBase> fe (FEBase::build(dim, fe_type));
  QGauss qrule (dim, fe_type.default_quadrature_order());
  fe->attach_quadrature_rule (&qrule);
  const std::vector<std::vector<RealGradient> >& dphi = fe->get_dphi();
  std::vector<dof_id_type> dof_indices;
  //========== Loop until converge ==========
  unsigned int unchecked_tip = 1;
  double iter = 0.0;
  while (iter < 500){
    //========== Generate bins ==========
    vector< list < Cell * > > Cell_Bins(total_bins);
    std::list<Cell>::iterator it;
    for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
      unsigned int ix = floor((*it).x/h_bin);
      unsigned int jy = floor((*it).y/h_bin);
      unsigned int xy = ix+jy*number_bins;
      if(jy>=number_bins || ix>=number_bins || xy >=total_bins){
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
    int n =1; int m = 1; double M = 1;
    double phi_x,phi_y,psi_x,psi_y;
    //========== Cell-cell force ==========
    for(unsigned int bin_i = 0; bin_i < Cell_Bins.size(); bin_i++){
      std::list<Cell*>::iterator cell_a;
      for(cell_a = Cell_Bins[bin_i].begin(); cell_a != Cell_Bins[bin_i].end(); ++cell_a){
	unsigned int jy = floor(bin_i/number_bins);
	unsigned int ix = bin_i%number_bins;
	unsigned int bin_j;
	if(unchecked_tip){
	  //===== If Activated Tip Cell =====*==//
	  if((*(*cell_a)).state==10){
	    int prev_st = 13;
	    for(int xx=-1; xx<=1 ; xx++){
	      if(ix+xx<number_bins && ix>0){
		for(int yy=-1; yy<=1 ; yy++){
		  if(jy+yy<number_bins && jy>0){
		    bin_j = (ix+xx)+(jy+yy)*number_bins;
		    std::list<Cell*>::iterator cell_b;
		    for(cell_b = Cell_Bins[bin_j].begin(); cell_b != Cell_Bins[bin_j].end(); ++cell_b){
		      double d;
			int count = -1;

//HERE IS WHERE I DID THE NEW WAYS TO SELECT STALK AND PREV_STATE.
//GAZE UPON THE BEAUTY
			std::list<Cell>::iterator caleb;
			for(caleb = Cells_local.begin(); caleb != Cells_local.end(); ++caleb){
					count++;
				if( (*caleb).state == 10 ){
					if( (*caleb).state == 10){
					(*caleb).state = 8;
					(*caleb).prev_state = 10;
					}

					advance(caleb,-2);
				//	Cells_vessel_start.push_back(*caleb);
					cout << "Cells_vessel_start.size() = " << Cells_vessel_start.size() << endl;
					if( (*caleb).x > 90){
					(*caleb).prev_state = 15;	
	                                }
					advance(caleb,1);	
					(*caleb).prev_state = 13;
					(*caleb).state = 9;
					(*caleb).time = time;
				
					advance(caleb,2);
					(*caleb).prev_state = 14;
					(*caleb).state = 9;
					(*caleb).time = time;
	
					advance(caleb,1);
					Cells_vessel_end.push_back((*caleb));
					if( (*caleb).x > 90){
					(*caleb).prev_state = 15;
					}
				}
			}
		      //=======* If Endothelial Cell =====**==//
		      if((*(*cell_b)).state==9){
			d = sqrt(pow((*(*cell_a)).x-(*(*cell_b)).x,2)+pow((*(*cell_a)).y-(*(*cell_b)).y,2));
			if(d < ((*(*cell_a)).A_radius+(*(*cell_b)).A_radius)){ // && (*(*cell_b)).prev_state==(*(*cell_a)).prev_state){
				if((*(*cell_b)).x > (*(*cell_a)).x && (*(*cell_b)).x >= 150  ){ // if stalk is closer to tumor skip first divide
					(*(*cell_b)).state = 12;
					cout << "this one is closer; endo " << (*(*cell_b)).x << " tip " << (*(*cell_a)).x << endl;
				} //else {
//			  (*(*cell_b)).state=9;
//				}
//			  (*(*cell_b)).time=time;
//			  (*(*cell_b)).prev_state=prev_st;
//			  prev_st=14;
                        }
			
			//	if(d > (*(*cell_a)).A_radius+(*(*cell_b)).A_radius && d < 2.0*((*(*cell_a)).A_radius+(*(*cell_b)).A_radius) && (*(*cell_b)).prev_state==(*(*cell_a)).prev_state){
			// if((*(*cell_b)).prev_state != 7)
//			    (*(*cell_b)).prev_state = 15;
			//	}
		      } // THIS ENDS THE ENDOTHELIAL CELL LOOP
		    } // end for cell_b
		  } // end if yy
		} // end for yy
	      } // end if xx
	    } // end for xx
	    // Change Activated Tip Cell to Tip Cell
	    (*(*cell_a)).state=8;
	    (*(*cell_a)).prev_state=10;
	  } // end if activated
	  //========== If Stalk Cell ==========
	  if((*(*cell_a)).state==9 || (*(*cell_a)).state==11 || (*(*cell_a)).state==12   ){
	    double minimum = 100;
	    double d = 0;
	    for(int xx=-1; xx<=1 ; xx++){
	      if(ix+xx<number_bins && ix>0){
		for(int yy=-1; yy<=1 ; yy++){
		  if(jy+yy<number_bins && jy>0){
		    bin_j = (ix+xx)+(jy+yy)*number_bins;
		    std::list<Cell*>::iterator cell_b;
		    for(cell_b = Cell_Bins[bin_j].begin(); cell_b != Cell_Bins[bin_j].end(); ++cell_b){
		      //========== If Tip Cell ==========
		      if((*(*cell_b)).state==8 || (*(*cell_b)).state==10){
			d = sqrt(pow((*(*cell_a)).x-(*(*cell_b)).x,2)+pow((*(*cell_a)).y-(*(*cell_b)).y,2));
			if(d < minimum) minimum = d;
		      }
		    } //end for cell_b
		  } // end if yy 
		} // end for yy
	      } // end if xx
	    } // end for xx
	    if(minimum > 1.55*max_c_radius){
	      (*(*cell_a)).state=7;
	      //(*(*cell_a)).prev_state=15;
	   } else 
	      continue;
	  } // end if stalk cell
 	}
	if((*(*cell_a)).state==8 || (*(*cell_a)).state==9 || (*(*cell_a)).state==11){
	  Point p((*(*cell_a)).x,(*(*cell_a)).y,0.);
	  // Get point locator structure
	  UniquePtr<PointLocatorBase> point_locator = mesh.sub_point_locator();
	  // Get element where point p is
	  const Elem* top_elem = (*point_locator)(p);
	  // Get the location (on the reference element) of the point p located in physical space
	  Point p_master = FEInterface::inverse_map(dim, fe_type, top_elem, p);
	  dof_map.dof_indices (top_elem, dof_indices);
	  fe->reinit (top_elem);
	  Gradient gradient_vegf;
	  double magnitude_vegf;
	  unsigned int qp=0;	
	  for(unsigned int l=0; l<dof_indices.size(); l++){
	    gradient_vegf.add_scaled(dphi[l][qp], system.old_solution(dof_indices[l]));
	  }
	  magnitude_vegf = sqrt(pow(gradient_vegf(0),2)+pow(gradient_vegf(1),2));
	  //======== Removed 0.8 scale ==========
		if( (*(*cell_a)).anasto == 0 ){
	  (*(*cell_a)).prev_F[0] = gradient_vegf(0)/magnitude_vegf;
	  (*(*cell_a)).prev_F[1] = gradient_vegf(1)/magnitude_vegf;
		}
			if( (*(*cell_a)).anasto == 1){
//cout << "(*(cell_a)).prev_F[0],F[1] = " << (*(*cell_a)).prev_F[0] << " " << (*(*cell_a)).prev_F[1] << endl;
}
	}
	if((*(*cell_a)).state!=15){
	  //========== Compute Forces ==========
	  int vy[5]={0, 1,1,1,0};
	  int vx[5]={0,-1,0,1,1};
	  for(unsigned int pxy = 0; pxy < 5; pxy++){
	    if(jy+vy[pxy]<number_bins && jy+vy[pxy]>=0){
	      if(ix+vx[pxy]<number_bins && ix+vx[pxy]>=0){
		bin_j = (ix+vx[pxy])+(jy+vy[pxy])*number_bins;
		if(bin_j<bin_i){cout << "Error in bin selection" << endl;getchar();}
		std::list<Cell*>::iterator cell_start = Cell_Bins[bin_j].begin();
		if(pxy==0) cell_start = next(cell_a);
		std::list<Cell*>::iterator cell_b;
		for(cell_b = cell_start; cell_b != Cell_Bins[bin_j].end(); ++cell_b){
		  double F_cca_i[2]={0.0,0.0},F_ccr_i[2]={0.0,0.0};
		  double r_x = (*(*cell_b)).x - (*(*cell_a)).x;
		  double r_y = (*(*cell_b)).y - (*(*cell_a)).y;
		  double R_A = (*(*cell_b)).A_radius + (*(*cell_a)).A_radius;
		  double R_N = (*(*cell_b)).N_radius + (*(*cell_a)).N_radius;
		  double R   = (*(*cell_b)).C_radius + (*(*cell_a)).C_radius;
		  potential_adh(phi_x,phi_y,r_x,r_y,R_A,n);
		  potential_rep(psi_x,psi_y,r_x,r_y,R_N,R,M,m);
		  if((*(*cell_a)).prev_state == 13 && (*(*cell_b)).prev_state == 14)
		    potential_rep(psi_x,psi_y,r_x,r_y,R_N,repulsion_distance*R,M,m);
		  if((*(*cell_a)).prev_state == 14 && (*(*cell_b)).prev_state == 13)
		    potential_rep(psi_x,psi_y,r_x,r_y,R_N,repulsion_distance*R,M,m); 
		  double c_cca = 0.0;
		  double repulsion = 0.0;
		  if((*(*cell_a)).state > 6 && (*(*cell_b)).state > 6){
		    repulsion = c_ccr;
		    c_cca = c_eea;
		  }
		  else if((*(*cell_a)).state == 6 || (*(*cell_b)).state == 6){
		    repulsion = c_ccr;
		    c_cca = 0;
		  }
		  else if((*(*cell_a)).state < 6 || (*(*cell_b)).state < 6){
		    repulsion = c_ccr;
		    c_cca = 0;
		  }		  
		  F_ccr_i[0] += -repulsion*psi_x;
		  F_ccr_i[1] += -repulsion*psi_y;
		  if(
		     ((*(*cell_a)).prev_state==13 && (*(*cell_b)).prev_state==14) ||
		     ((*(*cell_a)).prev_state==14 && (*(*cell_b)).prev_state==13)
		     ){ // Lummen repulsion
		    (*(*cell_a)).F[0] += 1.0*F_ccr_i[0];//changing to 1.0 * Caleb 1/26
		    (*(*cell_a)).F[1] += 1.0*F_ccr_i[1];
		    (*(*cell_b)).F[0] -= 1.0*F_ccr_i[0];
		    (*(*cell_b)).F[1] -= 1.0*F_ccr_i[1];
		  }
		  else if(
			  ((*(*cell_a)).state == 8 && (*(*cell_b)).state == 7) ||
			  ((*(*cell_a)).state == 7 && (*(*cell_b)).state == 8)
			  ){ // No tip-endothelial adhesion
		    F_cca_i[0] += -c_cca*phi_x;
		    F_cca_i[1] += -c_cca*phi_y;
		    double test = 0.75;
/*		    (*(*cell_a)).F[0] -= test*(F_cca_i[0]); //CHANGED TO MAKE 5 TIMES, HOPING TO FIX INWARD GROWTH 1/24
		    (*(*cell_a)).F[1] -= test*(F_cca_i[1]);
		    (*(*cell_b)).F[0] += test*(F_cca_i[0]);
		    (*(*cell_b)).F[1] += test*(F_cca_i[1]);
*/
		double distance;
		double d = sqrt(pow((*(*cell_a)).x-(*(*cell_b)).x,2) + pow((*(*cell_a)).y-(*(*cell_b)).y,2));
		double hold = 0.0;
		if( (*(*cell_a)).C_radius >= (*(*cell_b)).C_radius ){
			hold = (*(*cell_a)).C_radius;
		}
		if( (*(*cell_b)).C_radius >= (*(*cell_a)).C_radius ){
			hold = (*(*cell_b)).C_radius;
		}
		if( d < 2.25*hold ){
			(*(*cell_a)).state = 7;
			(*(*cell_b)).state = 7;
//This would mean that both cells are deactivated, these are the droids we're looking for!
			//if( (*(*cell_a)).anasto == 1 || (*(*cell_b)).anasto == 1) {
			if( (*(*cell_a)).time >= 200 && (*(*cell_b)).time >= 200)  {
cout << "HEY MAN WE ANASTOD SOME CELLS HERE THEY JUST GOT DONE ACTIVATED!@@@@@@@@@@@@@@@@@" << endl;
			(*(*cell_a)).activate = 1;
			(*(*cell_b)).activate = 1;
			}
//What if we just activate here...? This never really gets hit more than once also.
//Since we are aiming to check this every timestep, to see if it's still working... Maybe we just supplement the deactivation algo?
//Or we do something VERY similar and add them to 




		} // end if d <
		  }
		  else if(
			  ((*(*cell_a)).state == 8 && (*(*cell_b)).prev_state == 15) || // THIS NEVER RUNS BECAUSE PREV_STATE == 15 MEANS STATE == 7 WHICH RUNS ABOVE
			  ((*(*cell_a)).prev_state == 15 && (*(*cell_b)).state == 8)
			  ){
		    F_cca_i[0] += -c_cca*phi_x;
		    F_cca_i[1] += -c_cca*phi_y;
		    double test = 0.75;
		   // double test = 1;
		    (*(*cell_a)).F[0] -= test*(F_cca_i[0]); //CHANGED TO MAKE 5 TIMES, HOPING TO FIX INWARD GROWTH 1/24
		    (*(*cell_a)).F[1] -= test*(F_cca_i[1]);
		    (*(*cell_b)).F[0] += test*(F_cca_i[0]);
		    (*(*cell_b)).F[1] += test*(F_cca_i[1]);
		  }
                  else if(
                          ((*(*cell_a)).state == 8 && ((*(*cell_b)).state == 9 || (*(*cell_b)).state == 11 )) || // i've got some questions about this for sure. 6/21
                          (((*(*cell_a)).state == 11  || (*(*cell_b)).state == 9 )&& (*(*cell_b)).state == 8)
                          ){
                    F_cca_i[0] += -c_cca*phi_x;
                    F_cca_i[1] += -c_cca*phi_y;

                    double test;
		if( (*(*cell_a)).state == 11  || (*(*cell_b)).state == 11 ){
                    test = 0.2;
		}
		else{
		   test = 0.05;
		}
			if( (*(*cell_a)).state ==8 ){
                    (*(*cell_a)).F[0] += test*(F_cca_i[0]); //CHANGED TO MAKE 5 TIMES, HOPING TO FIX INWARD GROWTH 1/24
                    (*(*cell_a)).F[1] += test*(F_cca_i[1]);
			}
			if( (*(*cell_b)).state ==8 ){
                    (*(*cell_b)).F[0] -= test*(F_cca_i[0]);
                    (*(*cell_b)).F[1] -= test*(F_cca_i[1]);
			}
                  }
		  else if((*(*cell_a)).state == 8 && (*(*cell_b)).state == 8){
		    double d = sqrt(pow((*(*cell_a)).x-(*(*cell_b)).x,2) + pow((*(*cell_a)).y-(*(*cell_b)).y,2));
		    if( d < 2*(*(*cell_b)).A_radius){
		      F_cca_i[0] += -c_cca*phi_x;
		      F_cca_i[1] += -c_cca*phi_y;
		      (*(*cell_a)).F[0] -= (F_cca_i[0]);
		      (*(*cell_a)).F[1] -= (F_cca_i[1]);
		      (*(*cell_b)).F[0] += (F_cca_i[0]);
		      (*(*cell_b)).F[1] += (F_cca_i[1]);
		      (*(*cell_a)).state = 7;
		      (*(*cell_b)).state = 7;
		      (*(*cell_a)).prev_state = 15;
		      (*(*cell_b)).prev_state = 15;
		    }
		  }
            else if( (((*(*cell_a)).state > 6 && (*(*cell_b)).state < 6)  || ((*(*cell_a)).state < 6 && (*(*cell_b)).state > 6)) && unbreakable == 1){
                      F_cca_i[0] += -c_cca*phi_x;
                      F_cca_i[1] += -c_cca*phi_y;

                        if( (*(*cell_a)).state > 6 ) {
                      (*(*cell_a)).F[0] -= (F_ccr_i[0]);
                      (*(*cell_a)).F[1] -= (F_ccr_i[1]);
                        }
                        if( (*(*cell_b)).state > 6 ) {
                      (*(*cell_b)).F[0] += (F_ccr_i[0]);
                      (*(*cell_b)).F[1] += (F_ccr_i[1]);
                        }
                  }
		  else { // if its not 13 and 14 theres adhesion
		    F_cca_i[0] += -c_cca*phi_x;
		    F_cca_i[1] += -c_cca*phi_y;
		  }
		  (*(*cell_a)).F[0] += (F_cca_i[0] + F_ccr_i[0]);
		  (*(*cell_a)).F[1] += (F_cca_i[1] + F_ccr_i[1]);
		  (*(*cell_b)).F[0] -= (F_cca_i[0] + F_ccr_i[0]);
		  (*(*cell_b)).F[1] -= (F_cca_i[1] + F_ccr_i[1]);


		  if( ( (*(*cell_a)).state==7 && (*(*cell_b)).state==9 ) || ( (*(*cell_b)).state == 7 && (*(*cell_a)).state==9) ){//ernesto
		    /*
			cout << " forces = " << (F_cca_i[0] + F_ccr_i[0]) << endl;
			cout << "statea,x,y = " << (*(*cell_a)).state << "," <<  (*(*cell_a)).x << "," <<(*(*cell_a)).y << endl;
			cout << "stateb,x,y = " << (*(*cell_b)).state << "," <<  (*(*cell_b)).x << "," <<(*(*cell_b)).y << endl;
			cout << "R_A = " << R_A << " R = " << R << endl;
		    */
		}


		  if((*(*cell_a)).state==8 && ((*(*cell_b)).state==12 || (*(*cell_b)).state == 9 || (*(*cell_b)).state==11     )){//ernesto
		    double psi_vx = 0.0;
		    double psi_vy = 0.0;
			if( (*(*cell_a)).anasto == 0 ){
		    potential_vegf(psi_vx,psi_vy,r_x,r_y,(*(*cell_a)).prev_F[0],(*(*cell_a)).prev_F[1],R_N,R,M,m);
		    (*(*cell_a)).F[0] += 0.6*c_cvg*psi_vx;//changed back to 1 1/27
		    (*(*cell_a)).F[1] += 0.6*c_cvg*psi_vy;
		if( (*(*cell_b)).state == 12 ){
		    (*(*cell_b)).F[0] += 0.05*c_cvg*psi_vx; // CALEB CHANGED 2/4
		    (*(*cell_b)).F[1] += 0.05*c_cvg*psi_vy;
		}
			} //end anasto
//it's gotta go here or it'll break. Make it work.

			if( (*(*cell_a)).anasto == 1 ){
		    potential_vegf(psi_vx,psi_vy,r_x,r_y,(*(*cell_a)).prev_F[0],(*(*cell_a)).prev_F[1],R_N,R,M,m);
		    (*(*cell_a)).F[0] += 0.6*c_cvg*psi_vx;//changed back to 1 1/27
		    (*(*cell_a)).F[1] += 0.6*c_cvg*psi_vy;
		if( (*(*cell_b)).state == 12 ){
		    (*(*cell_b)).F[0] += 0.05*c_cvg*psi_vx; // CALEB CHANGED 2/4
		    (*(*cell_b)).F[1] += 0.05*c_cvg*psi_vy;
		}
			}


		  }
		  else if((*(*cell_b)).state==8 && ((*(*cell_a)).state==9 || (*(*cell_a)).state==11 || (*(*cell_a)).state == 12)){
		    double psi_vx = 0.0;
		    double psi_vy = 0.0;
			if( (*(*cell_b)).anasto == 0 ){
		    potential_vegf(psi_vx,psi_vy,r_x,r_y,(*(*cell_b)).prev_F[0],(*(*cell_b)).prev_F[1],R_N,R,M,m);
		    (*(*cell_b)).F[0] += 0.6*c_cvg*psi_vx;//changed back to 1 1/27
		    (*(*cell_b)).F[1] += 0.6*c_cvg*psi_vy;
		if( (*(*cell_a)).state == 12 ){
		    (*(*cell_a)).F[0] += 0.05*c_cvg*psi_vx;
		    (*(*cell_a)).F[1] += 0.05*c_cvg*psi_vy;
		} // end if state == 12
			} // end anasto

			if( (*(*cell_b)).anasto == 1 ){
		    potential_vegf(psi_vx,psi_vy,r_x,r_y,(*(*cell_b)).prev_F[0],(*(*cell_b)).prev_F[1],R_N,R,M,m);
		    (*(*cell_b)).F[0] += 0.6*c_cvg*psi_vx;//changed back to 1 1/27
		    (*(*cell_b)).F[1] += 0.6*c_cvg*psi_vy;
		if( (*(*cell_a)).state == 12 ){
		    (*(*cell_a)).F[0] += 0.05*c_cvg*psi_vx;
		    (*(*cell_a)).F[1] += 0.05*c_cvg*psi_vy;
		}
			}





		  } // else if state == 8 
		}
	      }
	    }
	  }
	  //========== Cell-boundary force ==========
	  /*
	    double F_ct[2]={0.0,0.0},F_rct[2]={0.0,0.0};
	    double N_x,N_y;
	    double K;
	    if(k_var) K = outside_cells/max_outside;
	    else K = k_val;
	    double c_ct = 10.0*K;
	    double c_rct = 4.88836*K;
	    normal(N_x,N_y,(*(*cell_a)).x,(*(*cell_a)).y,height);
	    potential_adh(phi_x,phi_y,N_x,N_y,(*(*cell_a)).A_radius,n);
	    potential_rep(psi_x,psi_y,N_x,N_y,(*(*cell_a)).N_radius,(*(*cell_a)).C_radius,M,m);
	    F_ct[0] = -c_ct*phi_x;
	    F_ct[1] = -c_ct*phi_y;
	    F_rct[0] = -c_rct*psi_x;
	    F_rct[1] = -c_rct*psi_y;
	    (*(*cell_a)).F[0] += (F_ct[0] + F_rct[0]);
	    (*(*cell_a)).F[1] += (F_ct[1] + F_rct[1]);  
	  */
	}
      }
    }
//Anastomosis actually might have to go right here. Since it can't go above because cell_b is within the surrounding bins and we need to check them all. Well, we need to check all the endothelial cells. Is that even possible? Everything is possible.
			std::list<Cell*>::iterator tip;
 for(tip = Cells_tip.begin(); tip != Cells_tip.end(); ++tip){
    if ( (**tip).anasto == 1 ){ // if already anastomosing
// 		    (**tip).F[0] -= 0.008*(**tip).prev_F[0];//changed back to 1 1/27
//		    (**tip).F[1] -= 0.008*(**tip).prev_F[1];//changed back to 1 1/27
//			cout << "F[0], F[1],x,y = " << (**tip).F[0] << "," << (**tip).F[1] << "," << (**tip).x << "," << (**tip).y << endl;
	break;
    } else if ( (**tip).anasto == 0 ){
			std::list<Cell>::iterator endo;
      for(endo = Cells_local.begin(); endo != Cells_local.end(); ++endo){
	if( (*endo).state > 6 ) {
		int pos_x,pos_y;
		pos_x = (**tip).x + 150*(**tip).prev_F[0];
		pos_y = (**tip).y + 150*(**tip).prev_F[1];
            if( (*endo).x < pos_x + 50 && (*endo).x > pos_x - 50 && (*endo).y > pos_y - 50 && (*endo).y < pos_y +50 ){
		(**tip).anasto = 1;
		cout << "tip position x,y  = " << (**tip).x << "," << (**tip).y << "####################################" << endl;
		cout << "endo position x,y  = " << (*endo).x << "," << (*endo).y << "####################################" << endl;
		break;
	    } // end if
	} // end if state 
      } // end for
    } // end else if
 } // end for

    //========== Compute cell speed ==========
    double v_max= 5.0e-3;
    for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
      (*it).v[0] = -0.5*((*it).F[0]);
      (*it).v[1] = -0.5*((*it).F[1]);
      double max_speed = euclid_norm((*it).v[0],(*it).v[1]);
      if (max_speed >= v_max) v_max = max_speed;
    }
    double delta_tt = (1.0/v_max);
    vector<int> remove_cells;
    //========== Update cell position ==========
    for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
      if((*it).state != 7){
	(*it).x += delta_tt * (*it).v[0];
	(*it).y += delta_tt * (*it).v[1];
      }
      else if((*it).prev_state==13 || (*it).prev_state==14  || (*it).prev_state ==15 ){
	(*it).x += delta_tt * (*it).v[0];
	(*it).y += delta_tt * (*it).v[1];
      }
if( (*it).state == 8 || (*it).state == 9 || (*it).state == 11 || (*it).prev_state > 12 ){
  //cout << "(*it).state = " << (*it).state << " x = " << (*it).x << " y = " << (*it).y << " A_radius = " << (*it).A_radius <<" F[0] = " << (*it).F[0] << " F[1] = " << (*it).F[1] << endl;
}
      (*it).F[0] = 0.0;
      (*it).F[1] = 0.0;
      double dist = distance((*it),height/2.0,height/2.0);
      if(dist >= height/2.0){
	//remove_cells.push_back(i);
	outside_cells++;
      }
    }
    //========== Remove outside cells ==========
    /*
    for(unsigned int i = 0; i < remove_cells.size(); i++){
      Cells_local.erase (Cells_local.begin() + remove_cells[i]-i); 
    }
    */
    iter+= delta_tt;
    unchecked_tip = 0;
  }
}

double euclid_norm(double x,double y){return sqrt(pow(x,2.0) + pow(y,2.0));}

void potential_adh(double& phi_x,double& phi_y,double r_x,double r_y,double R_A,int n){
  double r = euclid_norm(r_x,r_y);
  if((r >= 0) && (r <= R_A)){
    double var = pow((1.0 - r/R_A),n+1)/r;
    phi_x = var*r_x;
    phi_y = var*r_y;
  }
  else{
    phi_x = 0.0;
    phi_y = 0.0;
  }
}

void potential_rep(double& psi_x,double& psi_y,double r_x,double r_y,double R_N,double R,double M,int m){
  double r = euclid_norm(r_x,r_y);
  psi_x = 0;
  psi_y = 0;
  if((r >= 0) && (r < R_N)){
    double c = pow((1.0 - R_N/R),m+1) - M;
    double var = -(c*r/R_N + M)/r;
    psi_x = var*r_x;
    psi_y = var*r_y;
  }
  else if((r >= R_N) && (r <= R)){
    double var = -pow((1.0 - (r/R)),m+1)/r;
    psi_x = var*r_x;
    psi_y = var*r_y;
  }
  else if(r > R){
    psi_x = 0.0;
    psi_y = 0.0;
  }
}

void potential_vegf(double& psi_x,double& psi_y,double r_x,double r_y,double vegf_x,double vegf_y,double R_N,double R,double M,int m){
  double r = euclid_norm(r_x,r_y);
  if((r >= 0) && (r < R_N)){
    double c = pow((1.0 - R_N/R),m+1) - M;
    double var = -(c*r/R_N + M)/r;
    psi_x = var*vegf_x;
    psi_y = var*vegf_y;
  }
  else if((r >= R_N) && (r <= R)){
    double var = -pow((1.0 - (r/R)),m+1)/r;
    psi_x = var*vegf_x;
    psi_y = var*vegf_y;
  }
  else if(r > R){
    psi_x = 0;
    psi_y = 0;   
  }
}

//========== Restart save file ==========
void restart_save(list<Cell> &Cells_local,string name_s){
  ofstream out_file;
  out_file.open(name_s);
  std::list<Cell>::iterator it;
  for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
    out_file << (*it).x << " " << (*it).y << " " << (*it).N_radius << " " << (*it).C_radius << " " << (*it).A_radius << " " << (*it).uptake << " " << (*it).cal << " " << (*it).time << " " << (*it).state << " " << (*it).prev_state << " " << (*it).v[0] << " " << (*it).v[1] << " " << (*it).F[0] <<" " << (*it).F[1] << " " << (*it).prev_F[0] << " " << (*it).prev_F[1] << " " << (*it).anasto << " " << (*it).activate << endl;
  }
  out_file.close();
}

//========== Reading save file ==========
void restart_function(list<Cell> &Cells_local,list<Cell*> &Cells_tip,string name_s){
  double x,y,N_radius,C_radius,A_radius,uptake,cal,v[2],F[2],prev_F[2];
  int time,state,prev_state,anasto,activate;
  int i_count = -1;
  std::list<Cell>::iterator it;
  std::ifstream read;
  std::string line;
  read.open(name_s);
  if(!read.is_open())
    std::cout << "Error opening.\n";
  while(std::getline(read,line)){
    std::istringstream(line) >> x >> y >> N_radius >> C_radius >> A_radius >> uptake >> cal >> time >> state >> prev_state >> v[0] >> v[1] >> F[0] >> F[1] >> prev_F[0] >> prev_F[1] >> anasto >> activate;
    Cell a;
    a.set(x,y,N_radius,C_radius,A_radius,uptake,time,state);
    a.prev_state = prev_state;
    a.v[0] = v[0];
    a.v[1] = v[1];
    a.prev_F[0] = prev_F[0];
    a.prev_F[1] = prev_F[1];
    a.anasto = anasto;
    a.activate = activate;
    Cells_local.push_back(a);
    i_count++;
    if(i_count==0)it = Cells_local.begin();
    else ++it;
    if(a.state==8){
      Cells_tip.push_back(&(*it));
    }
  }
  read.close();
}

void init_cond_cells(list<Cell> &Cells_local,const Parameters& all_parameters,Ran& ran){
  const Real n_radius        = all_parameters.get<Real>("nucleus_radius");
  const Real c_radius        = all_parameters.get<Real>("cell_radius");
  const Real a_radius        = all_parameters.get<Real>("action_radius");
  const Real domain_diameter = all_parameters.get<Real>("domain_diameter");
  const Real lambda_cell     = all_parameters.get<Real>("lambda_cell");
  const Real initial_con     = all_parameters.get<Real>("initial_con");
  const int  total_t         = all_parameters.get<int> ("initial_tumor");
  const int ic_type          = all_parameters.get<int> ("ic_type");
  const Real tau_P           = all_parameters.get<Real>("cellc_time");
  const Real tau_g1          = all_parameters.get<Real>("g1_time");
  double pos[2];
  int number_t = 0;
  //int number_h = 0;
  double center[2];
  center[0] = 0.5*domain_diameter;
  center[1] = 0.5*domain_diameter;
  double confluence = 0.0;
  if(ic_type==3){
    //============*Initialize Tumor Cells==========**==//
    for(int i = 0; i<5; i++){
      for(int j = 0; j<5; j++){
	pos[0] = 600.0+2.0*c_radius*j;
	pos[1] = 500.0+2.0*c_radius*(i-2);
	double dist = sqrt(pow(pos[0]-center[0],2)+pow(pos[1]-center[1],2)); 
	if(dist < 0.5*domain_diameter){
	  Cell a;
	  a.set(pos[0],pos[1],n_radius,c_radius,a_radius,0.1,0,1);
	  if(ran.doub() < 0.5){
	    a.time = -ran.doub()*(tau_P-tau_g1);
	    a.prev_state = 1;
	    a.state = 2;
	    a.anasto = 0;
	    a.activate = 0;
	  }
	  Cells_local.push_back(a);
	}
      }
    }
    //============*Initialize Endothelial Cells=====*==// ERNESTO
    for(int j = 0; j<2; j++){
      for(int i = 0; i<180; i++){
	pos[0] = 80.0+20*j;
	pos[1] = 10+9.0*i;
	double dist = sqrt(pow(pos[0]-center[0],2)+pow(pos[1]-center[1],2)); 
	if(dist < 0.5*domain_diameter){
	  Cell a;
	  a.set(pos[0],pos[1],n_radius/2,c_radius/2,a_radius/2,0.1,0,7);
	  a.prev_state = 7;
	  a.closer = 0;
	  a.activate = 0;
	  a.anasto = 0;
	  Cells_local.push_back(a);
	}
      }
    }
  }
  else if(ic_type==2){
    while(confluence < initial_con){
      Point pos;
      unsigned int place_cell = 1;
      pos(0) = ran.doub()*domain_diameter;
      pos(1) = ran.doub()*domain_diameter;
      Point center;
      center(0) = 0.5*domain_diameter;
      center(1) = 0.5*domain_diameter;
      if(point_distance(pos,center)>0.5*domain_diameter){
	place_cell = 0;
      }
      else{
	std::list<Cell>::iterator it;
	for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
	  Point cell;
	  cell(0) = (*it).x;
	  cell(1) = (*it).y;
	  if(point_distance(pos,cell)<2.*n_radius){
	    place_cell = 0;
	    break;
	  }
	}
      }
      if(place_cell){
	Cell a;
	a.set(pos(0),pos(1),n_radius,c_radius,a_radius,lambda_cell,0,1);
	if(ran.doub() < 0.5){
	  a.time = -ran.doub()*(tau_P-tau_g1);
	  a.prev_state = 1;
	  a.state = 2;
	}
	Cells_local.push_back(a);
	confluence += ((M_PI*std::pow(c_radius,2))/(M_PI*std::pow(0.5*domain_diameter,2)));
      }
    }
  }
  else{
    while(number_t < total_t){
      Point pos;
      unsigned int place_cell = 1;
      pos(0) = ran.doub()*domain_diameter;
      pos(1) = ran.doub()*domain_diameter;
      Point center;
      center(0) = 0.5*domain_diameter;
      center(1) = 0.5*domain_diameter;
      if(point_distance(pos,center)>0.5*domain_diameter){
	place_cell = 0;
      }
      else{
	std::list<Cell>::iterator it;
	for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
	  Point cell;
	  cell(0) = (*it).x;
	  cell(1) = (*it).y;
	  if(point_distance(pos,cell)<2.*n_radius){
	    place_cell = 0;
	    break;
	  }
	}
      }
      if(place_cell){
	Cell a;
	a.set(pos(0),pos(1),n_radius,c_radius,a_radius,lambda_cell,0,1);
	if(ran.doub() < 0.5){
	  a.time = -ran.doub()*(tau_P-tau_g1);
	  a.prev_state = 1;
	  a.state = 2;
	}
	Cells_local.push_back(a);
	number_t++;
      }
    }
  }
}

//========== Latex order save file ==========
void print_order(list<Cell> &Cells_local,double domain_diameter,string s,int file_number,int t){
  //========== Creating a string with the correct name ==========//
  const char *c = s.c_str();
  char n[100],name[200];
  sprintf(n,"%d_%05d.txt",file_number,t);
  strcpy(name,c);
  strcat(name,n);
  stringstream ss;
  string name_s;
  ss << name;
  ss >> name_s;
  //========== Saving the data ==========//
  ofstream out_file;
  out_file.open (name_s);
  out_file << domain_diameter << " " << domain_diameter << endl;
  //out_file << Cells_local.size() << endl;
  out_file << Cells_local.size() << " " << t << endl;
  //out_file << 0 << " " << 0 << endl;
  std::list<Cell>::iterator it;
  for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
    if((*it).state==13){
      out_file << 27 << endl;
    }
    else{
      if( (*it).activate == 1 || (*it).activate == 2)
	out_file << 1227 << endl;
      else if( (*it).prev_state==13 || (*it).prev_state==14 || (*it).prev_state ==15 )
	out_file << (*it).prev_state << endl;
      else	
	out_file << (*it).state << endl;
    }
    out_file << scientific;
    out_file << (*it).x << " " << (*it).y << endl;
    out_file << (*it).N_radius << " " << (*it).C_radius << " " << (*it).cal << endl;
    out_file << (*it).v[0] << " " << (*it).v[1] << endl;
  }  
  out_file.close();
}

//========== Latex save file ==========
void save_latex(list<Cell> &Cells_local,double domain_diameter,string s,int file_number,int t){
  //========== Creating a string with the correct name ==========//
  const char *c = s.c_str();
  char n[100],name[200];
  sprintf(n,"%d_%05d.txt",file_number,t);
  strcpy(name,c);
  strcat(name,n);
  stringstream ss;
  string name_s;
  ss << name;
  ss >> name_s;
  //========== Saving the data ==========//
  ofstream out_file;
  out_file.open (name_s);
  out_file << domain_diameter << " " << domain_diameter << endl;
  out_file << Cells_local.size() << " " << t << endl;
  out_file << 0 << " " << 0 << endl;
  std::list<Cell>::iterator it;
  for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
    out_file << (*it).state << endl;
    out_file << scientific;
    out_file << (*it).x << " " << (*it).y << endl;
    out_file << (*it).N_radius << " " << (*it).C_radius << " " << (*it).cal << endl;
    out_file << (*it).v[0] << " " << (*it).v[1] << endl;
  }  
  out_file.close();
}

//========== Matlab save file ==========
void save_cells(list<Cell> &Cells_local,double domain_diameter,string s,int file_number,int t){
  //========== Creating a string with the correct name ==========//
  const char *c = s.c_str();
  char n[100],name[200];
  sprintf(n,"%d_%05d.m",file_number,t);
  strcpy(name,c);
  strcat(name,n);
  stringstream ss;
  string name_s;
  ss << name;
  ss >> name_s;
  //========== Saving the data ==========//
  ofstream out_file;
  out_file.open (name_s);
  out_file << "cells = zeros(" << 2*Cells_local.size()+1 << "," << 4 << ");" <<endl;
  out_file << "cells = [" << -1 << " " << domain_diameter*0.5 << " " << domain_diameter*0.5 << " " << domain_diameter*0.5 << endl;
  std::list<Cell>::iterator it;
  for(it = Cells_local.begin(); it != Cells_local.end(); ++it){
    out_file << (*it).state << " " << (*it).x << " " << (*it).y << " " << (*it).C_radius << endl;
    out_file << 0 << " " << (*it).x << " " << (*it).y << " " << (*it).N_radius << endl;
  }
  out_file << "];";
  out_file.close();
}

void Cell::print(){
  cout << "Position = ( " << x << " , " << y << " )" << endl;
  cout << "Radius | Nuclear = " << N_radius << ", Cell = " << C_radius << ", Max = " << A_radius << endl;
  cout << "Time = " << time << endl;
  cout << "V = (" << v[0] << " , " << v[1] << " )" << endl;
 
  switch (state) {
    case 0:
      cout << "State = Dead" << endl;
      break;
    case 1:
      cout << "State = Alive" << endl;
      break;
    case 2:
      cout << "State = Proliferative" << endl;
      break;
    case 3:
      cout << "State = Hypoxic" << endl;
      break;
    case 4:
      cout << "State = Dying" << endl;
      break;
    case 5:
      cout << "State = G1" << endl;
      break;
    case 6:
      cout << "State = Normoxic" << endl;
      break;
    default:
      cout << "State = " << state << endl;
  }
}

void Cell::set(double X,double Y,double RN,double R,double RA,double Uptake, int Time,int State){
  x = X;
  y = Y;
  N_radius = RN;
  C_radius = R;
  A_radius = RA;
  uptake = Uptake;
  time = Time;
  state = State;
  v[0] = 0.0;
  v[1] = 0.0;
  F[0] = 0.0;
  F[1] = 0.0;
  cal  = 0.0;
  prev_state = -1;
}
