#ifndef GENERAL_LIB_H_
#define GENERAL_LIB_H_

#include <iostream>
#include <math.h>
#include <sstream>
#include <algorithm>
#include <list>
#include "libmesh/libmesh.h"
#include "libmesh/mesh.h"
#include "libmesh/serial_mesh.h"
#include "libmesh/mesh_generation.h"
#include "libmesh/equation_systems.h"
#include "libmesh/linear_implicit_system.h"
#include "libmesh/transient_system.h"
#include "libmesh/fe.h"
#include "libmesh/quadrature_gauss.h"
#include "libmesh/sparse_matrix.h"
#include "libmesh/numeric_vector.h"
#include "libmesh/dense_matrix.h"
#include "libmesh/dense_vector.h"
#include "libmesh/elem.h"
#include "libmesh/dof_map.h"
#include "libmesh/dirichlet_boundaries.h"
#include "libmesh/zero_function.h"
#include "libmesh/exodusII_io.h"
#include "libmesh/gnuplot_io.h"
#include "libmesh/getpot.h"
#include "nr3.h"
#include "ran.h"

using namespace libMesh;
using namespace std;

#endif
