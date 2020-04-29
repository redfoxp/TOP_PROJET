#ifndef LBM_INIT_H
#define LBM_INIT_H

/********************  HEADERS  *********************/
#include "lbm_struct.h"
#include "lbm_comm.h"

/*******************  FUNCTION  *********************/
void init_cond_velocity_0_density_1(Mesh * mesh);
void setup_init_state_circle_obstacle(Mesh * mesh, lbm_mesh_type_t * mesh_type, const lbm_comm_t * mesh_comm);
void setup_init_state_global_poiseuille_profile(Mesh * mesh, lbm_mesh_type_t * mesh_type,const lbm_comm_t * mesh_comm);
void setup_init_state_border(Mesh * mesh, lbm_mesh_type_t * mesh_type, const lbm_comm_t * mesh_comm);
void setup_init_state(Mesh * mesh, lbm_mesh_type_t * mesh_type, const lbm_comm_t * mesh_comm);

#endif //LBM_INIT_H
