/********************  HEADERS  *********************/
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include <unistd.h>
#include "lbm_comm.h"

/*******************  FUNCTION  *********************/
int lbm_helper_pgcd(int a, int b)
{
	int c;
	while(b!=0)
	{
		c = a % b;
		a = b;
		b = c;
	}
	return a;
}

/*******************  FUNCTION  *********************/
/**
 * Affiche la configuation du lbm_comm pour un rank donné
 * @param mesh_comm Configuration à afficher
**/
void  lbm_comm_print( lbm_comm_t *mesh_comm )
{
	int rank ;
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	printf( " RANK %d ( LEFT %d RIGHT %d TOP %d BOTTOM %d CORNER %d, %d, %d, %d ) ( POSITION %d %d ) (WH %d %d ) \n", rank,
									    mesh_comm->left_id,
									    mesh_comm->right_id,
										mesh_comm->top_id,
									    mesh_comm->bottom_id,
										mesh_comm->corner_id[0],
		 								mesh_comm->corner_id[1],
		 								mesh_comm->corner_id[2],
		 		 						mesh_comm->corner_id[3],
									    mesh_comm->x,
									    mesh_comm->y,
									    mesh_comm->width,
									    mesh_comm->height );
}

/*******************  FUNCTION  *********************/
int helper_get_rank_id(int nb_x,int nb_y,int rank_x,int rank_y)
{
	if (rank_x < 0 || rank_x >= nb_x)
		return -1;
	else if (rank_y < 0 || rank_y >= nb_y)
		return -1;
	else
		return (rank_x + rank_y * nb_x);
}

/*******************  FUNCTION  *********************/
/**
 * Initialise un lbm_comm :
 * - Voisins
 * - Taille du maillage local
 * - Position relative
 * @param mesh_comm MeshComm à initialiser
 * @param rank Rank demandant l'initalisation
 * @param comm_size Taille totale du communicateur
 * @param width largeur du maillage
 * @param height hauteur du maillage
**/
void lbm_comm_init( lbm_comm_t * mesh_comm, int rank, int comm_size, int width, int height )
{
	//vars
	int nb_x;
	int nb_y;
	int rank_x;
	int rank_y;

	//compute splitting
	nb_y = lbm_helper_pgcd(comm_size,width);
	nb_x = comm_size / nb_y;
	//nb_x = lbm_helper_pgcd(comm_size,width);
	//nb_y = comm_size / nb_x;

	//check
	assert(nb_x * nb_y == comm_size);
	if (height % nb_y != 0)
		fatal("Can't get a 2D cut for current problem size and number of processes.");

	//calc current rank position (ID)
	rank_x = rank % nb_x;
	rank_y = rank / nb_x;

	//setup nb
	mesh_comm->nb_x = nb_x;
	mesh_comm->nb_y = nb_y;

	//setup size (+2 for ghost cells on border)
	mesh_comm->width = width / nb_x + 2;
	mesh_comm->height = height / nb_y + 2;

	//setup position
	mesh_comm->x = rank_x * width / nb_x;
	mesh_comm->y = rank_y * height / nb_y;

	// Compute neighbour nodes id
	mesh_comm->left_id  = helper_get_rank_id(nb_x,nb_y,rank_x - 1,rank_y);
	mesh_comm->right_id = helper_get_rank_id(nb_x,nb_y,rank_x + 1,rank_y);
	mesh_comm->top_id = helper_get_rank_id(nb_x,nb_y,rank_x,rank_y - 1);
	mesh_comm->bottom_id = helper_get_rank_id(nb_x,nb_y,rank_x,rank_y + 1);
	mesh_comm->corner_id[CORNER_TOP_LEFT] = helper_get_rank_id(nb_x,nb_y,rank_x - 1,rank_y - 1);
	mesh_comm->corner_id[CORNER_TOP_RIGHT] = helper_get_rank_id(nb_x,nb_y,rank_x + 1,rank_y - 1);
	mesh_comm->corner_id[CORNER_BOTTOM_LEFT] = helper_get_rank_id(nb_x,nb_y,rank_x - 1,rank_y + 1);
	mesh_comm->corner_id[CORNER_BOTTOM_RIGHT] = helper_get_rank_id(nb_x,nb_y,rank_x + 1,rank_y + 1);

	//if more than 1 on y, need transmission buffer
	if (nb_y > 1)
	{
		mesh_comm->buffer = malloc(sizeof(double) * DIRECTIONS * width / nb_x);
	} else {
		mesh_comm->buffer = NULL;
	}

	//if debug print comm
	#ifndef NDEBUG
	lbm_comm_print( mesh_comm );
	#endif
}


/*******************  FUNCTION  *********************/
/**
 * Libere un lbm_comm
 * @param mesh_comm MeshComm à liberer
**/
void lbm_comm_release( lbm_comm_t * mesh_comm )
{
	mesh_comm->x = 0;
	mesh_comm->y = 0;
	mesh_comm->width = 0;
	mesh_comm->height = 0;
	mesh_comm->right_id = -1;
	mesh_comm->left_id = -1;
	if (mesh_comm->buffer != NULL)
		free(mesh_comm->buffer);
	mesh_comm->buffer = NULL;
}

/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh_to_process Mesh a utiliser lors de l'échange des mailles fantomes
**/
void lbm_comm_sync_ghosts_horizontal( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int x )
{
	//vars
	//MPI_Status status;
	MPI_Request request[mesh->height-2];
	int flag = 0;

	//if target is -1, no comm
	if (target_rank == -1)
		return;

	int y;

	switch (comm_type)
	{
		case COMM_SEND:
			for( y = 0 ; y < mesh->height-2 ; y++ )
			{
				//MPI_Send( &Mesh_get_col( mesh_to_process, x )[y], DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
				MPI_Isend( &Mesh_get_col( mesh_to_process, x )[y], DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &request[y]);
			}

			break;
		case COMM_RECV:
			for( y = 0 ; y < mesh->height-2 ; y++ )
			{
				//MPI_Recv(  &Mesh_get_col( mesh_to_process, x )[y], DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,&status);
				MPI_Irecv( &Mesh_get_col( mesh_to_process, x )[y], DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &request[y]);
			}

			while (!flag) {
				MPI_Testall(mesh->height-2,request, &flag,MPI_STATUSES_IGNORE);
			}
			break;
		default:
			fatal("Unknown type of communication.");
	}
}

/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh_to_process Mesh a utiliser lors de l'échange des mailles fantomes
**/
void lbm_comm_sync_ghosts_diagonal( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int x ,int y)
{
	//vars
	MPI_Status status;
	//MPI_Request request;

	//if target is -1, no comm
	if (target_rank == -1)
		return;

	switch (comm_type)
	{
		case COMM_SEND:
			MPI_Send( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
			//MPI_Isend( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &request);

			break;
		case COMM_RECV:
			MPI_Recv( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &status);
			//MPI_Irecv( Mesh_get_cell( mesh_to_process, x, y ), DIRECTIONS, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD, &request);
			//MPI_Wait(&request, MPI_STATUS_IGNORE);
			break;
		default:
			fatal("Unknown type of communication.");
	}
}


/*******************  FUNCTION  *********************/
/**
 * Debut de communications asynchrones
 * @param mesh_comm MeshComm à utiliser
 * @param mesh_to_process Mesh a utiliser lors de l'échange des mailles fantomes
**/
void lbm_comm_sync_ghosts_vertical( lbm_comm_t * mesh, Mesh *mesh_to_process, lbm_comm_type_t comm_type, int target_rank, int y )
{
	//vars
	MPI_Status status;
	int x, k;

	//if target is -1, no comm
	if (target_rank == -1)
		return;

	switch (comm_type)
	{
		case COMM_SEND:
			for ( x = 1 ; x < mesh_to_process->width - 2 ; x++)
			{
			for ( k = 0 ; k < DIRECTIONS ; k++)
				{
				MPI_Send( &Mesh_get_cell(mesh_to_process, x, y)[k], 1, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD);
				}
			}

			break;
		case COMM_RECV:
			for ( x = 1 ; x < mesh_to_process->width - 2 ; x++)
			{
				for ( k = 0 ; k < DIRECTIONS ; k++)
				{
						MPI_Recv( &Mesh_get_cell(mesh_to_process, x, y)[k], 1, MPI_DOUBLE, target_rank, 0, MPI_COMM_WORLD,&status);
					}
				}


			break;
		default:
			fatal("Unknown type of communication.");
	}
}

/*******************  FUNCTION  *********************/
void lbm_comm_ghost_exchange(lbm_comm_t * mesh, Mesh *mesh_to_process )
{
	//vars
	int rank;

	//get rank
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	//Left to right phase : on reçoit à droite et on envoie depuis la gauche
	lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->right_id,mesh->width - 2);
	lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->left_id,0);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	// Right to left phase : on reçoit à gauche et on envoie depuis la droite
	lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->left_id,1);
	lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->right_id,mesh->width - 1);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	//top to bottom : on reçoit en bas et on envoie depuis le hauteur
	lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->bottom_id,mesh->height - 2);
	lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->top_id,0);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	// Right to left phase : on reçoit en haut et on envoie depuis le bas
	lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_SEND,mesh->top_id,1);
	lbm_comm_sync_ghosts_vertical(mesh,mesh_to_process,COMM_RECV,mesh->bottom_id,mesh->height - 1);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	//top left
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_LEFT],1,1);
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 1,mesh->height - 1);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	//bottom left
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_LEFT],1,mesh->height - 2);
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 1,0);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	//top right
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 2,1);
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_BOTTOM_LEFT],0,mesh->height - 1);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	//bottom left
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_LEFT],1,mesh->height - 2);
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_RIGHT],mesh->width - 1,0);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	//bottom right
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_SEND,mesh->corner_id[CORNER_BOTTOM_RIGHT],mesh->width - 2,mesh->height - 2);
	lbm_comm_sync_ghosts_diagonal(mesh,mesh_to_process,COMM_RECV,mesh->corner_id[CORNER_TOP_LEFT],0,0);

	//prevend comm mixing to avoid bugs
	MPI_Barrier(MPI_COMM_WORLD);

	// Right to left phase : on reçoit à gauche et on envoie depuis la droite
	lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_SEND,mesh->left_id,1);
	lbm_comm_sync_ghosts_horizontal(mesh,mesh_to_process,COMM_RECV,mesh->right_id,mesh->width - 1);

	//wait for IO to finish, VERY important, do not remove.
	//FLUSH_INOUT();
}

/*******************  FUNCTION  *********************/
/**
 * Rendu du mesh en effectuant une réduction a 0
 * @param mesh_comm MeshComm à utiliser
 * @param temp Mesh a utiliser pour stocker les segments
**/
void save_frame_all_domain( FILE * fp, Mesh *source_mesh, Mesh *temp )
{
	//vars
	int i = 0;
	int comm_size, rank ;
	MPI_Status status;
	//get rank and comm size
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );

	/* If whe have more than one process */
	if( 1 < comm_size )
	{
		if( rank == 0 )
		{
		unsigned long long int size = source_mesh->width  * source_mesh->height * DIRECTIONS;
			/* Rank 0 renders its local Mesh */
			save_frame(fp,source_mesh);
			/* Rank 0 receives & render other processes meshes */
			for( i = 1 ; i < comm_size ; i++ )
			{
//				MPI_Recv( temp->cells, source_mesh->width  * source_mesh->height * DIRECTIONS, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				MPI_Recv( temp->cells, size, MPI_DOUBLE, i, 0, MPI_COMM_WORLD, &status );
				save_frame(fp,temp);
			}
		} else {
			/* All other ranks send their local mesh */
			MPI_Send( source_mesh->cells, source_mesh->width * source_mesh->height * DIRECTIONS, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD );
		}
	} else {
		/* Only 0 renders its local mesh */
		save_frame(fp,source_mesh);
	}

}
