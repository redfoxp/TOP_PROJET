/********************  HEADERS  *********************/
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <stdbool.h>
#include <math.h>
#include <stdint.h>
#include "lbm_config.h"
#include "lbm_struct.h"
#include "lbm_phys.h"
#include "lbm_init.h"
#include "lbm_comm.h"

/*******************  FUNCTION  *********************/
/**
 * Ecrit l'en-tête du fichier de sortie. Cet en-tête sert essentiellement à fournir les informations
 * de taille du maillage pour les chargements.
 * @param fp Descripteur de fichier à utiliser pour l'écriture.
**/
void write_file_header(FILE * fp,lbm_comm_t * mesh_comm)
{
	//setup header values
	lbm_file_header_t header;
	header.magick      = RESULT_MAGICK;
	header.mesh_height = MESH_HEIGHT;
	header.mesh_width  = MESH_WIDTH;
	header.lines       = mesh_comm->nb_y;

	//write file
	fwrite(&header,sizeof(header),1,fp);
}

/*******************  FUNCTION  *********************/
FILE * open_output_file(lbm_comm_t * mesh_comm)
{
	//vars
	FILE * fp;

	//check if empty filename => so noout
	if (RESULT_FILENAME == NULL)
		return NULL;

	//open result file
	fp = fopen(RESULT_FILENAME,"w");

	//errors
	if (fp == NULL)
	{
		perror(RESULT_FILENAME);
		abort();
	}

	//write header
	write_file_header(fp,mesh_comm);

	return fp;
}

void close_file(FILE* fp){
	//wait all before closing
	MPI_Barrier(MPI_COMM_WORLD);
	//close file
	fclose(fp);
}

/*******************  FUNCTION  *********************/
/**
 * Sauvegarde le résultat d'une étape de calcul. Cette fonction peu être appelée plusieurs fois
 * lors d'une sauvegarde en MPI sur plusieurs processus pour sauvegarder les un après-les autres
 * chacun des domaines.
 * Ne sont écrit que les vitesses et densités macroscopiques sous forme de flotant simple.
 * @param fp Descripteur de fichier à utiliser pour l'écriture.
 * @param mesh Domaine à sauvegarder.
**/
void save_frame(FILE * fp,const Mesh * mesh)
{
	//write buffer to write float instead of double
	lbm_file_entry_t buffer[WRITE_BUFFER_ENTRIES];
	int i,j,cnt;
	double density;
	Vector v;
	double norm;

	//loop on all values
	cnt = 0;
	for ( i = 1 ; i < mesh->width - 1 ; i++)
	{
		for ( j = 1 ; j < mesh->height - 1 ; j++)
		{
			//compute macrospic values
			density = get_cell_density(Mesh_get_cell(mesh, i, j));
			get_cell_velocity(v,Mesh_get_cell(mesh, i, j),density);
			norm = sqrt(get_vect_norme_2(v,v));

			//fill buffer
			buffer[cnt].density = density;
			buffer[cnt].v = norm;
			cnt++;

			//errors
			assert(cnt <= WRITE_BUFFER_ENTRIES);
			
			//flush buffer if full
			if (cnt == WRITE_BUFFER_ENTRIES)
			{
				fwrite(buffer,sizeof(lbm_file_entry_t),cnt,fp);
				cnt = 0;
			}
		}
	}

	//final flush
	if (cnt != 0)
		fwrite(buffer,sizeof(lbm_file_entry_t),cnt,fp);
}

/*******************  FUNCTION  *********************/
int main(int argc, char * argv[])
{
	//vars
	Mesh mesh;
	Mesh temp;
	Mesh temp_render;
	lbm_mesh_type_t mesh_type;
	lbm_comm_t mesh_comm;
	int i, rank, comm_size;
	FILE * fp = NULL;
	const char * config_filename = NULL;

	//init MPI and get current rank and commuincator size.
	MPI_Init( &argc, &argv );
	MPI_Comm_rank( MPI_COMM_WORLD, &rank );
	MPI_Comm_size( MPI_COMM_WORLD, &comm_size );

	//get config filename
	if (argc >= 2)
		config_filename = argv[1];
	else
		config_filename = "config.txt";

	//load config file and display it on master node
	load_config(config_filename);
	if (rank == RANK_MASTER)
		print_config();

	//init structures, allocate memory...
	lbm_comm_init( &mesh_comm, rank, comm_size, MESH_WIDTH, MESH_HEIGHT);
	Mesh_init( &mesh, lbm_comm_width( &mesh_comm ), lbm_comm_height( &mesh_comm ) );
	Mesh_init( &temp, lbm_comm_width( &mesh_comm ), lbm_comm_height( &mesh_comm ) );
	Mesh_init( &temp_render, lbm_comm_width( &mesh_comm ), lbm_comm_height( &mesh_comm ) );
	lbm_mesh_type_t_init( &mesh_type, lbm_comm_width( &mesh_comm ), lbm_comm_height( &mesh_comm ));

	//master open the output file
	if( rank == RANK_MASTER )
		fp = open_output_file(&mesh_comm);

	//setup initial conditions on mesh
	setup_init_state( &mesh, &mesh_type, &mesh_comm);
	setup_init_state( &temp, &mesh_type, &mesh_comm);

	//write initial condition in output file
	if (lbm_gbl_config.output_filename != NULL)
		save_frame_all_domain(fp, &mesh, &temp_render );

	//barrier to wait all before start
	MPI_Barrier(MPI_COMM_WORLD);

	//time steps
	for ( i = 1 ; i < ITERATIONS ; i++ )
	{
		//print progress
		if( rank == RANK_MASTER )
			printf("Progress [%5d / %5d]\n",i,ITERATIONS);

		//compute special actions (border, obstacle...)
		special_cells( &mesh, &mesh_type, &mesh_comm);

		//need to wait all before doing next step
		MPI_Barrier(MPI_COMM_WORLD);

		//compute collision term
		collision( &temp, &mesh);

		//need to wait all before doing next step
		MPI_Barrier(MPI_COMM_WORLD);

		//propagate values from node to neighboors
		lbm_comm_ghost_exchange( &mesh_comm, &temp );
		propagation( &mesh, &temp);

		//need to wait all before doing next step
		MPI_Barrier(MPI_COMM_WORLD);

		//save step
		if ( i % WRITE_STEP_INTERVAL == 0 && lbm_gbl_config.output_filename != NULL )
			save_frame_all_domain(fp, &mesh, &temp_render );
	}

	if( rank == RANK_MASTER && fp != NULL)
	{
		close_file(fp);
	}

	//free memory
	lbm_comm_release( &mesh_comm );
	Mesh_release( &mesh );
	Mesh_release( &temp );
	Mesh_release( &temp_render );
	lbm_mesh_type_t_release( &mesh_type );

	//close MPI
	MPI_Finalize();

	return EXIT_SUCCESS;
}
