#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define TYPE int
#define MPI_TYPE MPI_INT
#define ROOT 0

/*
Este es el esquema:

Mesh 2D MPI ranks -> p x q
--------------------------

b-> algorithmic block size
n-> matrix dimension (n x n)

Buffers-> 2 buffers n x b each

for (k=0; k<n; k+=b) {
  
  m = min(n-k+1,b);


  k) Iteration
  k.1)  All processes within row k/b broadcast to processes in the same column

      Example k=0
      Process (0,0) broadcast to (1,0)      | Process (0,1) broadcast to (1,1)    | Process (0,2) broadcast to (1,2)   | ...
                                 (2,0)   	|                            (2,1)    |                            (2,2)   | ...
                                 ...  		|                            ...      |                            ...     | ...
                                 (p-1,0)    |                            (p-1,1)  |                            (p-1,2) | ...

      Message size: (m/q)*b

  k.2)  All processes within column k/b broadcast to processes in the same row

      Message size: (m/p)*b

}
*/

//adaptar el tamaño de la matriz con respecto al tamaño de los procesos 
void resize(size_t matrix_rows, size_t matrix_cols, size_t proc_rows, size_t proc_cols, size_t *new_matrix_rows, size_t *new_matrix_cols) {
        if(matrix_rows % proc_rows || matrix_rows < proc_rows) *new_matrix_rows = ((matrix_rows/proc_rows)+1)*proc_rows;
            else *new_matrix_rows = matrix_rows;
        if(matrix_cols % proc_cols || matrix_cols < proc_cols) *new_matrix_cols = ((matrix_cols/proc_cols)+1)*proc_cols;
            else *new_matrix_cols = matrix_cols;
}


int main(int argc, char *argv[])
{
    int rank, wsize, i;
    int reps = 1;
    int reps2;
    size_t m_rows, m_cols, proc_rows, proc_cols;
    

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    char processor_name[MPI_MAX_PROCESSOR_NAME];
    int name_len;
    MPI_Get_processor_name(processor_name, &name_len);

    // Print off a hello world message so we ensure the correct MPI mapping (one per node)
    printf("#Hello world from processor %s, rank %d out of %d processors\n",
            processor_name, rank, wsize);

    int total_wsize = wsize; // Do not modify it
    /* Warm-up zone */
    MPI_Bcast(&reps,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Allreduce(&reps,&reps2,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);
    //Reps2 is never used but we want to make the allreduce movement

    //Matriz de tamaño 1048576, repartidos de formas iguales entre filas y columnas
    m_rows = (argc > 1) ? atoi(argv[1]) : 131072; //num filas 
    m_cols = (argc > 2) ? atoi(argv[2]) : 131072; //num columnas 
    //pasar por parametros el nº de procesos
    proc_rows = (argc > 3) ? atoi(argv[3]) : 16; //pillar por parametros
    proc_cols = (argc > 4) ? atoi(argv[4]) : 8; //pillar por parametros
    
    size_t new_m_rows, new_m_cols;
    resize(m_rows, m_cols, proc_rows, proc_cols, &new_m_rows, &new_m_cols);

    size_t bp = new_m_rows / proc_rows; // tamaño bloque rows
    size_t bq = new_m_cols / proc_cols; // tamaño bloque cols
    size_t num_bloques = bp * bq; // numero de procesos p x q
    size_t bsize = proc_rows * proc_cols; //num elementos bloques
    //printf("bp %ld\t bq %ld\t num_bloques %ld\t bsize %ld\t prow_rows %ld\t prow_cols %ld\t new_w_rows %ld\t new_w_cols %ld\n",bp, bq, num_bloques, bsize, proc_rows, proc_cols, new_m_rows, new_m_cols);
    int max_size = (new_m_rows > new_m_cols) ? new_m_rows: new_m_cols;
    int *info = malloc(max_size*sizeof(int));
    
    MPI_Group * row_grs = malloc(proc_rows*sizeof(MPI_Group));
    MPI_Group * col_grs = malloc(proc_cols*sizeof(MPI_Group));

    MPI_Comm * row_comms = malloc(proc_rows*sizeof(MPI_Comm));
    MPI_Comm * col_comms = malloc(proc_cols*sizeof(MPI_Comm));
     

    //local_row = ((int) (rank / n)) % b;
    //local_col =  (rank % n) % b;
    //cada 0 de cada subloque. Ejemplo de 36 procesos, n = 6. subloques de 9 procs b = 3
    // 0    1   2   |   3   4   5           [0 1 2]  Proceso 0 tiene 2 grupos fila y col
    // 6    7   8   |   9   10  11          --
    // 12   13  14  |   15  16  17          0
    //----------------------------          6
    // 18   19  20  |   21  22  23          12
    // 24   25  26  |   27  28  29          --
    // 30   31  32  |   33  34  35

    //p*q procesos de filas p x columnas q
    MPI_Group world_group;
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);
    
    
    // Grupos filas
    int * row_ranks = malloc(wsize * sizeof(int));
    for(i = 0; i < wsize; i++){
        row_ranks[i] = i;
    }
    for(i = 0; i < proc_rows; i++){
        MPI_Group_incl(world_group, proc_rows, &row_ranks[i*proc_rows], &row_grs[i]);
        MPI_Comm_create_group(MPI_COMM_WORLD, row_grs[i], i, &row_comms[i]);
    }
    

    // Grupos columnas
    int * col_ranks = malloc(wsize * sizeof(int));
    for(i = 0; i < wsize; i++){
        col_ranks[i] = row_ranks[(i * proc_cols) % wsize + (int) (i / proc_rows)];
    } 
    for(i = 0; i < proc_cols; i++){
        MPI_Group_incl(world_group, proc_cols, &col_ranks[i*proc_cols], &col_grs[i]);
        MPI_Comm_create_group(MPI_COMM_WORLD, col_grs[i], i, &col_comms[i]);
    }


    // Envios
    int iter = (proc_rows > proc_cols) ? proc_rows: proc_cols;
    for(i = 0; i < iter; i++) { 
        //filas
        if (i < proc_rows) {
            if( (int)(rank / proc_rows) == i){
                MPI_Bcast(info, bp, MPI_INT, 0, row_comms[i]);
                //printf("%d MPI_Bcast info %ld, MPI_INT, 0, row_comms[i], %d\n", rank, bp, i);
            }
        }
        
        //barrier
        MPI_Barrier(MPI_COMM_WORLD);
        //columnas
        if (i < proc_cols) {
            if( (int)(rank % proc_cols) == i){
                MPI_Bcast(info, bq, MPI_INT, 0, col_comms[i]);
                //printf("%d MPI_Bcast info %ld, MPI_INT, 0, col_comms[i], %d\n", rank, bq, i);
            }
        }
        //barrier
        MPI_Barrier(MPI_COMM_WORLD);
    }
    

    free( info );
    free( row_ranks );
    free( col_ranks );   
    MPI_Group_free(&world_group);
    
    MPI_Finalize();
    
    return 0;
}