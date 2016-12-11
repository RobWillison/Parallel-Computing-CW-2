#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <mpi.h>

//The size of the grid
int size = 4;
int numberOfNodes = 2;
double precision = 0.01;

/**
* getMatrix
* width : the number of cells accross the matrix needs to be
* height : the number of cells down the matrix needs to be
*
* returns a double ** contatining a zeroed matrix width x height
**/
double **getMatrix(int width, int height)
{
  double** matrix;

  matrix = malloc(width * sizeof(double*));

  int i;
  for (i = 0; i < width; i++) {
    matrix[i] = calloc(height, sizeof(double));
  }

  return matrix;
}
/**
* printArray
* matrix : the matrix to print
*
* returns void
**/
void printArray(double **matrix)
{
  int x;
  int y;
  for (x = 0; x < size; x++) {
    for (y = 0; y < size; y++) {
      printf("%f ", matrix[x][y]);
    }
    printf("\n");
  }
}
/**
* relaxRow
* row : the intager relating to the index of the row in the matrix
*       to perform relaxation on
*
* sets cont to 1 or 0 depending on whether all cells met the precision
**/
int relaxRow(int row, double** readMatrix, double** writeMatrix, double precision)
{
  //a double to keep track of the largest diffrence
  int x;
  int cont = 0;

  //Loop through all cells in the row
  for (x = 1; x < size - 1; x++)
  {
      //calculate the average of all 4 neighbours
      double temp = readMatrix[row][x + 1];
      temp = temp + readMatrix[row][x - 1];
      temp = temp + readMatrix[row + 1][x];
      temp = temp + readMatrix[row - 1][x];
      temp = temp / 4.0f;
      //write the answer into the writeMatrix
      writeMatrix[row][x] = temp;
      //quick check if cont allready set then
      //don't bother checking precision for this cell
      if (cont) continue;
      //Work out the diffrence betoween new and old values
      double diffrence = readMatrix[row][x] - temp;
      //Check if the diffrence is bigger than the current max
      if (diffrence < 0)
      {
        diffrence = diffrence * -1.0f;
      }

      //If the precision is smaller than the diffrence set cont to true
      if (precision < diffrence) cont = 1;
  }

  return cont;
}
/**
* setupMatrix
* sets the read and write global matrixs with the boundary values
**/
void setupMatrix(double** readMatrix, double** writeMatrix)
{
  double boundingValues = 10;
  int x;
  //fill sides
  for (x = 0; x < size; x++)
  {
    readMatrix[x][0] = boundingValues;
    readMatrix[x][size - 1] = boundingValues;
    writeMatrix[x][0] = boundingValues;
    writeMatrix[x][size - 1] = boundingValues;
  }
  //fill top and bottom
  for (x = 0; x < size; x++)
  {
    readMatrix[0][x] = boundingValues;
    readMatrix[size - 1][x] = boundingValues;
    writeMatrix[0][x] = boundingValues;
    writeMatrix[size - 1][x] = boundingValues;
  }
}

/**
* getChunckSize
* Works out how many rows to give each node
* returns an array of the same length as there are nodes,
* each item is the number of rows to give that node
**/
int* getChunckSize(int matrixSize, int numberOfNodes)
{
  matrixSize = matrixSize - 2; //remove the two edges which are fixed
  int* chunkSize = (int*)malloc(sizeof(int) * numberOfNodes);
  int roundHigh = ceil((float) matrixSize / (float) numberOfNodes);
  int rowsLeft = matrixSize;

  int i;
  for (i = 0; i < numberOfNodes; i++)
  {
    if (rowsLeft - roundHigh > 0)
    {
      chunkSize[i] = roundHigh;
      rowsLeft = rowsLeft - roundHigh;
    } else if (rowsLeft > 0){
      chunkSize[i] = rowsLeft;
      rowsLeft = 0;
    } else {
      chunkSize[i] = 0;
    }
    //printf("%d\n", chunkSize[i]);
  }

  return chunkSize;
}
/**
* relaxChunk
* This section performs one relaxation iteration on the chunkSize rows from the
* offset in the matrix. it also sends the two outside rows to the other nodes
* which will require them for the next cycle, and also recives the outer most rows
* from those neighbor nodes as they will be needed in the next iteration
**/
int relaxChunk(double** readMatrix, double** writeMatrix, int chunkSize, int offset, double precision, int rank, int size)
{
  MPI_Request request, recieve;
  int cont = 0;
  //relax two outside rows, the ones required by other processes
  if(relaxRow(offset, readMatrix, writeMatrix, precision)) cont = 1;
  if(relaxRow(offset + (chunkSize - 1), readMatrix, writeMatrix, precision)) cont = 1;
  //If this isnt the first chunk send the leftmost row to the node doing the chunck to the left
  if (offset > 1)
  {
    int target = rank - 1;
    MPI_Isend(
      writeMatrix[offset],
      size,
      MPI_DOUBLE,
      target,
      offset,
      MPI_COMM_WORLD,
      &request
    );
  }
  //If this isnt the last chunk send the rightmost row to the node doing the chunck to the right
  if (offset + (chunkSize - 1) < size - 2)
  {
    int target = rank + 1;
    MPI_Isend(
      writeMatrix[offset + (chunkSize - 1)],
      size,
      MPI_DOUBLE,
      target,
      offset + (chunkSize - 1),
      MPI_COMM_WORLD,
      &request
    );
  }
  //Relax the rest of the rows
  int i;
  for (i = offset + 1; i < offset + (chunkSize - 1); i++)
  {
    if(relaxRow(i, readMatrix, writeMatrix, precision)) cont = 1;
  }
  //If this isnt the first chunk recieve the rightmost row from the node to the left
  if (offset > 1)
  {
    MPI_Status stat;
    int target = rank - 1;
    MPI_Irecv(
      writeMatrix[offset - 1],
      size,
      MPI_DOUBLE,
      target,
      offset - 1,
      MPI_COMM_WORLD,
      &recieve
    );
  }
  //If this isnt the last chunk recieve the leftmost row from the node to the right
  if (offset + (chunkSize - 1) < size - 2)
  {
    MPI_Status stat;
    int target = rank + 1;
    MPI_Irecv(
      writeMatrix[offset + (chunkSize - 1) + 1],
      size,
      MPI_DOUBLE,
      target,
      offset + (chunkSize - 1) + 1,
      MPI_COMM_WORLD,
      &recieve
    );
  }


  return cont;
}
/**
* flatternMatrixChunk
* Convert a matrix into a 1D array for sending to another node
**/
double* flatternMatrixChunk(double** matrix, int size, int offset, int chuckSize)
{
  double* array = (double*)malloc(sizeof(double) * chuckSize * size);
  int i;
  int j;
  int currentIndex = 0;
  for(i = offset; i < offset + chuckSize; i++)
  {
    for(j = 0; j < size; j++)
    {
      array[currentIndex] = matrix[i][j];
      currentIndex++;
    }
  }
  return array;
}
/**
* writeArrayIntoMatrix
* write a 1D array into the matrix, oposite of flatternMatrixChunk
*/
void writeArrayIntoMatrix(double* chunk, double** matrix, int offset, int chuckSize, int size)
{
  int i;
  int currentIndex = 0;
  for(i = offset; i < offset + chuckSize; i++)
  {
    int j;
    for(j = 0; j < size; j++)
    {
      matrix[i][j] = chunk[currentIndex];
      currentIndex++;
    }
  }
}

int main(int argc, char **argv)
{
  double starttime, endtime;

  if(argc <= 2) {
      printf("No Arguments");
      exit(1);
  }
  size = atoi(argv[1]);
  numberOfNodes = atoi(argv[2]);

  double** readMatrix = getMatrix(size, size);
  double** writeMatrix = getMatrix(size, size);

  setupMatrix(readMatrix, writeMatrix);

  printf("Running Size %d on %d\n", size, numberOfNodes);

  int rank;
  int rc = MPI_Init(NULL, NULL);
  //Record start time
  starttime = MPI_Wtime();
  //Check for error
  if (rc != MPI_SUCCESS) {
    printf ("Error\n");
    MPI_Abort(MPI_COMM_WORLD, rc);
  }
  //Get rank
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //Get the chunk sizes
  int* chunkSize = getChunckSize(size,numberOfNodes);
  //Work out what the chunk size and offset are for this node
  int offset = 1;
  int i;
  for (i = 0; i < rank; i++)
  {
    offset = offset + chunkSize[i];
  }
  //Which every process hasen't finished
  int cont = 1;
  while (cont)
  {
    //If I have a chunk size of more than 0, run one relaxation iteration
    if (chunkSize[rank] != 0){
      cont = relaxChunk(readMatrix, writeMatrix, chunkSize[rank], offset, precision, rank, size);
    } else {
      cont = 0;
    }

    int addedCont = 0;
    //Let the processor with rank 0 reduce all the cont variables by addition
    MPI_Reduce(&cont, &addedCont, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    cont = addedCont;
    //Broadcast the result of the reduce to all nodes and recive in cont
    MPI_Bcast(&cont, 1, MPI_INT, 0, MPI_COMM_WORLD);
    //Swap the matrixs round if the cont value is non zero
    if (cont) {
      double** temp = readMatrix;
      readMatrix = writeMatrix;
      writeMatrix = temp;
    }
  }
  //Once relaxation has reached the precision
  //Let the processor of rank 0 recieve all relaxed chunks and build the matrix
  if (rank == 0){
    int i;
    int offset = chunkSize[0] + 1;
    //For each node
    for(i = 1; i < numberOfNodes; i++)
    {
      MPI_Status stat;
      double* array = (double*)malloc(sizeof(double) * chunkSize[i] * size);
      //recieve the matrix chunk from that node
      MPI_Recv(
        array,
        chunkSize[i] * size,
        MPI_DOUBLE,
        i,
        0,
        MPI_COMM_WORLD,
        &stat
      );
      //Insert that chunk into this nodes writeMatrix
      writeArrayIntoMatrix(array, writeMatrix, offset, chunkSize[i], size);
      offset = offset + chunkSize[i];
    }
    printArray(writeMatrix);
  } else {
    //If this node hasen't got a rank of 0 send its chunk to the node of rank 0
    double* array = flatternMatrixChunk(writeMatrix, size, offset, chunkSize[rank]);
    MPI_Send(array, chunkSize[rank] * size, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD);
  }
  //Finish everything
  MPI_Finalize();
  //Get the run time
  endtime   = MPI_Wtime();
  printf("That took %f\n",endtime - starttime);

  return 0;
}
