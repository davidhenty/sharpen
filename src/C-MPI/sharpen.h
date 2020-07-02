void pgmsize(char *filename, int *nx, int *ny);
void pgmread(char *filename, void *vp, int nxmax, int nymax, int *nx, int *ny);
void pgmwrite(char *filename, void *vx, int nx, int ny);

void dosharpen(char *filename, int nx, int ny, MPI_Comm comm);
double filter(int d, int i, int j);

int **int2Dmalloc(int nx, int ny);
double **double2Dmalloc(int nx, int ny);
