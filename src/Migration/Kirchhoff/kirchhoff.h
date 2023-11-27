
void filterup(int *bandwx, int *bandwz, float *bandwt, int *csize, int msize);
void filterdown(int *bandwx, int *bandwz, float *bandwt, int *csize, int msize);
void bandw_insert(int insertz, int insertx, float insertt, int *bandwx, int *bandwz, float *bandwt, int *isroot, int *csize, int msize);
void bandw_remove(int *imin, int *jmin, float *tmin, int *bandwx, int *bandwz, float *bandwt, int *isroot, int *csize, int msize);
float ltifun(float x, float y, float z, float h);

void ltifmm(int n, int m, int sj, int si, float h, int msize, float **alive,float **v);



