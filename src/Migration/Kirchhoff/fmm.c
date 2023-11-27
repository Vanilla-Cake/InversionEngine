
#include "kirchhoff.h"
#include "main.h"

void bandw_insert(int insertz, int insertx, float insertt, int *bandwx, int *bandwz, float *bandwt, int *isroot, int *csize, int msize)
{
    if (*isroot == 1)
    {
        bandwz[0] = insertx;
        bandwx[0] = insertz;
        bandwt[0] = insertt;
        (*csize)++;
        *isroot = 0;
        // printf("infunc bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
        filterdown(bandwx, bandwz, bandwt, csize, msize);
        // printf("infunc bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
    }
    else
    {
        bandwz[*csize] = insertx;
        bandwx[*csize] = insertz;
        bandwt[*csize] = insertt;
        (*csize)++;
        filterup(bandwx, bandwz, bandwt, csize, msize);
    }
}

void bandw_remove(int *imin, int *jmin, float *tmin, int *bandwx, int *bandwz, float *bandwt, int *isroot, int *csize, int msize)
{
    if (*isroot == 1)
    {
        bandwx[0] = bandwx[*csize];
        bandwz[0] = bandwz[*csize];
        bandwt[0] = bandwt[*csize];
        bandwz[*csize] = 0;
        bandwx[*csize] = 0;
        bandwt[*csize] = 100.0;
        filterdown(bandwx, bandwz, bandwt, csize, msize);
    }
    *imin = bandwx[0];
    *jmin = bandwz[0];
    *tmin = bandwt[0];
    (*csize)--;
    *isroot = 1;
}

void filterup(int *bandwx, int *bandwz, float *bandwt, int *csize, int msize)
{
    int i, j, tempx, tempz;
    float tempt;

    j = *csize;
    i = j / 2;
    tempx = bandwx[j - 1];
    tempz = bandwz[j - 1];
    tempt = bandwt[j - 1];

    while (j > 1)
    {
        if (bandwt[i - 1] > tempt)
        {
            bandwx[j - 1] = bandwx[i - 1];
            bandwz[j - 1] = bandwz[i - 1];
            bandwt[j - 1] = bandwt[i - 1];
            j = i;
            i = j / 2;
            if (j > 1)
                continue;
        }
        break;
    }

    bandwx[j - 1] = tempx;
    bandwz[j - 1] = tempz;
    bandwt[j - 1] = tempt;
}

void filterdown(int *bandwx, int *bandwz, float *bandwt, int *csize, int msize)
{
    int i, j, tempx, tempz;
    float tempt;

    i = 1;
    j = 2 * i;
    tempx = bandwx[i-1];
    tempz = bandwz[i-1];
    tempt = bandwt[i-1];
// printf("bandwx[1]=%d, bandwz[1]=%d, bandwt[1]=%f\n", bandwx[25], bandwz[1], bandwt[1]);


    while (j <= *csize)
    {
        if (j + 1 <= *csize && bandwt[j-1] > bandwt[j])
            j = j + 1;
        if (tempt > bandwt[j-1])
        {
            bandwx[i-1] = bandwx[j-1];
            bandwz[i-1] = bandwz[j-1];
            bandwt[i-1] = bandwt[j-1];
            i = j;
            j = 2 * i;
            if (j <= *csize)
                continue;
        }
        break;
    }

    bandwx[i-1] = tempx;
    bandwz[i-1] = tempz;
    bandwt[i-1] = tempt;
   
}

float ltifun(float x, float y, float z, float h)
{
    if ((x > y) && (x < (y + sqrt(2.0) * z * h / 2.0)))
    {
        return x + sqrt((z * h) * (z * h) - (x - y) * (x - y));
    }
    else
    {
        if (x <= y)
        {
            return x + z * h;
        }
        else
        {
            return y + sqrt(2.0) * z * h;
        }
    }
}

void ltifmm(int n, int m, int sj, int si, float h, int msize, float **alive,float **v)
{
    int imin, jmin, i, j, csize, isroot, rcsize, kg;
    float bandmin, remin, tmin, isbandwt, txmin, tzmin, x, z;
    float recpt[n][m], bandwt[msize];
    int bandwx[msize], bandwz[msize];


    // Input the slowness parameter
    // FILE *file = fopen("vmodel.txt", "r");
    // for (i = 0; i < n; i++)
    // {
    //     for (j = 0; j < m; j++)
    //     {
    //         fscanf(file, "%f", &v[i][j]);
    //     }
    // }
    // fclose(file);

    // for (i = 0; i < n; i++)
    // {
    //     for (j = 0; j < m; j++)
    //     {
    //         v[i][j] = 1000;
    //     }
    // }

    // Initialize Alive, Narrow, Far point
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < m; j++)
        {
            recpt[i][j] = 100.0;
            alive[i][j] = 100.0;
        }
    }
    float tsurf = 100.0;
    for (i = 0; i < msize; i++)
    {
        bandwx[i] = 0;
        bandwz[i] = 0;
        bandwt[i] = 100.0;
    }
    csize = 0;
    isroot = 1;
    recpt[si][sj] = 0.0;
    alive[si][sj] = 0.0;

    if (si - 1 >= 0 && sj - 1 >= 0)
    {
        recpt[si - 1][sj - 1] = sqrt(2.0) * h / v[si - 1][sj - 1];
        alive[si - 1][sj - 1] = sqrt(2.0) * h / v[si - 1][sj - 1];
        bandw_insert(si - 1, sj - 1, recpt[si - 1][sj - 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }
    if (si - 1 >= 0)
    {
        recpt[si - 1][sj] = h / v[si - 1][sj];
        alive[si - 1][sj] = h / v[si - 1][sj];
        bandw_insert(si - 1, sj, recpt[si - 1][sj], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }
    if (si - 1 >= 0 && sj + 1 < m)
    {
        recpt[si - 1][sj + 1] = sqrt(2.0) * h / v[si - 1][sj + 1];
        alive[si - 1][sj + 1] = sqrt(2.0) * h / v[si - 1][sj + 1];
        bandw_insert(si - 1, sj + 1, recpt[si - 1][sj + 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }
    if (si + 1 < n && sj - 1 >= 0)
    {
        recpt[si + 1][sj - 1] = sqrt(2.0) * h / v[si + 1][sj - 1];
        alive[si + 1][sj - 1] = sqrt(2.0) * h / v[si + 1][sj - 1];
        bandw_insert(si + 1, sj - 1, recpt[si + 1][sj - 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }
    if (si + 1 < n)
    {
        recpt[si + 1][sj] = h / v[si + 1][sj];
        alive[si + 1][sj] = h / v[si + 1][sj];
        bandw_insert(si + 1, sj, recpt[si + 1][sj], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }
    if (si + 1 < n && sj + 1 < m)
    {
        recpt[si + 1][sj + 1] = sqrt(2.0) * h / v[si + 1][sj + 1];
        alive[si + 1][sj + 1] = sqrt(2.0) * h / v[si + 1][sj + 1];
        bandw_insert(si + 1, sj + 1, recpt[si + 1][sj + 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }
    if (sj - 1 >= 0)
    {
        recpt[si][sj - 1] = h / v[si][sj - 1];
        alive[si][sj - 1] = h / v[si][sj - 1];
        bandw_insert(si, sj - 1, recpt[si][sj - 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }
    if (sj + 1 < m)
    {
        recpt[si][sj + 1] = h / v[si][sj + 1];
        alive[si][sj + 1] = h / v[si][sj + 1];
        bandw_insert(si, sj + 1, recpt[si][sj + 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
    }

    while (csize > 0)
    {

        bandw_remove(&imin, &jmin, &tmin, bandwx, bandwz, bandwt, &isroot, &csize, msize);

        alive[imin][jmin] = tmin;

 
        if (imin - 1 >= 0)
        {
            if (alive[imin - 1][jmin] == 100.0)
            {
                isbandwt = recpt[imin - 1][jmin];

                if (jmin - 1 >= 0 && alive[imin][jmin - 1] != 100.0)
                {
                    if (jmin + 1 < m && alive[imin][jmin + 1] != 100.0)
                    {
                        txmin = fmin(alive[imin][jmin - 1], alive[imin][jmin + 1]);
                    }
                    else
                    {
                        txmin = alive[imin][jmin - 1];
                    }
                    recpt[imin - 1][jmin] = ltifun(alive[imin][jmin], txmin, 1.0 / v[imin - 1][jmin], h);
                }
                else
                {
                    if (jmin + 1 < m && alive[imin][jmin + 1] != 100.0)
                    {
                        recpt[imin - 1][jmin] = ltifun(alive[imin][jmin], alive[imin][jmin + 1], 1.0 / v[imin - 1][jmin], h);
                    }
                    else
                    {
                        recpt[imin - 1][jmin] = alive[imin][jmin] + h / v[imin - 1][jmin];
                    }
                }

                if (isbandwt == 100.0)
                {
                //    printf("bi00 imin=%d, jmin=%d, tmin=%f,isroot=%d,csize=%d,msize=%d\n", imin, jmin, tmin, isroot, csize, msize);
                    bandw_insert(imin - 1, jmin, recpt[imin - 1][jmin], bandwx, bandwz, bandwt, &isroot, &csize, msize);
                    // printf("bi bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
                }
                else
                {
                    if (recpt[imin - 1][jmin] < isbandwt)
                    {
                        rcsize = csize + 1;
                        for (i = 0; i < rcsize; i++)
                        {
                            if (bandwx[i] == imin - 1 && bandwz[i] == jmin)
                            {
                                bandwt[i] = recpt[imin - 1][jmin];
                                // Call filterup function
                                filterup(bandwx, bandwz, bandwt, &i, msize);
                                // printf("fu bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
                            }
                        }
                    }
                    else
                    {
                        recpt[imin - 1][jmin] = isbandwt;
                    }
                }
            }
        }
// printf("11bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
        if (imin + 1 < n)
        {
            if (alive[imin + 1][jmin] == 100.0)
            {
                isbandwt = recpt[imin + 1][jmin];

                if (jmin - 1 >= 0 && alive[imin][jmin - 1] != 100.0)
                {
                    if (jmin + 1 < m && alive[imin][jmin + 1] != 100.0)
                    {
                        txmin = fmin(alive[imin][jmin - 1], alive[imin][jmin + 1]);
                    }
                    else
                    {
                        txmin = alive[imin][jmin - 1];
                    }
                    recpt[imin + 1][jmin] = ltifun(alive[imin][jmin], txmin, 1.0 / v[imin + 1][jmin], h);
                }
                else
                {
                    if (jmin + 1 < m && alive[imin][jmin + 1] != 100.0)
                    {
                        recpt[imin + 1][jmin] = ltifun(alive[imin][jmin], alive[imin][jmin + 1], 1.0 / v[imin + 1][jmin], h);
                    }
                    else
                    {
                        recpt[imin + 1][jmin] = alive[imin][jmin] + h / v[imin + 1][jmin];
                    }
                }

                if (isbandwt == 100.0)
                {
                    bandw_insert(imin + 1, jmin, recpt[imin + 1][jmin], bandwx, bandwz, bandwt, &isroot, &csize, msize);
                }
                else
                {
                    if (recpt[imin + 1][jmin] < isbandwt)
                    {
                        rcsize = csize + 1;
                        for (i = 0; i < rcsize; i++)
                        {
                            if (bandwx[i] == imin + 1 && bandwz[i] == jmin)
                            {
                                bandwt[i] = recpt[imin + 1][jmin];
                                // Call filterup function
                                filterup(bandwx, bandwz, bandwt, &i, msize);
                            }
                        }
                    }
                    else
                    {
                        recpt[imin + 1][jmin] = isbandwt;
                    }
                }
            }
        }

// printf("22bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
 
        if (jmin + 1 < m)
        {
            if (alive[imin][jmin + 1] == 100.0)
            {
                isbandwt = recpt[imin][jmin + 1];

                if (imin - 1 >= 0 && alive[imin - 1][jmin] != 100.0)
                {
                    if (imin + 1 < n && alive[imin + 1][jmin] != 100.0)
                    {
                        tzmin = fmin(alive[imin - 1][jmin], alive[imin + 1][jmin]);
                    }
                    else
                    {
                        tzmin = alive[imin - 1][jmin];
                    }
                    recpt[imin][jmin + 1] = ltifun(alive[imin][jmin], tzmin, 1.0 / v[imin][jmin + 1], h);
                }
                else
                {
                    if (imin + 1 < n && alive[imin + 1][jmin] != 100.0)
                    {
                        recpt[imin][jmin + 1] = ltifun(alive[imin][jmin], alive[imin + 1][jmin], 1.0 / v[imin][jmin + 1], h);
                    }
                    else
                    {
                        recpt[imin][jmin + 1] = alive[imin][jmin] + h / v[imin][jmin + 1];
                    }
                }

                if (isbandwt == 100.0)
                {
                    bandw_insert(imin, jmin + 1, recpt[imin][jmin + 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
                }
                else
                {
                    if (recpt[imin][jmin + 1] < isbandwt)
                    {
                        rcsize = csize + 1;
                        for (i = 0; i < rcsize; i++)
                        {
                            if (bandwx[i] == imin && bandwz[i] == jmin + 1)
                            {
                                bandwt[i] = recpt[imin][jmin + 1];
                                // Call filterup function
                                filterup(bandwx, bandwz, bandwt, &i, msize);
                            }
                        }
                    }
                    else
                    {
                        recpt[imin][jmin + 1] = isbandwt;
                    }
                }
            }
        }
// printf("33bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
 
        if (jmin - 1 >= 0)
        {
            if (alive[imin][jmin - 1] == 100.0)
            {
                isbandwt = recpt[imin][jmin - 1];

                if (imin - 1 >= 0 && alive[imin - 1][jmin] != 100.0)
                {
                    if (imin + 1 < n && alive[imin + 1][jmin] != 100.0)
                    {
                        tzmin = fmin(alive[imin - 1][jmin], alive[imin + 1][jmin]);
                    }
                    else
                    {
                        tzmin = alive[imin - 1][jmin];
                    }
                    recpt[imin][jmin - 1] = ltifun(alive[imin][jmin], tzmin, 1.0 / v[imin][jmin - 1], h);
                }
                else
                {
                    if (imin + 1 < n && alive[imin + 1][jmin] != 100.0)
                    {
                        recpt[imin][jmin - 1] = ltifun(alive[imin][jmin], alive[imin + 1][jmin], 1.0 / v[imin][jmin - 1], h);
                    }
                    else
                    {
                        recpt[imin][jmin - 1] = alive[imin][jmin] + h / v[imin][jmin - 1];
                    }
                }

                if (isbandwt == 100.0)
                {
                    bandw_insert(imin, jmin - 1, recpt[imin][jmin - 1], bandwx, bandwz, bandwt, &isroot, &csize, msize);
                }
                else
                {
                    if (recpt[imin][jmin - 1] < isbandwt)
                    {
                        rcsize = csize + 1;
                        for (i = 0; i < rcsize; i++)
                        {
                            if (bandwx[i] == imin && bandwz[i] == jmin - 1)
                            {
                                bandwt[i] = recpt[imin][jmin - 1];
                                // Call filterup function
                                filterup(bandwx, bandwz, bandwt, &i, msize);
                            }
                        }
                    }
                    else
                    {
                        recpt[imin][jmin - 1] = isbandwt;
                    }
                }
            }
        }
        // printf("44bandwx[0]=%d, bandwz[0]=%d, bandwt[0]=%f\n", bandwx[0], bandwz[0], bandwt[0]);
 
    }
}
