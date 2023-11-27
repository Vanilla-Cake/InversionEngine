#include "kirchhoff.h"
#include "main.h"
#include <time.h>
#include "par.h"

int xargc;
char **xargv;

typedef struct headInformation
{
      int is;
      float isp;
      float firec;
      int ntr;
      int tataltr;
} headI;

int main(int argc, char *argv[])
{

      int i, j, k;
      int is, ir, xs, zs, xr, zr;
      int ix, iz, it;
      int itp1;

      int sp, firp, ntr, totaltr;

      int nx, nz; // x方向和z方向的网格数
      float dxz;  // 网格间距
      int msize;  // 一维数组的长度
      int nshot, nt;
      float dtr, dt;
      float aa6;

      int ns1, ns2; // 添加参数定义
      float t;      // 添加参数定义
      float d;

      int n1r, n2r; // 添加参数定义

      float tsamp; // 添加参数定义
      int ntstep;  // 添加参数定义

      float **v, **tts, **ttr, **agrid;
      float **f;
      float **ff;

      char *headInfo, *vfile, *seifile;
      FILE *headInfofp, *vfp, *seifp; // 道头文件，速度文件，地震记录文件

      char *imgfile;
      FILE *imgfp;

      // Start---Input the parameter and the model

      int nxp1 = nx;
      int nzp1 = nz;

      initargs(argc, argv);

      if (!getparint("nx", &nx))
            nx = 100;
      if (!getparint("nz", &nz))
            nz = 100;
      if (!getparfloat("dxz", &dxz))
            dxz = 10.0;
      if (!getparint("msize", &msize))
            msize = 10 * nx;

      if (!getparint("ntr", &ntr))
            ntr = 100;
      if (!getparint("nshot", &nshot))
            nshot = 1;
      if (!getparint("nt", &nt))
            nt = 750;
      if (!getparfloat("dtr", &dtr))
            dtr = dxz;
      if (!getparfloat("dt", &dt))
            dt = 0.004;

      if (!getparstring("headInfo", &headInfo))
            headInfo = "headInfo";
      if (!getparstring("vfile", &vfile))
            vfile = "v.bin";
      if (!getparstring("seifile", &seifile))
            seifile = "sei.bin";
      if (!getparstring("imgfile", &imgfile))
            imgfile = "img.bin";

      ntstep = nt;
      tsamp = dt; // 历史遗留问题

      v = alloc2float(nz, nx);
      tts = alloc2float(nz, nx);
      ttr = alloc2float(nz, nx);
      agrid = alloc2float(nz, nx);

      // read head info
      headInfofp = fopen(headInfo, "r");
      headI *headInformation;
      headInformation = (headI *)malloc(sizeof(headI) * nshot);
      for (i = 0; i < nshot; i++)
      {
            fscanf(headInfofp, "%d,%f,%f,%d,%d\n", &headInformation[i].is, &headInformation[i].isp,
                   &headInformation[i].firec, &headInformation[i].ntr, &headInformation[i].tataltr);
      }
      fclose(headInfofp);

      // read velocity model
      vfp = fopen(vfile, "rb");
      for (i = 0; i < nx; i++)
      {
            for (j = 0; j < nz; j++)
            {
                  fread(&v[i][j], sizeof(float), 1, vfp);
            }
      }
      fclose(vfp);

      // open other file
      seifp = fopen(seifile, "rb");
      imgfp = fopen(imgfile, "wb");

      // initialize tts ,ttr and agrid
      for (i = 0; i < nx; i++)
      {
            for (j = 0; j < nz; j++)
            {
                  tts[i][j] = 0.0;
                  ttr[i][j] = 0.0;
                  agrid[i][j] = 0.0;
            }
      }

      // Start---summation loops

      for (is = 0; is < nshot; is++)
      {

            printf("\n%d------Beginning of computation......\n", is);

            sp = headInformation[is].isp / dxz;
            firp = headInformation[is].firec / dxz;
            ntr = headInformation[is].ntr;
            totaltr = headInformation[is].tataltr;

            // Loop over source positions
            xs = sp;
            zs = 0;
            if(xs>nx-1)
            {
                  continue;
            }
            ltifmm(nx, nz, zs, xs, dxz, msize, tts, v);

            f = alloc2float(nt, ntr);
            ff = alloc2float(nt, ntr);

            fseek(seifp, sizeof(float) * totaltr * nt, SEEK_SET);
            for (i = 0; i < ntr; i++)
            {
                  fread(ff[i], sizeof(float), nt, seifp);
            }
            fseek(seifp, 0, SEEK_SET);

            for (ir = 0; ir < ntr; ir++)
            {
                  f[ir][0] = (ff[ir][1] - ff[ir][0]) / tsamp;
                  for (it = 1; it < nt; it++)
                  {
                        f[ir][it] = (ff[ir][it] - ff[ir][it - 1]) / tsamp;
                  }
            }


            // Loop over receiver positions

            for (ir = 0; ir < ntr; ir++)
            {
                  xr = (firp*dxz + ir * dtr) / dxz;
                  zr = 0;
                  if(xr>nx-1)
                  {
                        continue;
                  }

                  // Call ltifmm function
                  ltifmm(nx, nz, zr, xr, dxz, msize, ttr, v);

                  // printf("         summation\n");

                  for (ix = 0; ix < nx; ix++)
                  {
                        for (iz = 0; iz < nz; iz++)
                        {
                              t = tts[ix][iz] + ttr[ix][iz];
                              it = (int)(t / tsamp);
                              if (it > ntstep)
                                    continue;
                              itp1 = it + 1;
                              if (itp1 > ntstep)
                                    itp1 = it;
                              aa6 = f[ir][it] + (f[ir][itp1] - f[ir][it]) * (t - it * tsamp) / tsamp;
                              agrid[ix][iz] += aa6;
                        }
                  }
            }
            free2float(f);
            free2float(ff);
      }

      // Write agrid to image file
      for (i = 0; i < nx; i++)
      {
            for (j = 0; j < nz; j++)
            {
                  fwrite(&agrid[i][j], sizeof(float), 1, imgfp);
            }
      }

      fclose(seifp);
      fclose(imgfp);

      return 0;
}
