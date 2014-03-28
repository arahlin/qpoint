#include "qpoint.h"
#include <getdata.h>
#include <stdio.h>
#include <stdlib.h>

#define SSAMP 250000
#define NSAMP 8000000

int main(int argc, char *argv[]) {
  
  printf("init\n");
  
  int i;
  // int i,j;
  double *az;
  double *el;
  double *pitch;
  double *roll;
  double *lon;
  double *lat;
  double *ctime;
  double *ra;
  double *dec;
  double *sin2psi;
  double *cos2psi;
  // double cra, cdec, cpsi;
  quat_t *q_bore;
  double **pmap;
  double delta_az = -3.314806;
  double delta_el = -8.010956;
  double delta_psi = 22.5;
  double deltas_az[9] = {-8,-6,-4,-2,0,2,4,6,8};
  double deltas_el[9] = {-8,-6,-4,-2,0,2,4,6,8};
  double deltas_psi[9] = {22.5,-22.5,22.5,-22.5,22.5,-22.5,22.5,-22.5,22.5};
  int nside = 256;
  long npix = 12*nside*nside;
  // double framerate;
  // qp_set_offset(-3.314806,-8.010956,22.5);
  
  int mode = 0;
  if (argc>1) mode = atoi(argv[1]);
  
  if (mode == 0) {
    az = malloc(sizeof(double)*NSAMP);
    el = malloc(sizeof(double)*NSAMP);
    pitch = malloc(sizeof(double)*NSAMP);
    roll = malloc(sizeof(double)*NSAMP);
    lon = malloc(sizeof(double)*NSAMP);
    lat = malloc(sizeof(double)*NSAMP);
  } else if (mode == 1) {
    ra = malloc(sizeof(double)*NSAMP);
    dec = malloc(sizeof(double)*NSAMP);
    sin2psi = malloc(sizeof(double)*NSAMP);
    cos2psi = malloc(sizeof(double)*NSAMP);
  } else if (mode == 2) {
    pmap = malloc(6*sizeof(double*));
    for (i=0; i<6; i++)
      pmap[i] = calloc(npix, sizeof(double));
  }
  ctime= malloc(sizeof(double)*NSAMP);
  q_bore = malloc(sizeof(quat_t)*NSAMP);
  
  printf("read data\n");
  DIRFILE *D = gd_open("/Users/sasha/code/scan_sim/scan_sim_data_uldc",GD_RDONLY);
  if (mode == 0) {
    gd_getdata(D, "AZ", 0, SSAMP, 0, NSAMP, GD_DOUBLE, (void *)az);
    gd_getdata(D, "EL", 0, SSAMP, 0, NSAMP, GD_DOUBLE, (void *)el);
    gd_getdata(D, "PITCH", 0, SSAMP, 0, NSAMP, GD_DOUBLE, (void *)pitch);
    // for (i = 0; i < NSAMP; i++) pitch[i] += 0.1;
    gd_getdata(D, "ROLL", 0, SSAMP, 0, NSAMP, GD_DOUBLE, (void *)roll);
    // for (i = 0; i < NSAMP; i++) roll[i] += 0.2;
    gd_getdata(D, "LON",     0, SSAMP, 0, NSAMP, GD_DOUBLE, (void *)lon);
    gd_getdata(D, "LAT",     0, SSAMP, 0, NSAMP, GD_DOUBLE, (void *)lat);
  }
  gd_getdata(D, "t",           0, SSAMP, 0, NSAMP, GD_DOUBLE, (void *)ctime);
  // gd_get_constant(D, "FRAMERATE", GD_DOUBLE, &framerate);
  // for (i=0; i<NSAMP/framerate; i++) {
  //   for (j=0; j<framerate; j++)
  //     ctime[(int)(i*framerate+j)] += 1.0/framerate;
  // }
  
  /*
  gd_getdata(D, "RA_T00_TESTA",  0, SSAMP, 0, 1, GD_DOUBLE, (void *)&cra);
  if (cra>180) cra-=360;
  gd_getdata(D, "DEC_T00_TESTA", 0, SSAMP, 0, 1, GD_DOUBLE, (void *)&cdec);
  if (cdec>180) cdec-=360;
  gd_getdata(D, "PSI_T00_TESTA", 0, SSAMP, 0, 1, GD_DOUBLE, (void *)&cpsi);
  if (cpsi>180) cpsi-=360;
  */
  
  // gd_getdata(D, "RA_T00_X1T1R1C8A",  0, SSAMP, 0, 1, GD_DOUBLE, (void *)&cra);
  // if (cra>180) cra-=360;
  // gd_getdata(D, "DEC_T00_X1T1R1C8A", 0, SSAMP, 0, 1, GD_DOUBLE, (void *)&cdec);
  // if (cdec>180) cdec-=360;
  // gd_getdata(D, "PSI_T00_X1T1R1C8A", 0, SSAMP, 0, 1, GD_DOUBLE, (void *)&cpsi);
  // if (cpsi>180) cpsi-=360;
  
  /*
  qp_azel2quat(az[0],el[0],pitch[0],roll[0],lon[0],lat[0],ctime[0],q,0,0);
  qp_azel2radec(delta_az,delta_el,delta_psi,
		az,el,pitch,roll,lon,lat,ctime,
		ra,dec,psi,NSAMP,0,0,0);
  */
  
  qp_memory_t *mem = qp_init_memory();
  qp_set_opt_accuracy(mem, 1);
  qp_set_opt_fast_math(mem, 1);
  
  if (mode == 0) {
    printf("azel2bore\n");
    qp_azel2bore(mem, az,el,pitch,roll,lon,lat,ctime,q_bore,NSAMP);
    
    printf("save\n");
    FILE *f = fopen("qbore.dat","w");
    for (i = 0; i < NSAMP; i++) {
      fprintf(f, "%.12f\t%.12f\t%.12f\t%.12f\n", q_bore[i][0], q_bore[i][1],
	      q_bore[i][2], q_bore[i][3]);
    }
    fclose(f);
  } else if (mode == 1) {
    printf("load\n");
    FILE *f = fopen("qbore.dat","r");
    for (i = 0; i < NSAMP; i++) {
      fscanf(f, "%lf\t%lf\t%lf\t%lf\n", &q_bore[i][0], &q_bore[i][1],
	     &q_bore[i][2], &q_bore[i][3]);
    }
    fclose(f);
    
    printf("bore2rasindec\n");
    qp_bore2rasindec(mem, delta_az,delta_el,delta_psi,
		     ctime,q_bore,ra,dec,sin2psi,cos2psi,NSAMP);

    printf("save\n");
    FILE *of = fopen("test.dat","w");
    for (i = 0; i < NSAMP; i++) {
      fprintf(of, "%f\t%f\t%f\t%f\n", ra[i], dec[i], sin2psi[i], cos2psi[i]);
    }
    fclose(of);
  } else if (mode == 2) {
    printf("load\n");
    FILE *f = fopen("qbore.dat","r");
    for (i = 0; i < NSAMP; i++) {
      fscanf(f, "%lf\t%lf\t%lf\t%lf\n", &q_bore[i][0], &q_bore[i][1],
	     &q_bore[i][2], &q_bore[i][3]);
    }
    fclose(f);
    
    printf("bore2map\n");
    qp_bore2map_omp(mem, deltas_az, deltas_el, deltas_psi, 9,
		    ctime, q_bore, NSAMP, pmap, nside);
    
    printf("save\n");
    FILE *of = fopen("map_omp.dat","w");
    for (i = 0; i < npix; i++) {
      fprintf(of, "%d\t%e\t%e\t%e\t%e\t%e\n", (int)pmap[0][i], pmap[1][i], pmap[2][i],
	      pmap[3][i], pmap[4][i], pmap[5][i]);
    }
    fclose(of);
  }
  /*
  printf("quat: %.6f %.6f %.6f %.6f\n",q[0],q[1],q[2],q[3]);
  printf("FS: ra %.6f dec %.6f psi %.6f\n",cra,cdec,cpsi);
  printf("QP: ra %.6f dec %.6f psi %.6f\n",ra[0],dec[0],psi[0]);
  printf("DI: ra %.6f dec %.6f psi %.6f\n",
	 (ra[0]-cra)*3600,
	 (dec[0]-cdec)*3600,
	 (psi[0]-cpsi)*3600);
  */
  // printf("%f\t%f\t%f\t%f\t%f\n",az[0],el[0],lon[0],lat[0],ctime[0]);
  
  if (mode == 0) {
    free(az);
    free(el);
    free(pitch);
    free(roll);
    free(lon);
    free(lat);
  } else if (mode == 1) {
    free(ra);
    free(dec);
    free(sin2psi);
    free(cos2psi);
  } else if (mode == 2) {
    for (i=0; i<6; i++)
      free(pmap[i]);
    free(pmap);
  }
  free(ctime);
  free(q_bore);
  
  qp_free_memory(mem);
  
  return 0;
};
