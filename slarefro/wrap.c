#include "wrap.h"

void slaf_refro ( double zobs, double hm, double tdk, double pmb,
		  double rh, double wl, double phi, double tlr,
		  double eps, double *ref )
{
	sla_refro_(&zobs, &hm, &tdk, &pmb, &rh, &wl, &phi, &tlr, &eps, ref);
}
