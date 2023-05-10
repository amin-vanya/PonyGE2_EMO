/* hv_wfg.c  wfg algorithm for calculating hypervolume indicator (see \cite{While12})
 *
 * Lyndon While, Lucas Bradstreet, Luigi Barone 
 * Email lyndon.while@uwa.edu.au for any queries regarding usage/bugs/improvements. 
 *
 * This code includes a high performance implementation of the WFG algorithm, 
 * used to calculate the hypervolume indicator for a set of non-dominated points. 
 * 
 * COPYRIGHT: 
 * This program is free software (software libre); you can redistribute
 *  it and/or modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2 of the
 *  License, or (at your option) any later version.
 * 
 *  This program is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  General Public License for more details.
 * 
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, you can obtain a copy of the GNU
 *  General Public License at:
 *                  http://www.gnu.org/copyleft/gpl.html
 *  or by writing to:
 *            Free Software Foundation, Inc., 59 Temple Place,
 *                  Suite 330, Boston, MA 02111-1307 USA
 * 
 * ----------------------------------------------------------------------
 * Updates
 * Original source code taken in april 2, 2016:
 * http://www.wfg.csse.uwa.edu.au/hypervolume/#code
 * version: metric WFG implementation (22 November 2015)
 * Raquel Hernandez Gomez
 *
 */

#include <stdio.h>
//RHG #include <stdbool.h>
#include <math.h>
#include <float.h>  // RHG
//RHG #include <sys/time.h>
//RHG #include <sys/resource.h>
#include "hv_wfg.h" // RHG
#include "sort.h"   // RHG

#define BEATS(x,y)   (x > y) 
#define WORSE(x,y)   (BEATS(y,x) ? (x) : (y)) 

// RHG int n;     // the number of objectives 
// RHG POINT ref; // the reference point 

// RHG FRONT *fs;    // memory management stuff 
// RHG int fr = 0;   // current depth 
// RHG int maxm = 0; // identify the biggest fronts in the file 
// RHG int maxn = 0; 
// RHG int safe;     // the number of points that don't need sorting 

// RHG double totaltime;

// RHG double hv(FRONT);


double hv_wfg(FRONT ps, FRONT *fs, int *safe, int *fr, int n);

int greater(const void *v1, const void *v2, int start, int n)
// this sorts points worsening in the last objective
{
  int i;

  POINT p = *(POINT*)v1;
  POINT q = *(POINT*)v2;
  for (i = n - 1; i >= start; i--)
    if BEATS(p.objectives[i],q.objectives[i]) return -1;
    else
    if BEATS(q.objectives[i],p.objectives[i]) return  1;
  return 0;
}


int greaterabbrev(const void *v1, const void *v2, int start, int n)
// this sorts points worsening in the penultimate objective
{
  int i;

  POINT p = *(POINT*)v1;
  POINT q = *(POINT*)v2;
  for (i = n - 2; i >= start; i--)
    if BEATS(p.objectives[i],q.objectives[i]) return -1;
    else
    if BEATS(q.objectives[i],p.objectives[i]) return  1;
  return 0;
}


int dominates2way(POINT p, POINT q, int k)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise 
// k is the highest index inspected 
{
  int i, j;

  for (i = k; i >= 0; i--)
    if BEATS(p.objectives[i],q.objectives[i]) 
      {for (j = i - 1; j >= 0; j--) 
         if BEATS(q.objectives[j],p.objectives[j]) return 0; 
       return -1;}
    else
    if BEATS(q.objectives[i],p.objectives[i]) 
      {for (j = i - 1; j >= 0; j--) 
         if BEATS(p.objectives[j],q.objectives[j]) return 0; 
       return  1;}
  return 2;
}


// bool, false, true
int dominates1way(POINT p, POINT q, int k)
// returns true if p dominates q or p == q, false otherwise 
// the assumption is that q doesn't dominate p 
// k is the highest index inspected 
{
  int i;

  for (i = k; i >= 0; i--)
    if BEATS(q.objectives[i],p.objectives[i]) 
      return 0;
  return 1;
}


void makeDominatedBit(FRONT ps, FRONT *fs, int *safe, int *fr, int p, int n)
// creates the front ps[0 .. p-1] in fs[fr], with each point bounded by ps[p] and dominated points removed 
{
  int l = 0;
  int u = p - 1;
  int i, j;

  for (i = p - 1; i >= 0; i--)
    if (BEATS(ps.points[p].objectives[n - 1],ps.points[i].objectives[n - 1]))
      {fs[*fr].points[u].objectives[n - 1] = ps.points[i].objectives[n - 1]; 
       for (j = 0; j < n - 1; j++) 
	 fs[*fr].points[u].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]); 
       u--;}
    else
      {fs[*fr].points[l].objectives[n - 1] = ps.points[p].objectives[n - 1]; 
       for (j = 0; j < n - 1; j++) 
 	 fs[*fr].points[l].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]); 
       l++;}
  POINT t;
  // points below l are all equal in the last objective; points above l are all worse 
  // points below l can dominate each other, and we don't need to compare the last objective 
  // points above l cannot dominate points that start below l, and we don't need to compare the last objective 
  fs[*fr].nPoints = 1;
  for (i = 1; i < l; i++)
     {j = 0;
      while (j < fs[*fr].nPoints)
	switch (dominates2way(fs[*fr].points[i], fs[*fr].points[j], n-2))
	  {case  0: j++; break;
	   case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
	            // SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
	            t = fs[*fr].points[j];
		    fs[*fr].points[j] = fs[*fr].points[i]; 
		    fs[*fr].points[i] = t; 
		    while(j < fs[*fr].nPoints - 1 && dominates1way(fs[*fr].points[j], fs[*fr].points[fs[*fr].nPoints - 1], n-1))
		      fs[*fr].nPoints--;
		    int k = j+1; 
		    while (k < fs[*fr].nPoints)
		      if(dominates1way(fs[*fr].points[j], fs[*fr].points[k], n-2))
			{t = fs[*fr].points[k];
			 fs[*fr].nPoints--;
			 fs[*fr].points[k] = fs[*fr].points[fs[*fr].nPoints]; 
			 fs[*fr].points[fs[*fr].nPoints] = t; 
			}
		      else
			k++;
	   default: j = fs[*fr].nPoints + 1;
	  }
      if (j == fs[*fr].nPoints) 
	{t = fs[*fr].points[fs[*fr].nPoints]; 
	 fs[*fr].points[fs[*fr].nPoints] = fs[*fr].points[i]; 
	 fs[*fr].points[i] = t; 
	 fs[*fr].nPoints++;}
     }
  *safe = WORSE(l,fs[*fr].nPoints);
  for (i = l; i < p; i++)
     {j = 0;
      while (j < *safe)
	if(dominates1way(fs[*fr].points[j], fs[*fr].points[i], n-2))
	  j = fs[*fr].nPoints + 1;
	else
	  j++;
      while (j < fs[*fr].nPoints)
	switch (dominates2way(fs[*fr].points[i], fs[*fr].points[j], n-1))
	  {case  0: j++; break;
	   case -1: // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
	            // SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
	            t = fs[*fr].points[j];
		    fs[*fr].points[j] = fs[*fr].points[i]; 
		    fs[*fr].points[i] = t; 
		    while(j < fs[*fr].nPoints - 1 && dominates1way(fs[*fr].points[j], fs[*fr].points[fs[*fr].nPoints - 1], n-1))
		      fs[*fr].nPoints--;
		    int k = j+1; 
		    while (k < fs[*fr].nPoints)
		      if(dominates1way(fs[*fr].points[j], fs[*fr].points[k], n-1))
			{t = fs[*fr].points[k];
			 fs[*fr].nPoints--;
			 fs[*fr].points[k] = fs[*fr].points[fs[*fr].nPoints]; 
			 fs[*fr].points[fs[*fr].nPoints] = t; 
			}
		      else
			k++;
	   default: j = fs[*fr].nPoints + 1;
	  }
      if (j == fs[*fr].nPoints) 
	{t = fs[*fr].points[fs[*fr].nPoints]; 
	  fs[*fr].points[fs[*fr].nPoints] = fs[*fr].points[i]; 
	  fs[*fr].points[i] = t; 
	  fs[*fr].nPoints++;}
     }
  (*fr)++;
}


double hv2_wfg(FRONT ps, int k)
// returns the hypervolume of ps[0 .. k-1] in 2D 
// assumes that ps is sorted improving
{
  int i;
  double volume = ps.points[0].objectives[0] * ps.points[0].objectives[1]; 

  for (i = 1; i < k; i++) 
    volume += ps.points[i].objectives[1] * 
              (ps.points[i].objectives[0] - ps.points[i - 1].objectives[0]);
  return volume;
}


double inclhv(POINT p, int n)
// returns the inclusive hypervolume of p
{
  double volume = 1;
  int i;

  for (i = 0; i < n; i++) 
    volume *= p.objectives[i];
  return volume;
}


double inclhv2(POINT p, POINT q, int n)
// returns the hypervolume of {p, q}
{
  double vp  = 1; double vq  = 1;
  double vpq = 1;
  int i;

  for (i = 0; i < n; i++) 
    {
      vp  *= p.objectives[i];
      vq  *= q.objectives[i];
      vpq *= WORSE(p.objectives[i],q.objectives[i]);
    }
  return vp + vq - vpq;
}


double inclhv3(POINT p, POINT q, POINT r, int n)
// returns the hypervolume of {p, q, r}
{
  double vp   = 1; double vq   = 1; double vr   = 1;
  double vpq  = 1; double vpr  = 1; double vqr  = 1;
  double vpqr = 1;
  int i;

  for (i = 0; i < n; i++) 
    {
      vp *= p.objectives[i];
      vq *= q.objectives[i];
      vr *= r.objectives[i];
      if (BEATS(p.objectives[i],q.objectives[i]))
	if (BEATS(q.objectives[i],r.objectives[i]))
	{
	  vpq  *= q.objectives[i];
	  vpr  *= r.objectives[i];
	  vqr  *= r.objectives[i];
	  vpqr *= r.objectives[i];
	}
	else
	{
	  vpq  *= q.objectives[i];
	  vpr  *= WORSE(p.objectives[i],r.objectives[i]);
	  vqr  *= q.objectives[i];
	  vpqr *= q.objectives[i];
	}
      else
	if (BEATS(p.objectives[i],r.objectives[i]))
	{
	  vpq  *= p.objectives[i];
	  vpr  *= r.objectives[i];
	  vqr  *= r.objectives[i];
	  vpqr *= r.objectives[i];
	}
	else
	{
	  vpq  *= p.objectives[i];
	  vpr  *= p.objectives[i];
	  vqr  *= WORSE(q.objectives[i],r.objectives[i]);
	  vpqr *= p.objectives[i];
	}
    }
  return vp + vq + vr - vpq - vpr - vqr + vpqr;
}


double inclhv4(POINT p, POINT q, POINT r, POINT s, int n)
// returns the hypervolume of {p, q, r, s}
{
  double vp    = 1; double vq   = 1; double vr   = 1; double vs   = 1;
  double vpq   = 1; double vpr  = 1; double vps  = 1; double vqr  = 1; double vqs = 1; double vrs = 1; 
  double vpqr  = 1; double vpqs = 1; double vprs = 1; double vqrs = 1; 
  double vpqrs = 1; 
  int i;

  for (i = 0; i < n; i++) 
    {
      vp *= p.objectives[i];
      vq *= q.objectives[i];
      vr *= r.objectives[i];
      vs *= s.objectives[i];
      if (BEATS(p.objectives[i],q.objectives[i]))
	if (BEATS(q.objectives[i],r.objectives[i]))
	  if (BEATS(r.objectives[i],s.objectives[i]))
	    {
	      vpq *= q.objectives[i];
	      vpr *= r.objectives[i];
	      vps *= s.objectives[i];
	      vqr *= r.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= s.objectives[i];
	      vpqr *= r.objectives[i];
	      vpqs *= s.objectives[i];
	      vprs *= s.objectives[i];
	      vqrs *= s.objectives[i];
	      vpqrs *= s.objectives[i];
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(q.objectives[i],s.objectives[i]);
	      vpq *= q.objectives[i];
	      vpr *= r.objectives[i];
	      vps *= WORSE(p.objectives[i],s.objectives[i]);
	      vqr *= r.objectives[i];
	      vqs *= z1;
	      vrs *= r.objectives[i];
	      vpqr *= r.objectives[i];
	      vpqs *= z1;
	      vprs *= r.objectives[i];
	      vqrs *= r.objectives[i];
	      vpqrs *= r.objectives[i];
	    }
	else
	  if (BEATS(q.objectives[i],s.objectives[i]))
	    {
	      vpq *= q.objectives[i];
	      vpr *= WORSE(p.objectives[i],r.objectives[i]);
	      vps *= s.objectives[i];
	      vqr *= q.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= s.objectives[i];
	      vpqr *= q.objectives[i];
	      vpqs *= s.objectives[i];
	      vprs *= s.objectives[i];
	      vqrs *= s.objectives[i];
	      vpqrs *= s.objectives[i];
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
	      vpq *= q.objectives[i];
	      vpr *= z1;
	      vps *= WORSE(p.objectives[i],s.objectives[i]);
	      vqr *= q.objectives[i];
	      vqs *= q.objectives[i];
	      vrs *= WORSE(r.objectives[i],s.objectives[i]);
	      vpqr *= q.objectives[i];
	      vpqs *= q.objectives[i];
	      vprs *= WORSE(z1,s.objectives[i]);
	      vqrs *= q.objectives[i];
	      vpqrs *= q.objectives[i];
	    }
      else
	if (BEATS(q.objectives[i],r.objectives[i]))
	  if (BEATS(p.objectives[i],s.objectives[i]))
	    {
	      OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
	      OBJECTIVE z2 = WORSE(r.objectives[i],s.objectives[i]);
	      vpq *= p.objectives[i];
	      vpr *= z1;
	      vps *= s.objectives[i];
	      vqr *= r.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= z2;
	      vpqr *= z1;
	      vpqs *= s.objectives[i];
	      vprs *= z2;
	      vqrs *= z2;
	      vpqrs *= z2;
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(p.objectives[i],r.objectives[i]);
	      OBJECTIVE z2 = WORSE(r.objectives[i],s.objectives[i]);
	      vpq *= p.objectives[i];
	      vpr *= z1;
	      vps *= p.objectives[i];
	      vqr *= r.objectives[i];
	      vqs *= WORSE(q.objectives[i],s.objectives[i]);
	      vrs *= z2;
	      vpqr *= z1;
	      vpqs *= p.objectives[i];
	      vprs *= z1;
	      vqrs *= z2;
	      vpqrs *= z1;
	    }
	else
	  if (BEATS(p.objectives[i],s.objectives[i]))
	    {
	      vpq *= p.objectives[i];
	      vpr *= p.objectives[i];
	      vps *= s.objectives[i];
	      vqr *= q.objectives[i];
	      vqs *= s.objectives[i];
	      vrs *= s.objectives[i];
	      vpqr *= p.objectives[i];
	      vpqs *= s.objectives[i];
	      vprs *= s.objectives[i];
	      vqrs *= s.objectives[i];
	      vpqrs *= s.objectives[i];
	    }
	  else
	    {
	      OBJECTIVE z1 = WORSE(q.objectives[i],s.objectives[i]);
	      vpq *= p.objectives[i];
	      vpr *= p.objectives[i];
	      vps *= p.objectives[i];
	      vqr *= q.objectives[i];
	      vqs *= z1;
	      vrs *= WORSE(r.objectives[i],s.objectives[i]);
	      vpqr *= p.objectives[i];
	      vpqs *= p.objectives[i];
	      vprs *= p.objectives[i];
	      vqrs *= z1;
	      vpqrs *= p.objectives[i];
	    }
    }
  return vp + vq + vr + vs - vpq - vpr - vps - vqr - vqs - vrs + vpqr + vpqs + vprs + vqrs - vpqrs;
}


double exclhv(FRONT ps, FRONT *fs, int *safe, int *fr, int p, int n)
// returns the exclusive hypervolume of ps[p] relative to ps[0 .. p-1] 
{

  makeDominatedBit(ps, fs, safe, fr, p, n);
  double volume = inclhv(ps.points[p], n) - hv_wfg(fs[*fr - 1], fs, safe, fr, n);
  (*fr)--;
  return volume;
}


double hv_wfg(FRONT ps, FRONT *fs, int *safe, int *fr, int n)
// returns the hypervolume of ps[0 ..] 
{
  int i;

  // process small fronts with the IEA 
  switch (ps.nPoints)
    {case 1: return inclhv (ps.points[0], n); 
     case 2: return inclhv2(ps.points[0], ps.points[1], n); 
     case 3: return inclhv3(ps.points[0], ps.points[1], ps.points[2], n); 
     case 4: return inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3], n); 
    }

  // these points need sorting 
  //RHG qsort(&ps.points[*safe], ps.nPoints - *safe, sizeof(POINT), greater); 
  EMO_quicksort(&(ps.points[*safe]), ps.nPoints - *safe, 0, n, sizeof(POINT), greater);

  // n = 2 implies that safe = 0 
  if (n == 2) return hv2_wfg(ps, ps.nPoints); 
  // these points don't NEED sorting, but it helps 
  //RHG qsort(ps.points, *safe, sizeof(POINT), greaterabbrev); 
  EMO_quicksort(ps.points, *safe, 0, n, sizeof(POINT), greaterabbrev);

  if (n == 3 && *safe > 0) 
    {
      double volume = ps.points[0].objectives[2] * hv2_wfg(ps, *safe); 
      n--;
      for (i = *safe; i < ps.nPoints; i++)
	// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit 
	volume += ps.points[i].objectives[n] * exclhv(ps, fs, safe, fr, i, n);
      n++; 
      return volume;
    }
  else
    {
      double volume = inclhv4(ps.points[0], ps.points[1], ps.points[2], ps.points[3], n); 
      n--;
      for (i = 4; i < ps.nPoints; i++)
	// we can ditch dominated points here, but they will be ditched anyway in makeDominatedBit 
	volume += ps.points[i].objectives[n] * exclhv(ps, fs, safe, fr, i, n);
      n++; 
      return volume;
    }
}

// RHG
void EMO_HV_alloc(EMO_HV *hv, int max_row, int col) {
  int i, j, maxdepth = col - 2;

  hv->fs = malloc(sizeof(FRONT) * maxdepth);

  for(i = 0; i < maxdepth; i++) {
    hv->fs[i].points = (POINT *) malloc(sizeof(POINT) * max_row);
    hv->fs[i].nPoints = max_row;

    for (j = 0; j < max_row; j++)
      hv->fs[i].points[j].objectives = (OBJECTIVE *) malloc(sizeof(OBJECTIVE) * (col - i - 1));
  }

  hv->f.points = (POINT *) malloc(sizeof(POINT) * max_row);

  for(i = 0; i < max_row; i++)
    hv->f.points[i].objectives = (OBJECTIVE *) malloc(sizeof(OBJECTIVE) * col);

  hv->row = max_row;
  hv->col = col;
}

void EMO_HV_free(EMO_HV *hv) {
  int i, j, maxdepth = hv->col - 2;

  for(i = 0; i < hv->row; i++)
    free(hv->f.points[i].objectives);

  free(hv->f.points);

  for (i = 0; i < maxdepth; i++) {
    for (j = 0; j < hv->row; j++)
       free(hv->fs[i].points[j].objectives);

    free(hv->fs[i].points);
  }

  free(hv->fs);
}

double EMO_HV_run(EMO_HV *hv, double *data, int *enable, int row, const double *ref) {
  int i, j, k, c;

  if(row > hv->row) {
    printf("Error, mismatch dimensions of row in EMO_HV_run (%d vs %d).\n", hv->row, row);
    exit(1);
  }

  hv->fr = hv->safe = 0;
  c = 0;

  for(k = i = 0; i < row; i++) {
    if(enable == NULL || enable[i]) {
      for(j = 0; j <  hv->col; j++) {
        // modify the objective values relative to the reference point 
        if(data[i* hv->col + j] >= ref[j]) {
          c++;
          k--;
          break;
        }
        else {
          hv->f.points[k].objectives[j] = fabs(data[i* hv->col + j] - ref[j]);
        }
      }
      k++;
    }
  }

  if(c > 0)
    printf("Warning: omitting %d/%d points that are outside the reference point (total %d)\n", c, row, k);

  if(k == 0) return 0;

  hv->f.nPoints = k;

  //double hv_wfg(FRONT ps, FRONT *fs, int *safe, int *fr, int n)
  return hv_wfg(hv->f, hv->fs, &(hv->safe), &(hv->fr), hv->col);
}


/* Compute the normalized hypervolume */
double EMO_HV_run2(EMO_HV *hv, double *data, int *enable, int row, const double *ref) {
  double v = 1.0;
  int i;

  for(i = hv->col - 1; i > -1; i--)
    v *= ref[i];

  if(v == 0) {
    printf("Error, one coordinate of the reference point is zero in hv_wfg.c.\n");
    exit(1);
  }
  return EMO_HV_run(hv, data, enable, row, ref) / v;
}

/* Hipervolume contribution, it only considers enabled solutions, RHG */
double EMO_HV_contribution(EMO_HV *hv, double *deltahv, double *data, int *enable, int row, const double *ref, int col) {
  double totalhv;
  int i;

  totalhv = EMO_HV_run(hv, data, enable, row, ref);

  for(i = 0; i < row; i++) { 
    if(enable[i]) {
      enable[i] = 0;
      deltahv[i] = totalhv - EMO_HV_run(hv, data, enable, row, ref); 
      //printf("%d: hv %f\n", i, deltahv[i]);
      enable[i] = 1;
    }
    else {
      deltahv[i] = DBL_MAX;
    }
  }

  return totalhv;
}

