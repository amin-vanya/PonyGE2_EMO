/* hv_iwfg.c  incremental IWFG algorithm (see \cite{Cox16})
 *
 * Lyndon While, Lucas Bradstreet, Wesley Cox 
 * Email lyndon.while@uwa.edu.au for any queries regarding usage/bugs/improvements. 
 *
 * This code includes a high performance implementation of the IWFG algorithm, 
 * used to identify the smallest hypervolume-contributor of a set of non-dominated points. 
 *
 * COPYRIGHT: 
 * This software is Copyright (C) 2015 Lyndon While, Lucas Bradstreet, Wesley Cox. 
 * This program is free software (software libre). You can redistribute it and/or modify it under 
 * the terms of the GNU General Public License as published by the Free Software Foundation; 
 * either version 2 of the License, or (at your option) any later version. 
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; 
 * without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
 * See the GNU General Public License for more details. 
 *
 * Taken from version: IWFG_1.01 (24 November 2015)
 *
 * http://www.wfg.csse.uwa.edu.au/hypervolume/#code
 *
 * Updates:
 * 
 * 26 Jul 2016  Angel Cesar Lara Bernal UPIITA-IPN
 *              Integration with EMO Project, global variables are removed.
 *
 * 20 May 2017  Raquel Hernandez Gomez
 *              Fixed one bug:  binarySearch returns -1 when a pointer is not found.
 *
 */

#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include "vector.h"
#include "hv_iwfg.h"
#include "hv_sort.h"

#define BEATS(x,y) (x > y) 
#define WORSE(x,y) (BEATS(y,x) ? (x) : (y))
#define MIN(a,b)   (a < b ? (a) : (b))
#define MAX(a,b)   (a > b ? (a) : (b))
#define SLICELIMIT 5

double fhv(EMO_IWFG *hv, FRONT ps);


int sorter(const void *a, const void *b, int start, EMO_IWFG *hv)
//sort point indexes into groups of last objective
{
   int i = *(int*)a;
   int j = *(int*)b;

   if(hv->torder[i][hv->n-1] == hv->torder[j][hv->n-1]) 
      return hv->tcompare[j][hv->torder[j][hv->n-1]] - hv->tcompare[i][hv->torder[i][hv->n-1]];
   else 
      return hv->torder[i][hv->n-1] - hv->torder[j][hv->n-1];
}


int iwfg_greater(const void *v1, const void *v2, int start, EMO_IWFG *hv)
// this sorts points worsening in the last objective
{
   int i;

   POINT p = *(POINT*)v1;
   POINT q = *(POINT*)v2;

   for (i = hv->n - 1; i >= start; i--)
      if BEATS(p.objectives[i],q.objectives[i]) 
         return -1;
      else if BEATS(q.objectives[i],p.objectives[i]) 
         return  1;

   return 0;
}


int greaterorder(const void *v1, const void *v2, int start, EMO_IWFG *hv)
// this sorts points worsening in the last objective for a certain objective ordering
{
   int i;

   POINT p = *(POINT*)v1;
   POINT q = *(POINT*)v2;

   for (i = hv->n - 1; i >= start; i--)
      if BEATS(p.objectives[hv->gorder[i]],q.objectives[hv->gorder[i]]) 
         return -1;
      else if BEATS(q.objectives[hv->gorder[i]],p.objectives[hv->gorder[i]]) 
         return  1;

   return 0;
}


int iwfg_greaterabbrev(const void *v1, const void *v2, int start, EMO_IWFG *hv)
// this sorts points worsening in the penultimate objective
{
   int i;

   POINT p = *(POINT*)v1;
   POINT q = *(POINT*)v2;

   for (i = hv->n - 2; i >= start; i--)
      if BEATS(p.objectives[i],q.objectives[i]) 
         return -1;
      else if BEATS(q.objectives[i],p.objectives[i]) 
         return  1;

   return 0;
}


int greaterabbrevorder(const void *v1, const void *v2, int start, EMO_IWFG *hv)
// this sorts points worsening in the penultimate objective for a certain objective ordering
{
   int i;

   POINT p = *(POINT*)v1;
   POINT q = *(POINT*)v2;

   for (i = hv->n - 2; i >= start; i--)
      if BEATS(p.objectives[hv->gorder[i]],q.objectives[hv->gorder[i]]) 
         return -1;
      else if BEATS(q.objectives[hv->gorder[i]],p.objectives[hv->gorder[i]]) 
         return  1;

   return 0;
}


int iwfg_dominates2way(POINT p, POINT q, int k)
// returns -1 if p dominates q, 1 if q dominates p, 2 if p == q, 0 otherwise 
// k is the highest index inspected 
{
   int i, j;

   for (i = k; i >= 0; i--)
      if BEATS(p.objectives[i],q.objectives[i])
      {
         for (j = i - 1; j >= 0; j--)
            if BEATS(q.objectives[j],p.objectives[j]) 
               return 0; 
         return -1;
      }
      else if BEATS(q.objectives[i],p.objectives[i])
      {
         for (j = i - 1; j >= 0; j--)
            if BEATS(p.objectives[j],q.objectives[j])
               return 0;
         return  1;
      }

   return 2;
}


int iwfg_dominates1way(POINT p, POINT q, int k)
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


int iwfg_dominates1wayOrder(POINT p, POINT q, int k, int* order)
// returns true if p dominates q or p == q, false otherwise 
// the assumption is that q doesn't dominate p 
// k is the highest index inspected 
{
   int i;

   for (i = k; i >= 0; i--)
      if BEATS(q.objectives[order[i]],p.objectives[order[i]]) 
         return 0;

   return 1;
}


void removeDominated(EMO_IWFG *hv, int l, int limit)
{
   int i, j, k;
   POINT t;

   // points below l are all equal in the last objective; points above l are all worse 
   // points below l can dominate each other, and we don't need to compare the last objective 
   // points above l cannot dominate points that start below l, and we don't need to compare the last objective 

   hv->fs[hv->fr].nPoints = 1;

   for (i = 1; i < l; i++)
   {
      j = 0;
      while (j < hv->fs[hv->fr].nPoints)
      {
         switch (iwfg_dominates2way(hv->fs[hv->fr].points[i], hv->fs[hv->fr].points[j], hv->n-2))
         {
            case  0: 
               j++; 
               break;
            case -1:
               // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
               // SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
               t = hv->fs[hv->fr].points[j];
               hv->fs[hv->fr].points[j] = hv->fs[hv->fr].points[i]; 
               hv->fs[hv->fr].points[i] = t; 
               while(j < hv->fs[hv->fr].nPoints - 1 && iwfg_dominates1way(hv->fs[hv->fr].points[j], hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints - 1], hv->n-1)) 
                  hv->fs[hv->fr].nPoints--;
               k = j+1; 
               while (k < hv->fs[hv->fr].nPoints)
                  if(iwfg_dominates1way(hv->fs[hv->fr].points[j], hv->fs[hv->fr].points[k], hv->n-2))
                  {
                     t = hv->fs[hv->fr].points[k];
                     hv->fs[hv->fr].nPoints--;
                     hv->fs[hv->fr].points[k] = hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints]; 
                     hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints] = t; 
                  }
                  else
                     k++;
            default: 
               j = hv->fs[hv->fr].nPoints + 1;
         }
      }
      if (j == hv->fs[hv->fr].nPoints)
      {
         t = hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints]; 
         hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints] = hv->fs[hv->fr].points[i]; 
         hv->fs[hv->fr].points[i] = t; 
         hv->fs[hv->fr].nPoints++;
      }
   }

   hv->safe = WORSE(l,hv->fs[hv->fr].nPoints);

   for (i = l; i < limit; i++)
   {
      j = 0;
      while(j < hv->safe)
         if(iwfg_dominates1way(hv->fs[hv->fr].points[j], hv->fs[hv->fr].points[i], hv->n-2))
            j = hv->fs[hv->fr].nPoints + 1;
         else
            j++;
      while(j < hv->fs[hv->fr].nPoints)
      {
         switch(iwfg_dominates2way(hv->fs[hv->fr].points[i], hv->fs[hv->fr].points[j], hv->n-1))
         {
            case  0: 
               j++; 
               break;
            case -1:
               // AT THIS POINT WE KNOW THAT i CANNOT BE DOMINATED BY ANY OTHER PROMOTED POINT j 
               // SWAP i INTO j, AND 1-WAY DOM FOR THE REST OF THE js 
               t = hv->fs[hv->fr].points[j];
               hv->fs[hv->fr].points[j] = hv->fs[hv->fr].points[i]; 
               hv->fs[hv->fr].points[i] = t; 
               while(j < hv->fs[hv->fr].nPoints - 1 && iwfg_dominates1way(hv->fs[hv->fr].points[j], hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints - 1], hv->n-1))
                   hv->fs[hv->fr].nPoints--;
               k = j+1; 
               while (k < hv->fs[hv->fr].nPoints)
                  if(iwfg_dominates1way(hv->fs[hv->fr].points[j], hv->fs[hv->fr].points[k], hv->n-1))
                  {
                     t = hv->fs[hv->fr].points[k];
                     hv->fs[hv->fr].nPoints--;
                     hv->fs[hv->fr].points[k] = hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints]; 
                     hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints] = t; 
                  }
                  else
                     k++;
            default: 
               j = hv->fs[hv->fr].nPoints + 1;
         }
      }
      if (j == hv->fs[hv->fr].nPoints)
      {
         t = hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints]; 
         hv->fs[hv->fr].points[hv->fs[hv->fr].nPoints] = hv->fs[hv->fr].points[i]; 
         hv->fs[hv->fr].points[i] = t; 
         hv->fs[hv->fr].nPoints++;
      }
   }

   hv->fr++;
}


void iwfg_makeDominatedBit(EMO_IWFG *hv, FRONT ps, int p)
// creates the front ps[0 .. p-1] in fs[fr], with each point bounded by ps[p] and dominated points removed 
{
   int i, j;
   int l = 0;
   int u = p - 1;

   for (i = p - 1; i >= 0; i--)
      if (BEATS(ps.points[p].objectives[hv->n - 1],ps.points[i].objectives[hv->n - 1]))
      {
         hv->fs[hv->fr].points[u].objectives[hv->n - 1] = ps.points[i].objectives[hv->n - 1]; 
         for(j = 0; j < hv->n - 1; j++)
            hv->fs[hv->fr].points[u].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]); 
         u--;
      }
      else
      {
         hv->fs[hv->fr].points[l].objectives[hv->n - 1] = ps.points[p].objectives[hv->n - 1]; 
         for (j = 0; j < hv->n - 1; j++) 
            hv->fs[hv->fr].points[l].objectives[j] = WORSE(ps.points[p].objectives[j],ps.points[i].objectives[j]); 
         l++;
      }

   removeDominated(hv, l, p);
}


double hv2(FRONT ps, int k)
// returns the hypervolume of ps[0 .. k-1] in 2D 
// assumes that ps is sorted improving
{
   int i;

   double volume = ps.points[0].objectives[0] * ps.points[0].objectives[1]; 

   for (i = 1; i < k; i++) 
      volume += ps.points[i].objectives[1] * (ps.points[i].objectives[0] - ps.points[i - 1].objectives[0]);
	
   return volume;
}


double iwfg_inclhv(EMO_IWFG *hv, POINT p)
// returns the inclusive hypervolume of p
{
   int i;

   double volume = 1;

   for (i = 0; i < hv->n; i++) 
      volume *= p.objectives[i];
	
   return volume;
}


double inclhvOrder(EMO_IWFG *hv, POINT p, int* order)
// returns the inclusive hypervolume of p
{
   double volume = 1;
   int i;

   for (i = 0; i < hv->n; i++) 
      volume *= p.objectives[order[i]];
	
   return volume;
}


double iwfg_inclhv2(EMO_IWFG *hv, POINT p, POINT q)
// returns the hypervolume of {p, q}
{
   double vp  = 1; 
   double vq  = 1;
   double vpq = 1;
   int i;

   for (i = 0; i < hv->n; i++)
   {
      vp  *= p.objectives[i];
      vq  *= q.objectives[i];
      vpq *= WORSE(p.objectives[i],q.objectives[i]);
   }

   return vp + vq - vpq;
}


double iwfg_inclhv3(EMO_IWFG *hv, POINT p, POINT q, POINT r)
// returns the hypervolume of {p, q, r}
{
   double vp = 1;
   double vq = 1; 
   double vr = 1;
   double vpq = 1; 
   double vpr = 1; 
   double vqr = 1;
   double vpqr = 1;
   int i;

   for (i = 0; i < hv->n; i++)
   {
      vp *= p.objectives[i];
      vq *= q.objectives[i];
      vr *= r.objectives[i];
      if (BEATS(p.objectives[i],q.objectives[i]))
      {
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
      }
      else if (BEATS(p.objectives[i],r.objectives[i]))
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


double iwfg_inclhv4(EMO_IWFG *hv, POINT p, POINT q, POINT r, POINT s)
// returns the hypervolume of {p, q, r, s}
{
   double vp = 1; 
   double vq = 1; 
   double vr = 1; 
   double vs = 1;
   double vpq = 1; 
   double vpr = 1; 
   double vps = 1; 
   double vqr = 1; 
   double vqs = 1; 
   double vrs = 1; 
   double vpqr = 1; 
   double vpqs = 1; 
   double vprs = 1; 
   double vqrs = 1; 
   double vpqrs = 1; 
   int i;

   for (i = 0; i < hv->n; i++)
   {
      vp *= p.objectives[i];
      vq *= q.objectives[i];
      vr *= r.objectives[i];
      vs *= s.objectives[i];
      if (BEATS(p.objectives[i],q.objectives[i]))
      {
         if (BEATS(q.objectives[i],r.objectives[i]))
         {
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
         }
         else if (BEATS(q.objectives[i],s.objectives[i]))
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
      }
      else if (BEATS(q.objectives[i],r.objectives[i])) {
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
      }
      else if (BEATS(p.objectives[i],s.objectives[i]))
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


/* Not used
double inclhv5(POINT p, POINT q, POINT r, POINT s, POINT t)
// returns the hypervolume of {p, q, r, s, t}
{
   double vp = 1; 
   double vq = 1; 
   double vr = 1; 
   double vs = 1;
   double vt = 1;

   double vpq = 1; 
   double vpr = 1; 
   double vps = 1; 
   double vpt = 1; 
   double vqr = 1; 
   double vqs = 1; 
   double vqt = 1; 
   double vrs = 1; 
   double vrt = 1; 
   double vst = 1; 

   double vpqr = 1; 
   double vpqs = 1; 
   double vpqt = 1; 
   double vprs = 1; 
   double vprt = 1; 
   double vpst = 1; 
   double vqrs = 1; 
   double vqrt = 1; 
   double vqst = 1; 
   double vrst = 1; 

   double vpqrs = 1; 
   double vpqrt = 1; 
   double vpqst = 1; 
   double vprst = 1; 
   double vqrst = 1; 

   double vpqrst = 1; 
   int i;

   for (i = 0; i < n; i++)
   {
      vp *= p.objectives[i];
      vq *= q.objectives[i];
      vr *= r.objectives[i];
      vs *= s.objectives[i];
      vt *= t.objectives[i];
      vpq *= WORSE(p.objectives[i],q.objectives[i]);
      vpr *= WORSE(p.objectives[i],r.objectives[i]);
      vps *= WORSE(p.objectives[i],s.objectives[i]);
      vpt *= WORSE(p.objectives[i],t.objectives[i]);
      vqr *= WORSE(q.objectives[i],r.objectives[i]);
      vqs *= WORSE(q.objectives[i],s.objectives[i]);
      vqt *= WORSE(q.objectives[i],t.objectives[i]);
      vrs *= WORSE(r.objectives[i],s.objectives[i]);
      vrt *= WORSE(r.objectives[i],t.objectives[i]);
      vst *= WORSE(s.objectives[i],t.objectives[i]);
      vpqr *= WORSE(p.objectives[i],WORSE(q.objectives[i],r.objectives[i]));
      vpqs *= WORSE(p.objectives[i],WORSE(q.objectives[i],s.objectives[i]));
      vpqt *= WORSE(p.objectives[i],WORSE(q.objectives[i],t.objectives[i]));
      vprs *= WORSE(p.objectives[i],WORSE(r.objectives[i],s.objectives[i]));
      vprt *= WORSE(p.objectives[i],WORSE(r.objectives[i],t.objectives[i]));
      vpst *= WORSE(p.objectives[i],WORSE(s.objectives[i],t.objectives[i]));
      vqrs *= WORSE(q.objectives[i],WORSE(r.objectives[i],s.objectives[i]));
      vqrt *= WORSE(q.objectives[i],WORSE(r.objectives[i],t.objectives[i]));
      vqst *= WORSE(q.objectives[i],WORSE(s.objectives[i],t.objectives[i]));
      vrst *= WORSE(r.objectives[i],WORSE(s.objectives[i[,t.objectives[i]));
      vpqrs *= WORSE(WORSE(p.objectives[i],q.objectives[i]),WORSE(r.objectives[i],s.objectives[i]));
      vpqrt *= WORSE(WORSE(p.objectives[i],q.objectives[i]),WORSE(r.objectives[i],t.objectives[i]));
      vpqst *= WORSE(WORSE(p.objectives[i],q.objectives[i]),WORSE(s.objectives[i],t.objectives[i]));
      vprst *= WORSE(WORSE(p.objectives[i],r.objectives[i]),WORSE(s.objectives[i],t.objectives[i]));
      vqrst *= WORSE(WORSE(q.objectives[i],r.objectives[i]),WORSE(s.objectives[i],t.objectives[i]));
      vpqrst *= WORSE(WORSE(p.objectives[i],WORSE(q.objectives[i],r.objectives[i])),WORSE(s.objectives[i],t.objectives[i]));
   }

   return vp + vq + vr + vs + vt 
   - vpq  - vpr  - vps  - vpt  - vqr  - vqs  - vqt  - vrs  - vrt  - vst 
   + vpqr + vpqs + vpqt + vprs + vprt + vpst + vqrs + vqrt + vqst + vrst 
   - vpqrs - vpqrt - vpqst - vprst - vqrst 
   + vpqrst;
}
*/


double iwfg_exclhv(EMO_IWFG *hv, FRONT ps, int p)
// returns the exclusive hypervolume of ps[p] relative to ps[0 .. p-1] 
{
   double volume;

   iwfg_makeDominatedBit(hv, ps, p);
   volume = iwfg_inclhv(hv, ps.points[p]) - fhv(hv, hv->fs[hv->fr - 1]);
   hv->fr--;

   return volume;
}


double fhv(EMO_IWFG *hv, FRONT ps)
// returns the hypervolume of ps[0 ..] 
{
   double volume;
   int i;

   // process small fronts with the IEA 
   switch (ps.nPoints)
   {
      case 1: 
         return iwfg_inclhv(hv, ps.points[0]); 
      case 2: 
         return iwfg_inclhv2(hv, ps.points[0], ps.points[1]); 
      case 3: 
         return iwfg_inclhv3(hv, ps.points[0], ps.points[1], ps.points[2]); 
      case 4: 
         return iwfg_inclhv4(hv, ps.points[0], ps.points[1], ps.points[2], ps.points[3]); 
   }

   // these points need sorting
   IWFG_quicksort(&ps.points[hv->safe], ps.nPoints - hv->safe, 0, hv, sizeof(POINT), iwfg_greater); 

   // n = 2 implies that safe = 0 
   if (hv->n == 2) 
      return hv2(ps, ps.nPoints); 

   // these points don't NEED sorting, but it helps 
   //qsort(ps.points, hv->safe, sizeof(POINT), greaterabbrev); 
   IWFG_quicksort(ps.points, hv->safe, 0, hv, sizeof(POINT), iwfg_greaterabbrev);

   if (hv->n == 3 && hv->safe > 0)
   {
      volume = ps.points[0].objectives[2] * hv2(ps, hv->safe); 
      hv->n--;
      for (i = hv->safe; i < ps.nPoints; i++)
         // we can ditch dominated points here, but they will be ditched anyway in iwfg_makeDominatedBit 
         volume += ps.points[i].objectives[hv->n] * iwfg_exclhv(hv, ps, i);
      hv->n++; 
      return volume;
   }
   else
   {
      volume = iwfg_inclhv4(hv, ps.points[0], ps.points[1], ps.points[2], ps.points[3]); 
      hv->n--;
      for (i = 4; i < ps.nPoints; i++)
         // we can ditch dominated points here, but they will be ditched anyway in hv, iwfg_makeDominatedBit 
         volume += ps.points[i].objectives[hv->n] * iwfg_exclhv(hv, ps, i);
      hv->n++; 
      return volume;
   }
}

void iwfg_makeDominatedBitPoint(EMO_IWFG *hv, FRONT ps, POINT p, int* order)
// creates the front ps in fs[fr], with each point bounded by p and dominated points removed 
{
   int i, j;

   int l = 0;
   int u = ps.nPoints - 1;
   for (i = ps.nPoints - 1; i >= 0; i--)
      if (BEATS(p.objectives[order[hv->n - 1]],ps.points[i].objectives[order[hv->n - 1]]))
      {
         hv->fs[hv->fr].points[u].objectives[hv->n - 1] = ps.points[i].objectives[order[hv->n - 1]]; 
         for (j = 0; j < hv->n - 1; j++)
            hv->fs[hv->fr].points[u].objectives[j] = WORSE(p.objectives[order[j]],ps.points[i].objectives[order[j]]); 
         u--;
      }
      else
      {
         hv->fs[hv->fr].points[l].objectives[hv->n - 1] = p.objectives[order[hv->n - 1]]; 
         for (j = 0; j < hv->n - 1; j++)
            hv->fs[hv->fr].points[l].objectives[j] = WORSE(p.objectives[order[j]],ps.points[i].objectives[order[j]]); 
         l++;
      }

   removeDominated(hv, l, ps.nPoints);
}


double iwfg_exclhvPoint(EMO_IWFG *hv, FRONT ps, POINT p, int* order)
// returns the exclusive hypervolume of p relative to ps
{
   double volume;

   iwfg_makeDominatedBitPoint(hv, ps, p, order);
   volume = inclhvOrder(hv, p, order) - fhv(hv, hv->fs[hv->fr - 1]);
   hv->fr--;

   return volume;
}


void heapify(EMO_IWFG *hv, int location, int index)
// restores heap property starting at location and working downwards to place index in heap
{
   int left, right;

   while(2*location+2 < hv->heapsize)
   {
      left = 0;
      right = 0;

      if (hv->partial[hv->heap[2*location+1]] < hv->partial[index]) 
         left = 1;
      if (hv->partial[hv->heap[2*location+2]] < hv->partial[index]) 
         right = 1;
      if (left) {
         if (right && hv->partial[hv->heap[2*location+2]] < hv->partial[hv->heap[2*location+1]])
         {
            hv->heap[location] = hv->heap[2*location+2];
            location = 2*location+2;
         }
         else
         {
            hv->heap[location] = hv->heap[2*location+1];
            location = 2*location+1;
         }
      }
      else if (right)
      {
         hv->heap[location] = hv->heap[2*location+2];
         location = 2*location+2;
      }
      else 
         break;
   }
   if (2*location+1<hv->heapsize && hv->partial[hv->heap[2*location+1]] < hv->partial[index])
   {
      hv->heap[location] = hv->heap[2*location+1];
      location = 2*location+1;
   }
   hv->heap[location] = index;
}


int peekFromHeap(EMO_IWFG *hv)
{
   return hv->heap[0];
}


void initialiseHeap(EMO_IWFG *hv, int capacity)
// creates the heap with the indexes 0..(capacity-1) 
{
   int i;

   hv->heapsize = capacity;

   for (i=hv->heapsize-1; i>=0; i--) 
      heapify(hv, i, i);
}


void insert(EMO_IWFG *hv, POINT p, int k, FRONT pl, int i, int j, int* order)
// inserts p into pl with the result in stacks[i][j]
{
   int place, placeNext;

   place = 0;
   while(place < pl.nPoints && pl.points[place].objectives[order[k]] > p.objectives[order[k]])
   {
      hv->stacks[i][j].front.points[place] = pl.points[place];
      place++;
   }
   POINT pp = pl.points[place];
   hv->stacks[i][j].front.points[place] = p;
   placeNext = place + 1;
   POINT ppn = pl.points[place+1];

   while (place < pl.nPoints)
   {
      if (!iwfg_dominates1wayOrder(p,pp,k,order))
      {
         hv->stacks[i][j].front.points[placeNext] = pp;
         placeNext++;
      }	
      place++;
      pp = ppn;
      ppn = pl.points[place+1];
   }
   hv->stacks[i][j].front.nPoints = placeNext;
}


void sliceOrder(EMO_IWFG *hv, int nPoints)
// slice using a separate objective ordering per point
{
   int i, p, j, k, l, seen, pos, start, end;

   int sorder[nPoints];

   for (i=0; i<nPoints; i++) 
      sorder[i] = i;

   IWFG_quicksort(sorder, nPoints, 0, hv, sizeof(int), sorter);

   seen = 0;

   for (p=0; p<nPoints; p++)
   {
      i = sorder[p];
      if (p==0 || hv->torder[i][hv->n-1]!=hv->torder[sorder[p-1]][hv->n-1])
      {
         seen = 0;
         hv->stacks[i][1].front.nPoints = 0;
      }
      else 
      {
         for(j=0; j<hv->stacks[sorder[p-1]][1].front.nPoints; j++)
            hv->stacks[i][1].front.points[j] = hv->stacks[sorder[p-1]][1].front.points[j];
         hv->stacks[i][1].front.nPoints = hv->stacks[sorder[p-1]][1].front.nPoints;
      }
      pos = nPoints - 1 - hv->tcompare[i][hv->torder[i][hv->n-1]];

      for(j=seen; j<pos; j++) {
         hv->stacks[i][1].front.points[hv->stacks[i][1].front.nPoints+j-seen] = hv->stacks[i][0].front.points[j];
      }

      start = hv->stacks[i][1].front.nPoints;
      end = hv->stacks[i][1].front.nPoints+pos-seen;
      seen = pos;
      POINT temp;

      for(j=start; j<end; j++)
      {
         k = 0;

         while(k < hv->stacks[i][1].front.nPoints)
            if(iwfg_dominates1wayOrder(hv->stacks[i][1].front.points[j],hv->stacks[i][1].front.points[k],hv->n-2,hv->torder[i]))
            {
               temp = hv->stacks[i][1].front.points[k];
               hv->stacks[i][1].front.points[k] = hv->stacks[i][1].front.points[j];
               hv->stacks[i][1].front.points[j] = temp;
               while(k < hv->stacks[i][1].front.nPoints-1 && 
               iwfg_dominates1wayOrder(hv->stacks[i][1].front.points[k],hv->stacks[i][1].front.points[hv->stacks[i][1].front.nPoints-1],hv->n-2,hv->torder[i]))
                  hv->stacks[i][1].front.nPoints--;
               l = k+1; 

               while(l < hv->stacks[i][1].front.nPoints)
                  if(iwfg_dominates1wayOrder(hv->stacks[i][1].front.points[k],hv->stacks[i][1].front.points[l],hv->n-2,hv->torder[i]))
                  {
                     temp = hv->stacks[i][1].front.points[l];
                     hv->stacks[i][1].front.nPoints--;
                     hv->stacks[i][1].front.points[l] = hv->stacks[i][1].front.points[hv->stacks[i][1].front.nPoints]; 
                     hv->stacks[i][1].front.points[hv->stacks[i][1].front.nPoints] = temp; 
                  }
                  else
                     l++;

               k = hv->stacks[i][1].front.nPoints + 1;
            }
            else
               k++;
         
         if(k == hv->stacks[i][1].front.nPoints)
         {
            temp = hv->stacks[i][1].front.points[hv->stacks[i][1].front.nPoints];
            hv->stacks[i][1].front.points[hv->stacks[i][1].front.nPoints] = hv->stacks[i][1].front.points[j];
            hv->stacks[i][1].front.points[j] = temp;
            hv->stacks[i][1].front.nPoints++;
         }
      }

      hv->stacks[i][1].index = pos+1;

      if(pos<nPoints-1) {
         hv->stacks[i][1].width = fabs(hv->stacks[i][0].front.points[pos].objectives[hv->torder[i][hv->n-1]]-
                                  hv->stacks[i][0].front.points[pos+1].objectives[hv->torder[i][hv->n-1]]);
      }
      else {
         hv->stacks[i][1].width = hv->stacks[i][0].front.points[pos].objectives[hv->torder[i][hv->n-1]];
      }
   }

   for(i=0; i<nPoints; i++)
   {
      hv->gorder = hv->torder[i];
      IWFG_quicksort(hv->stacks[i][1].front.points, hv->stacks[i][1].front.nPoints, 0, hv, sizeof(POINT), greaterabbrevorder);
   }
}


/* Not used
void slice(FRONT pl)
// slice in the last objective
{
   int i;

   stacks[0][1].front.nPoints = 0;

   for (i=0; i<pl.nPoints-1; i++)
   {
      stacks[i][1].width = fabs(pl.points[i].objectives[n-1]-pl.points[i+1].objectives[n-1]);
      stacks[i][1].index = i+1;
      insert(pl.points[i],n-2,stacks[i][1].front,i+1,1,torder[i]);
   }

   stacks[pl.nPoints-1][1].width = pl.points[pl.nPoints-1].objectives[n-1];
   stacks[pl.nPoints-1][1].index = pl.nPoints;
}
*/


int binarySearch(EMO_IWFG *hv, POINT p, int d)
{
   int min, max, mid, r, i;

   min = 0;
   max = hv->fsorted[d].nPoints-1;
   hv->gorder = hv->torder[d];

   while (min<=max)
   {
      mid = (max+min)/2;

      if(p.objectives == hv->fsorted[d].points[mid].objectives) {
         return mid;
      }
      else if((r = greaterorder(&p, &hv->fsorted[d].points[mid], 0, hv)) == -1) {
         max = mid-1;
      }
      else if(r == 1) {
         min = mid+1;
      }
      else { /* r == 0 */

        i = mid - 1;

        while(i >= min && greaterorder(&p, &hv->fsorted[d].points[i], 0, hv) == 0) {
          if(p.objectives == hv->fsorted[d].points[i].objectives) {
            return i;
          }
          i--;
        }
 
        i = mid + 1;

        while(i <= max && greaterorder(&p, &hv->fsorted[d].points[i], 0, hv) == 0) {
          if(p.objectives == hv->fsorted[d].points[i].objectives) {
            return i;
          }
          i++;
        }
     }
   }

   return -1;
}


void runHeuristic(EMO_IWFG *hv, FRONT ps)
{
   int i, j, k, x;

   for (i=0; i<hv->n-1; i++)
   {
      hv->torder[i][hv->n-1] = i;
      hv->torder[i][i] = hv->n-1;
   }

   for (i=hv->n-1; i>=0; i--)
   {
      for (j=0; j<ps.nPoints; j++)
      {
         hv->fsorted[i].points[j] = ps.points[j];
         hv->tcompare[j][i] = 0;
      }
      hv->fsorted[i].nPoints = ps.nPoints;
      hv->gorder = hv->torder[i];

      IWFG_quicksort(hv->fsorted[i].points, ps.nPoints, 0, hv, sizeof(POINT), greaterorder);
   }

   for (i=0; i<ps.nPoints; i++) { 
      for (k=0; k<hv->n; k++) {
         hv->tcompare[i][k] = ps.nPoints - 1 - binarySearch(hv, ps.points[i], k);
      }
   }

   for (i=0; i<ps.nPoints; i++)
      for (j=1; j<hv->n; j++)
      {
         x = hv->torder[i][j];
         k = j;

         while (k>0 && hv->tcompare[i][x] < hv->tcompare[i][hv->torder[i][k-1]])
         {
            hv->torder[i][k] = hv->torder[i][k-1];
            k--;
         }
         hv->torder[i][k] = x;
      }
}


int slicingDepth(int d)
{
   if (d <=  5) return 1;
   if (d <=  7) return 2;
   if (d <= 12) return 3;

   return 4;
}


int ihv2(EMO_IWFG *hv, FRONT ps, double *min)
{
   double vol, kvol;
   int k, sm;

   // returns the minimum exclusive hypervolume of points in ps for the 2D case
   
   IWFG_quicksort(ps.points, ps.nPoints, 0, hv, sizeof(POINT), iwfg_greater);

   vol = ps.points[0].objectives[0] * (ps.points[0].objectives[1] - ps.points[1].objectives[1]);
   sm = 0;

   for(k = 1; k<ps.nPoints-1; k++)
   {
      kvol = (ps.points[k].objectives[0] - ps.points[k-1].objectives[0]) * (ps.points[k].objectives[1] - ps.points[k+1].objectives[1]);
      if (kvol <= 0)
      {
         min[0] = ps.points[k].objectives[0];
         min[1] = ps.points[k].objectives[1];
         min[2] = 0;
         return ps.points[k].idx;
      }
      else if (kvol < vol) 
      {
         vol = kvol;
         sm = k;
      }
   }

   kvol = (ps.points[ps.nPoints - 1].objectives[0] - ps.points[ps.nPoints - 2].objectives[0]) * ps.points[ps.nPoints - 1].objectives[1];

   if (kvol < vol)
   {
      vol = kvol;
      sm = ps.nPoints - 1;
   }

   min[0] = ps.points[sm].objectives[0];
   min[1] = ps.points[sm].objectives[1];
   min[2] = vol;

   return ps.points[sm].idx;
}


int ihv(EMO_IWFG *hv, FRONT ps, double *min)
// returns the minimum exclusive hypervolume of points in ps
{
   int i, j, k, z, index;
   double width;
   int min_idx;

   for (i=0; i<MAX(ps.nPoints,hv->n); i++) 
      for (j=0; j<hv->n; j++) 
         hv->torder[i][j] = j;

   runHeuristic(hv, ps);

   for (i=0; i<ps.nPoints; i++)
   {
      hv->stacks[i][0].front = hv->fsorted[hv->torder[i][hv->n-1]];
      hv->stacks[i][0].width = 1;
      hv->stacksize[i] = 2;
   }

   sliceOrder(hv, ps.nPoints);
   hv->n--;

   for (i=0; i<ps.nPoints; i++)
   {
      SLICE top = hv->stacks[i][hv->stacksize[i]-1];

      while (hv->stacksize[i] < hv->maxStackSize && top.front.nPoints > SLICELIMIT)
      {
         hv->stacks[i][hv->stacksize[i]].front.nPoints = 0;
         index = 0;

         while (index < top.front.nPoints && ps.points[i].objectives[hv->torder[i][hv->n-1]] < top.front.points[index].objectives[hv->torder[i][hv->n-1]])
         {
            insert(hv, top.front.points[index], hv->n-2, hv->stacks[i][hv->stacksize[i]].front, i, hv->stacksize[i], hv->torder[i]);
            index++;
         }

         if (index < top.front.nPoints)
            hv->stacks[i][hv->stacksize[i]].width = ps.points[i].objectives[hv->torder[i][hv->n-1]]-top.front.points[index].objectives[hv->torder[i][hv->n-1]];
         else
            hv->stacks[i][hv->stacksize[i]].width = ps.points[i].objectives[hv->torder[i][hv->n-1]];

         hv->stacks[i][hv->stacksize[i]].index = index;
         top = hv->stacks[i][hv->stacksize[i]];
         hv->stacksize[i]++;
         hv->n--;
      }

      width = 1;

      for (j=0; j<hv->stacksize[i]; j++) 
         width *= hv->stacks[i][j].width;

      if (top.front.nPoints == 0) 
         hv->partial[i] = width * inclhvOrder(hv, ps.points[i],hv->torder[i]);
      else 
         hv->partial[i] = width * iwfg_exclhvPoint(hv, top.front,ps.points[i],hv->torder[i]);

      hv->n += hv->stacksize[i]-2;

      while (hv->stacksize[i]>1 && (top.index==hv->stacks[i][hv->stacksize[i]-2].front.nPoints || 
      iwfg_dominates1wayOrder(hv->stacks[i][hv->stacksize[i]-2].front.points[top.index],ps.points[i],hv->n - hv->stacksize[i]+1,hv->torder[i])))
      {
         hv->stacksize[i]--;
         top = hv->stacks[i][hv->stacksize[i]-1];
      }
   }

   initialiseHeap(hv, ps.nPoints);

   while (true)
   {
      // int i = removeFromHeap();
      i = peekFromHeap(hv);

      if (hv->stacksize[i]<=1)
      {
         for (z=0; z<ps.n; z++)
            min[z] = ps.points[i].objectives[z];

         min[ps.n] = hv->partial[i];
         min_idx = ps.points[i].idx;
         break;
      }

      hv->n -= hv->stacksize[i]-2;
      j = hv->stacks[i][hv->stacksize[i]-1].index;

      if (j<hv->stacks[i][hv->stacksize[i]-2].front.nPoints-1)
         hv->stacks[i][hv->stacksize[i]-1].width = hv->stacks[i][hv->stacksize[i]-2].front.points[j].objectives[hv->torder[i][hv->n]] -
                                                   hv->stacks[i][hv->stacksize[i]-2].front.points[j+1].objectives[hv->torder[i][hv->n]];
      else
         hv->stacks[i][hv->stacksize[i]-1].width = hv->stacks[i][hv->stacksize[i]-2].front.points[j].objectives[hv->torder[i][hv->n]];

      insert(hv, hv->stacks[i][hv->stacksize[i]-2].front.points[j], hv->n-1,hv->stacks[i][hv->stacksize[i]-1].front, i, hv->stacksize[i]-1,hv->torder[i]);
      hv->stacks[i][hv->stacksize[i]-1].index = j+1;
      SLICE top = hv->stacks[i][hv->stacksize[i]-1];

      width = 1;

      for (k=0; k<hv->stacksize[i]; k++) 
         width *= hv->stacks[i][k].width;

      if (top.front.nPoints == 0) 
         hv->partial[i] += width * inclhvOrder(hv, ps.points[i],hv->torder[i]);
      else 
         hv->partial[i] += width * iwfg_exclhvPoint(hv, top.front,ps.points[i],hv->torder[i]);

      hv->n += hv->stacksize[i]-2;

      while (hv->stacksize[i]>1 && (top.index==hv->stacks[i][hv->stacksize[i]-2].front.nPoints || 
      iwfg_dominates1wayOrder(hv->stacks[i][hv->stacksize[i]-2].front.points[top.index],ps.points[i],hv->n - hv->stacksize[i]+1,hv->torder[i])))
      {
         hv->stacksize[i]--;
         top = hv->stacks[i][hv->stacksize[i]-1];
      }

      heapify(hv, 0, i);
   }

   hv->n++;
   return min_idx;
}


//Allocate Memory, maxm = rows, maxn = objectives
void EMO_IWFG_alloc(EMO_IWFG *hv, int max_row, int col) {
   int i, j, maxdepth = col - 2; 

   hv->fs = (FRONT *) malloc(sizeof(FRONT)*maxdepth); 
   hv->fr = 0;  // RHG
   hv->n = col; // RHG
   hv->maxm = max_row;
   hv->maxn = col;

   printf("IWFG Algorithm\n");

   for(i = 0; i < maxdepth; i++) {
      hv->fs[i].points = (POINT *) malloc(sizeof(POINT)*max_row); 
      
      for(j=0; j<max_row; j++)
        hv->fs[i].points[j].objectives = (OBJECTIVE *) malloc(sizeof(OBJECTIVE)*(col - i - 1));  
   }

   hv->f.points = (POINT *) malloc(sizeof(POINT) * max_row);
   hv->f.nPoints = max_row;
   hv->f.n = col;
    
   for(i = 0; i < max_row; i++)
     hv->f.points[i].objectives = (OBJECTIVE *) malloc(sizeof(OBJECTIVE) * col);

   hv->partial = malloc(sizeof(double)*max_row);
   hv->heap = malloc(sizeof(int)*max_row);
   hv->stacksize = malloc(sizeof(int)*max_row); 
   hv->stacks = malloc(sizeof(SLICE*)*max_row);  //RHG stacks is a matrix of SLICE of size max_row x maxStackSize

   hv->maxStackSize = MIN(slicingDepth(hv->maxn),maxdepth)+1; 

   for(i=0; i<max_row; i++)
   {
      hv->stacks[i] = malloc(sizeof(SLICE)*hv->maxStackSize);
      
      for(j=1; j<hv->maxStackSize; j++)
         hv->stacks[i][j].front.points = malloc(sizeof(POINT)*max_row);
   }

   hv->fsorted = malloc(sizeof(FRONT)*col);

   for(i=0; i<col; i++)
      hv->fsorted[i].points = malloc(sizeof(POINT)*max_row);

   hv->torder = malloc(sizeof(int*)*MAX(max_row,col));
   hv->tcompare = malloc(sizeof(int*)*max_row);

   for(i=0; i<MAX(max_row,col); i++)
      hv->torder[i] = malloc(sizeof(int)*col);

   for(i=0; i<max_row; i++)
      hv->tcompare[i] = malloc(sizeof(int)*col); 
}

void EMO_IWFG_free(EMO_IWFG *hv)
{
  int i, j, maxdepth = hv->maxn-2;

  for(i = 0; i < hv->maxm; i++)
    free(hv->f.points[i].objectives);

  free(hv->f.points);

   for(i=0; i<maxdepth; i++)
   {
      for(j=0; j<hv->maxm; j++)
         free(hv->fs[i].points[j].objectives);

      free(hv->fs[i].points);
   }
   free(hv->fs);

   free(hv->partial);
   free(hv->heap);
   free(hv->stacksize);

   for(i=0; i<hv->maxm; i++) {
      for(j=1; j<hv->maxStackSize; j++) 
         free(hv->stacks[i][j].front.points);
      free(hv->stacks[i]);
   }

   free(hv->stacks);

   for(i=0; i<hv->maxn; i++)
      free(hv->fsorted[i].points);

   free(hv->fsorted);

   for(i=0; i<MAX(hv->maxm,hv->maxn); i++)
      free(hv->torder[i]);

   free(hv->torder);

   for(i=0; i<hv->maxm; i++)
      free(hv->tcompare[i]);

   free(hv->tcompare);
}

// Identify the smallest hypervolume-contributor of a set of non-dominated points.
int EMO_IWFG_run(EMO_IWFG *hv, double *data, int *enable, int row, const double *ref, double *hv_worst) {

  int i, j, k, c, idx;
  double eh[hv->n + 1];  // number of objectives + hypervolume
  double v;

  if(row > hv->maxm) {
    printf("Error, mismatch dimensions of row in EMO_HV_run (%d vs %d).\n", hv->maxm, row);
    exit(1);
  }

  c = 0;

  for(k = i = 0; i < row; i++) {
    if(enable == NULL || enable[i]) {
      for(j = 0; j < hv->maxn; j++) {
        // modify the objective values relative to the reference point 
        if(data[i* hv->maxn + j] >= ref[j]) {
          c++;
          k--;
          break;
        }
        else {
          v = fabs(data[i* hv->maxn + j] - ref[j]);
          hv->f.points[k].objectives[j] = v;
          hv->f.points[k].idx = i;
        }
      }
      k++;
    }
  }

  if(c > 0)
    printf("Warning: omitting %d/%d points that are outside the reference point (total %d)\n", c, row, k);

  hv->f.nPoints = k;

  if(hv->maxn == 2)
   idx = ihv2(hv, hv->f, eh);
  else 
   idx = ihv(hv, hv->f, eh);

  *hv_worst = eh[hv->n];

  return idx;
}

