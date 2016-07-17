/*
  Copyright (C) 2012 Alexander (Polyakov) Peletskyi

  This file is part of FFLCU.

  FFLCU is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  FFLCU is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

/*
 * fflComputeBH.cu
 *
 *  Created on: Jan 31, 2012
 *      Author: alexander
 */

#include "fflComputeBH.cuh"

#define SQ(x) ((x)*(x))
#define IND (blockDim.x * blockIdx.x + threadIdx.x)
#define DEBUG
const float myu0 = 4e-7 * M_PI;
const float g = 9.81;

/*__device__ template <typename T> int sgn(T val)
{
    return (val > T(0)) - (val < T(0));
}
*/
__device__ float sgn(float x){
	return (x > 0.0f) ? 1.0f : ((x < 0.0f) ? -1.0f : 0.0f);
}


using namespace std;

__constant__ int nnodesd, nbodiesd;
__constant__ volatile float epssqd, itolsqd;

__device__ volatile int bottomd, maxdepthd, blkcntd;
__device__ volatile float radiusd;

__constant__ ConstParams devCParBH;
__constant__ ChangableParams devChParBH;

__constant__ volatile float* xd;
__constant__ volatile float* yd;
__constant__ volatile float* zd;
__constant__ volatile float* phyd;
__constant__ volatile float* thetad;
__constant__ volatile float* uxd;
__constant__ volatile float* uyd;
__constant__ volatile float* uzd;

__constant__ volatile float* phy1d; //=d(phy)/(dt)
__constant__ volatile float* theta1d; //=d(theta)/(dt)

__constant__ volatile float* deltaPhyd;
__constant__ volatile float* deltaThetad;

__constant__ volatile float* x1d;
__constant__ volatile float* y1d;
__constant__ volatile float* z1d;

__constant__ volatile float* deltaXd;
__constant__ volatile float* deltaYd;
__constant__ volatile float* deltaZd;

__constant__ volatile float* massd;

__constant__ volatile int *errd;
__constant__ volatile int *sortd;
__constant__ volatile int *childd;
__constant__ volatile int *countd;
__constant__ volatile int *startd;

__constant__ volatile float *minxd;
__constant__ volatile float *maxxd;
__constant__ volatile float *minyd;
__constant__ volatile float *maxyd;
__constant__ volatile float *minzd;
__constant__ volatile float *maxzd;

__constant__ curandState *rndStatesd;

__global__ __launch_bounds__(THREADS0, FACTOR0)
void averageUKernel(float3* u)
{
	register int i, inc, cntr;
	float sinTheta, sinPhy, cosTheta, cosPhy;
	__shared__ float ux[THREADS0], uy[THREADS0], uz[THREADS0];
	__shared__ int shCntr;

	if (threadIdx.x + blockIdx.x * blockDim.x == 0) {
		u->x = 0.0f;
		u->y = 0.0f;
		u->z = 0.0f;
	}

	if (threadIdx.x == 0) {
		shCntr = 0;
	}

	ux[threadIdx.x] = 0.0f;
	uy[threadIdx.x] = 0.0f;
	uz[threadIdx.x] = 0.0f;

	// iterate over all nodes assigned to thread
	inc = blockDim.x * gridDim.x;
	cntr = 0;
	for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc) {
		sincos(thetad[i], &sinTheta, &cosTheta);
		sincos(phyd[i], &sinPhy, &cosPhy);
		ux[threadIdx.x] += sinTheta * cosPhy;
		uy[threadIdx.x] += sinTheta * sinPhy;
		uz[threadIdx.x] += cosTheta;
		cntr++;
	}

	atomicAdd(&shCntr, cntr);

	if (cntr > 0) {
		ux[threadIdx.x] = ux[threadIdx.x] / ((float)cntr);
		uy[threadIdx.x] = uy[threadIdx.x] / ((float)cntr);
		uz[threadIdx.x] = uz[threadIdx.x] / ((float)cntr);
	}

	__syncthreads();
	//cntr = 0;
	for (cntr = THREADS0 / 2; cntr > 0; cntr /= 2 ) {
		if (threadIdx.x < cntr) {
			if(threadIdx.x + blockIdx.x * blockDim.x + cntr < nbodiesd){
				ux[threadIdx.x] = 0.5 * (ux[threadIdx.x] + ux[threadIdx.x + cntr]);
				uy[threadIdx.x] = 0.5 * (uy[threadIdx.x] + uy[threadIdx.x + cntr]);
				uz[threadIdx.x] = 0.5 * (uz[threadIdx.x] + uz[threadIdx.x + cntr]);
			}
		}
		__syncthreads();
	}

	if (threadIdx.x == 0) {
		if (shCntr > 0){
			ux[0] *= (float)shCntr / (float)nbodiesd;
			uy[0] *= (float)shCntr / (float)nbodiesd;
			uz[0] *= (float)shCntr / (float)nbodiesd;
			atomicAdd(&(u->x), ux[0]);
			atomicAdd(&(u->y), uy[0]);
			atomicAdd(&(u->z), uz[0]);
		}
	}
}

/******************************************************************************/
/*** initialize memory ********************************************************/
/******************************************************************************/

__global__ void initializationKernel(curandState *state, unsigned long seed)
{
	register int ind;
	ind = IND;
	if (ind == 0) {
		*errd = 0;
		maxdepthd = 1;
		blkcntd = 0;
	}

	curand_init (seed, ind, 0, &state[ind]);
}



__global__
__launch_bounds__(THREADS0, FACTOR0)
void calcUKernel()
{
	register int i, inc;
	float sinTheta, sinPhy, cosTheta, cosPhy;
	// iterate over all nodes assigned to thread
	inc = blockDim.x * gridDim.x;
	for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc) {
		sincos(thetad[i], &sinTheta, &cosTheta);
		sincos(phyd[i], &sinPhy, &cosPhy);
		uxd[i] = sinTheta * cosPhy;
		uyd[i] = sinTheta * sinPhy;
		uzd[i] = cosTheta;
	}
}

/******************************************************************************/
/*** compute center and radius ************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS1, FACTOR1)
void boundingBoxKernel()
{
	register int i, j, k, inc;
	register float val, minx, maxx, miny, maxy, minz, maxz;
	__shared__ volatile float sminx[THREADS1], smaxx[THREADS1], sminy[THREADS1], smaxy[THREADS1], sminz[THREADS1], smaxz[THREADS1];

	// initialize with valid data (in case #bodies < #threads)
	minx = maxx = xd[0];
	miny = maxy = yd[0];
	minz = maxz = zd[0];

	// scan all bodies
	i = threadIdx.x;
	inc = THREADS1 * gridDim.x;
	for (j = i + blockIdx.x * THREADS1; j < nbodiesd; j += inc) {
		val = xd[j];
		minx = min(minx, val);
		maxx = max(maxx, val);
		val = yd[j];
		miny = min(miny, val);
		maxy = max(maxy, val);
		val = zd[j];
		minz = min(minz, val);
		maxz = max(maxz, val);
	}

	// reduction in shared memory
	sminx[i] = minx;
	smaxx[i] = maxx;
	sminy[i] = miny;
	smaxy[i] = maxy;
	sminz[i] = minz;
	smaxz[i] = maxz;

	for (j = THREADS1 / 2; j > 0; j /= 2) {
		__syncthreads();
		if (i < j) {
			k = i + j;
			sminx[i] = minx = min(minx, sminx[k]);
			smaxx[i] = maxx = max(maxx, smaxx[k]);
			sminy[i] = miny = min(miny, sminy[k]);
			smaxy[i] = maxy = max(maxy, smaxy[k]);
			sminz[i] = minz = min(minz, sminz[k]);
			smaxz[i] = maxz = max(maxz, smaxz[k]);
		}
	}

	// write block result to global memory
	if (i == 0) {
		k = blockIdx.x;
		minxd[k] = minx;
		maxxd[k] = maxx;
		minyd[k] = miny;
		maxyd[k] = maxy;
		minzd[k] = minz;
		maxzd[k] = maxz;

		inc = gridDim.x - 1;
		if (inc == atomicInc((unsigned int *)&blkcntd, inc)) {
		// I'm the last block, so combine all block results
			for (j = 0; j <= inc; j++) {
				minx = min(minx, minxd[j]);
				maxx = max(maxx, maxxd[j]);
				miny = min(miny, minyd[j]);
				maxy = max(maxy, maxyd[j]);
				minz = min(minz, minzd[j]);
				maxz = max(maxz, maxzd[j]);
			}

			// compute 'radius'
			val = max(maxx - minx, maxy - miny);
			radiusd = max(val, maxz - minz) * 0.5f;

			// create root node
			k = nnodesd;
			bottomd = k;

			massd[k] = -1.0f;
			startd[k] = 0;
			xd[k] = (minx + maxx) * 0.5f;
			yd[k] = (miny + maxy) * 0.5f;
			zd[k] = (minz + maxz) * 0.5f;
			k *= 8;
			for (i = 0; i < 8; i++) childd[k + i] = -1;
		}
	}
}


/******************************************************************************/
/*** build tree ***************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS2, FACTOR2)
void treeBuildingKernel()
{
	register int i, j, k, depth, localmaxdepth, skip, inc;
	register float x, y, z, r;
	register float px, py, pz;
	register int ch, n, cell, locked, patch;
	register float radius, rootx, rooty, rootz;

	// cache root data
	radius = radiusd;
	rootx = xd[nnodesd];
	rooty = yd[nnodesd];
	rootz = zd[nnodesd];

	localmaxdepth = 1;
	skip = 1;
	inc = blockDim.x * gridDim.x;
	i = threadIdx.x + blockIdx.x * blockDim.x;

	// iterate over all bodies assigned to thread
	while (i < nbodiesd) {
		if (skip != 0) {
			// new body, so start traversing at root
			skip = 0;
			px = xd[i];
			py = yd[i];
			pz = zd[i];
			n = nnodesd;
			depth = 1;
			r = radius;
			j = 0;
			// determine which child to follow
			if (rootx < px) j = 1;
			if (rooty < py) j += 2;
			if (rootz < pz) j += 4;
		}

		// follow path to leaf cell
		ch = childd[n * 8 + j];
		while (ch >= nbodiesd) {
			n = ch;
			depth++;
			r *= 0.5f;
			j = 0;
			// determine which child to follow
			if (xd[n] < px) j = 1;
			if (yd[n] < py) j += 2;
			if (zd[n] < pz) j += 4;
			ch = childd[n * 8 + j];
		}

		if (ch != -2) {	// skip if child pointer is locked and try again later
			locked = n  * 8 + j;
			if (ch == atomicCAS((int *)&childd[locked], ch, -2)) {	// try to lock
				if (ch == -2) {
					printf("!!!!!!!Error: ch = -2\n");
					break;
				}
				if (ch == -1) {
					// if null, just insert the new body
					childd[locked] = i;
				} else {	// there already is a body in this position
					patch = -1;
					// create new cell(s) and insert the old and new body
					do {
						depth++;

						cell = atomicSub((int *)&bottomd, 1) - 1;
						if (cell <= nbodiesd) {
							*errd = 1;
							bottomd = nnodesd;
						}
						patch = max(patch, cell);

						x = (j & 1) * r;
						y = ((j >> 1) & 1) * r;
						z = ((j >> 2) & 1) * r;
						r *= 0.5f;

						massd[cell] = -1.0f;
						startd[cell] = -1;
						x = xd[cell] = xd[n] - r + x;
						y = yd[cell] = yd[n] - r + y;
						z = zd[cell] = zd[n] - r + z;
						for (k = 0; k < 8; k++)
							childd[cell* 8 + k] = -1;

						if (patch != cell) {
							childd[n * 8 + j] = cell;
						}

						j = 0;
						if (x < xd[ch]) j = 1;
						if (y < yd[ch]) j += 2;
						if (z < zd[ch]) j += 4;
						childd[cell * 8 + j] = ch;

						n = cell;
						j = 0;
						if (x < px) j = 1;
						if (y < py) j += 2;
						if (z < pz) j += 4;

						ch = childd[n * 8 + j];
						// repeat until the two bodies are different children
					} while (ch >= 0);

					childd[n * 8 + j] = i;
					__threadfence();	// push out subtree
					childd[locked] = patch;
				}

				localmaxdepth = max(depth, localmaxdepth);
				i += inc;	// move on to next body
				skip = 1;
			}
		}
		__syncthreads();	// throttle
	}
	// record maximum tree depth
	atomicMax((int *)&maxdepthd, localmaxdepth);
#ifdef DEBUG
//	if ( maxdepthd == localmaxdepth ) printf ("Max tree depth: %d\n", maxdepthd);
#endif
}


__global__ void __treeBuildingKernel()
{
  register int i, j, k, depth, localmaxdepth, skip, inc;
  register float x, y, z, r;
  register float px, py, pz;
  register int ch, n, cell, locked, patch;
  __shared__ float radius, rootx, rooty, rootz;

  i = threadIdx.x;
  if (i == 0) {
    // cache root data
    radius = radiusd;
    rootx = xd[nnodesd];
    rooty = yd[nnodesd];
    rootz = zd[nnodesd];
  }
  __syncthreads();

  localmaxdepth = 1;
  skip = 1;
  inc = blockDim.x * gridDim.x;
  i += blockIdx.x * blockDim.x;

  // iterate over all bodies assigned to thread
  while (i < nbodiesd) {
    if (skip != 0) {
      // new body, so start traversing at root
      skip = 0;
      px = xd[i];
      py = yd[i];
      pz = zd[i];
      n = nnodesd;
      depth = 1;
      r = radius;
      j = 0;
      // determine which child to follow
      if (rootx < px) j = 1;
      if (rooty < py) j += 2;
      if (rootz < pz) j += 4;
    }

    ch = childd[n*8+j];
    // follow path to leaf cell
    while (ch >= nbodiesd) {
      n = ch;
      depth++;
      r *= 0.5f;
      j = 0;
      // determine which child to follow
      if (xd[n] < px) j = 1;
      if (yd[n] < py) j += 2;
      if (zd[n] < pz) j += 4;
      ch = childd[n*8+j];
    }

    if (ch != -2) {  // skip if child pointer is locked and try again later
      locked = n*8+j;
      if (ch == atomicCAS((int*)&childd[locked], ch, -2)) {  // try to lock
        if (ch == -1) {
          // if null, just insert the new body
          childd[locked] = i;
        } else {  // there already is a body in this position
          patch = -1;
          // create new cell(s) and insert the old and new body
          do {
            depth++;

            cell = atomicSub((int*)&bottomd, 1) - 1;
            if (cell <= nbodiesd) {
              *errd = 1;
              bottomd = nnodesd;
            }
            patch = max(patch, cell);

            x = (j & 1) * r;
            y = ((j >> 1) & 1) * r;
            z = ((j >> 2) & 1) * r;
            r *= 0.5f;

            massd[cell] = -1.0f;
            startd[cell] = -1;
            x = xd[cell] = xd[n] - r + x;
            y = yd[cell] = yd[n] - r + y;
            z = zd[cell] = zd[n] - r + z;
#pragma unroll 8
            for (k = 0; k < 8; k++) childd[cell*8+k] = -1;

            if (patch != cell) {
              childd[n*8+j] = cell;
            }

            j = 0;
            if (x < xd[ch]) j = 1;
            if (y < yd[ch]) j += 2;
            if (z < zd[ch]) j += 4;
            childd[cell*8+j] = ch;

            n = cell;
            j = 0;
            if (x < px) j = 1;
            if (y < py) j += 2;
            if (z < pz) j += 4;

            ch = childd[n*8+j];
            // repeat until the two bodies are different children
          } while (ch >= 0);
          childd[n*8+j] = i;
          __threadfence();
          childd[locked] = patch;
        }

        localmaxdepth = max(depth, localmaxdepth);
        i += inc;  // move on to next body
        skip = 1;
      }
    }
    __syncthreads();
  }
  atomicMax((int*)&maxdepthd, localmaxdepth);
}


/******************************************************************************/
/*** compute center of mass ***************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS3, FACTOR3)
void summarizationKernel()
{
	register int i, j, k, ch, inc, missing, cnt, bottom;
	register float m, cm, px, py, pz, ux, uy, uz;
	__shared__ volatile int child[THREADS3 * 8];

	bottom = bottomd;
	inc = blockDim.x * gridDim.x;
	k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;	// align to warp size
	if (k < bottom) k += inc;

	missing = 0;
	//int iteration = 0;
	// iterate over all cells assigned to thread
	while (k <= nnodesd) {
		//iteration++;
		if (missing == 0) {
			// new cell, so initialize
			cm = 0.0f;
			px = 0.0f;
			py = 0.0f;
			pz = 0.0f;
			ux = 0.0f;
			uy = 0.0f;
			uz = 0.0f;
			cnt = 0;
			j = 0;
			for (i = 0; i < 8; i++) {
				ch = childd[k * 8 + i];
				if (ch >= 0) {
			if (i != j) {
				// move children to front (needed later for speed)
				childd[k* 8 + i] = -1;
				childd[k* 8 + j] = ch;
			}
			child[missing * THREADS3 + threadIdx.x] = ch;	// cache missing children
			m = massd[ch];
			missing++;
			if (m >= 0.0f) {
				// child is ready
				missing--;
				if (ch >= nbodiesd) {	// count bodies (needed later)
					cnt += countd[ch] - 1;
				}
				// add child's contribution
				cm += m;
				px += xd[ch] * m;
				py += yd[ch] * m;
				pz += zd[ch] * m;
				ux += uxd[ch];
				uy += uyd[ch];
				uz += uzd[ch];
			}
			j++;
				}
			}
			cnt += j;
		}

		if (missing != 0) {
			do {
				// poll missing child
				ch = child[(missing - 1) * THREADS3 + threadIdx.x];
				m = massd[ch];
				if (m >= 0.0f) {
					// child is now ready
					missing--;
					if (ch >= nbodiesd) {
						// count bodies (needed later)
						cnt += countd[ch] - 1;
					}
					// add child's contribution
					cm += m;
					px += xd[ch] * m;
					py += yd[ch] * m;
					pz += zd[ch] * m;
					ux += uxd[ch];
					uy += uyd[ch];
					uz += uzd[ch];
				}
				// repeat until we are done or child is not ready
			} while ((m >= 0.0f) && (missing != 0));
		}

		if (missing == 0) {
			// all children are ready, so store computed information
	//		if (m > 1E-5){
				countd[k] = cnt;
				m = 1.0f / cm;
				xd[k] = px * m;
				yd[k] = py * m;
				zd[k] = pz * m;

				uxd[k] = ux;
				uyd[k] = uy;
				uzd[k] = uz;
	/*		} else {
				printf ("Warning - zero cell in summation");
				countd[k] = cnt;
				m = 0.0f;
				xd[k] = px * m;
				yd[k] = py * m;
				zd[k] = pz * m;

				uxd[k] = 0.0f;
				uyd[k] = 0.0f;
				uzd[k] = 0.0f;
			}*/
			__threadfence();	// make sure data are visible before setting mass
			massd[k] = cm;
			k += inc;	// move on to next cell
		}
	/*	if (iteration > 100000){
			printf ("To much iterations");
			printf ("k = %d, nnodes = %d", k, nnodesd);
			break;
		}*/
	}	//while
}


__global__ void __summarizationKernel()
{
  register int i, j, k, ch, inc, missing, cnt;
  register float m, cm, px, py, pz, ux, uy, uz;
  __shared__ int bottom, child[THREADS3 * 8];

  if (0 == threadIdx.x) {
    bottom = bottomd;
  }
  __syncthreads();

  inc = blockDim.x * gridDim.x;
  k = (bottom & (-WARPSIZE)) + threadIdx.x + blockIdx.x * blockDim.x;  // align to warp size
  if (k < bottom) k += inc;

  missing = 0;
  // iterate over all cells assigned to thread
  while (k <= nnodesd) {
    if (missing == 0) {
      // new cell, so initialize
		cm = 0.0f;
		px = 0.0f;
		py = 0.0f;
		pz = 0.0f;
		ux = 0.0f;
		uy = 0.0f;
		uz = 0.0f;
		cnt = 0;
		j = 0;
#pragma unroll 8
      for (i = 0; i < 8; i++) {
        ch = childd[k*8+i];
        if (ch >= 0) {
          if (i != j) {
            // move children to front (needed later for speed)
            childd[k*8+i] = -1;
            childd[k*8+j] = ch;
          }
          child[missing*THREADS3+threadIdx.x] = ch;  // cache missing children
          m = massd[ch];
          missing++;
          if (m >= 0.0f) {
            // child is ready
            missing--;
            if (ch >= nbodiesd) {  // count bodies (needed later)
              cnt += countd[ch] - 1;
            }
            // add child's contribution
			cm += m;
			px += xd[ch] * m;
			py += yd[ch] * m;
			pz += zd[ch] * m;
			ux += uxd[ch];
			uy += uyd[ch];
			uz += uzd[ch];;
          }
          j++;
        }
      }
      cnt += j;
    }

    if (missing != 0) {
      do {
        // poll missing child
        ch = child[(missing-1)*THREADS3+threadIdx.x];
        m = massd[ch];
        if (m >= 0.0f) {
          // child is now ready
          missing--;
          if (ch >= nbodiesd) {
            // count bodies (needed later)
            cnt += countd[ch] - 1;
          }
          // add child's contribution
			cm += m;
			px += xd[ch] * m;
			py += yd[ch] * m;
			pz += zd[ch] * m;
			ux += uxd[ch];
			uy += uyd[ch];
			uz += uzd[ch];
        }
        // repeat until we are done or child is not ready
      } while ((m >= 0.0f) && (missing != 0));
    }

    if (missing == 0) {
      // all children are ready, so store computed information
		countd[k] = cnt;
		m = 1.0f / cm;
		xd[k] = px * m;
		yd[k] = py * m;
		zd[k] = pz * m;

		uxd[k] = ux;
		uyd[k] = uy;
		uzd[k] = uz;
      __threadfence();
      massd[k] = cm;
      k += inc;  // move on to next cell
    }
  }
}


/******************************************************************************/
/*** sort bodies **************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS4, FACTOR4)
void sortKernel()
{
	register int i, k, ch, dec, start, bottom;

	bottom = bottomd;
	dec = blockDim.x * gridDim.x;
	k = nnodesd + 1 - dec + threadIdx.x + blockIdx.x * blockDim.x;

	// iterate over all cells assigned to thread
	while (k >= bottom) {
		start = startd[k];
		if (start >= 0) {
			for (i = 0; i < 8; i++) {
				ch = childd[k* 8 + i];
				if (ch >= nbodiesd) {
			// child is a cell
			startd[ch] = start;	// set start ID of child
			start += countd[ch];	// add #bodies in subtree
				} else if (ch >= 0) {
			// child is a body
			sortd[start] = ch;	// record body in 'sorted' array
			start++;
				}
			}
			k -= dec;	// move on to next cell
		}
		__syncthreads();	// throttle
	}
}


/******************************************************************************/
/*** compute force ************************************************************/
/******************************************************************************/

__global__
__launch_bounds__(THREADS5, FACTOR5)
void forceCalculationKernel()
{
	register int i, j, k, n, depth, base, sbase, diff, t;
	register float px, py, pz, dx, dy, dz, tmp, fx, fy, fz, hx, hy, hz, ux, uy, uz;
	register float ucx, ucy, ucz;
	__shared__ volatile int pos[MAXDEPTH * THREADS5/WARPSIZE], node[MAXDEPTH * THREADS5/WARPSIZE];
	__shared__ float dq[MAXDEPTH * THREADS5/WARPSIZE];
	curandState localState;
	//localState = rndStatesd[IND];

	float b, b2, d1, dd5;
	float bb2d7, umd5;
	float pow3s2d2, flj;

	if (0 == threadIdx.x) {
		tmp = radiusd;
		// precompute values that depend only on tree level
		dq[0] = tmp * tmp * itolsqd;
		for (i = 1; i < maxdepthd; i++) {
			dq[i] = dq[i - 1] * 0.25f;
			dq[i - 1] += epssqd;
		}
		dq[i - 1] += epssqd;

		if (maxdepthd > MAXDEPTH) {
			*errd = maxdepthd;
		}
	}
	__syncthreads();

	if (maxdepthd <= MAXDEPTH) {
		// figure out first thread in each warp (lane 0)
		base = threadIdx.x / WARPSIZE;
		sbase = base * WARPSIZE;
		j = base * MAXDEPTH;

		diff = threadIdx.x - sbase;
		// make multiple copies to avoid index calculations later
		if (diff < MAXDEPTH) {
			dq[diff+j] = dq[diff];
		}
		__syncthreads();

		// iterate over all bodies assigned to thread
		for (k = threadIdx.x + blockIdx.x * blockDim.x; k < nbodiesd; k += blockDim.x * gridDim.x) {
			i = sortd[k];	// get permuted/sorted index
			// cache position info
			px = xd[i];
			py = yd[i];
			pz = zd[i];

			ux = uxd[i];
			uy = uyd[i];
			uz = uzd[i];

			hx = 0.0f;
			hy = 0.0f;
			hz = 0.0f;

			fx = 0.0f;
			fy = 0.0f;
			fz = 0.0f;

			// initialize iteration stack, i.e., push root node onto stack
			depth = j;
			if (sbase == threadIdx.x) {
				node[j] = nnodesd;
				pos[j] = 0;
			}

			while (depth >= j) {
				// stack is not empty
				while ((t = pos[depth]) < 8) {
					// node on top of stack has more children to process
					n = childd[node[depth] * 8 + t];	// load child pointer
					if (sbase == threadIdx.x) {
						// I'm the first thread in the warp
						pos[depth] = t + 1;
					}
					if (n >= 0) {
						dx = -(xd[n] - px);
						dy = -(yd[n] - py);
						dz = -(zd[n] - pz);
						tmp = dx * dx + (dy * dy + dz * dz);	// compute distance squared
						if ((n < nbodiesd) || __all(tmp >= dq[depth])) {	// check if all threads agree that cell is far enough away (or is a body)
							if (n != i) {

								ucx = uxd[n];
								ucy = uyd[n];
								ucz = uzd[n];

								b = ucx * dx + ucy * dy + ucz * dz;
								b2 = ux * dx + uy * dy + uz * dz;

								d1 = sqrtf(tmp/*, 0.5f*/);
								dd5 = __fdividef(1.0f, tmp * tmp * d1);
								bb2d7 = 15.0f * b * b2 * __fdividef(dd5, tmp);
								umd5 = - 3.0f * (ux*ucx + uy*ucy + uz*ucz) * dd5;

								hx += (b * 3.0f * dx - tmp * ucx) * dd5;
								hy += (b * 3.0f * dy - tmp * ucy) * dd5;
								hz += (b * 3.0f * dz - tmp * ucz) * dd5;


								fx += -dx * (umd5 + bb2d7)
								+ 3.0f * (b * ux + b2 * ucx) * dd5;

								fy += -dy * (umd5  +  bb2d7)
								+ 3.0f * (b * uy + b2 * ucy) * dd5;

								fz += -dz * (umd5  +  bb2d7)
								+ 3.0f * (b * uz + b2 * ucz) * dd5;

								/*if (fx != fx || fy != fy || fz != fz) {	//nan
									printf("NAN in particle iteraction (Before LJ)[%d] and meta:[%d]\n", i, n);
									printf("x = %f, y = %f, z = %f,\n", px, py, pz);
									printf("x = %f, y = %f, z = %f,\n", xd[n], xd[n], xd[n]);
									printf("fx = %f, fy = %f, fz = %f,\n", fx, fy, fz);
									printf("hx = %f, hy = %f, hz = %f,\n", hx, hy, hz);
									if (fx != fx) fx = 0.0f;
									if (fy != fy) fy = 0.0f;
									if (fz != fz) fz = 0.0f;
								}*/


								/* Lennard Jonnes force is equal:
								 * flj = 24.0 * E * (2 * __powf(SIGMA2/D2, 6) - __powf(SIGMA2/D2, 3)) / d1;
								 * where SIGMA2 is SIGMA*SIGMA;
								 * But each component is equal to (for x) fx = dx  * flj / D1;
								 * the last D1 in were included to force flj to decrease calculations;
								 * */

								if (d1 < devCParBH.ljCutOffR * devCParBH.sigma) {
									pow3s2d2 = __powf(SQ(devCParBH.sigma)/tmp, 3);
									flj = 24.0f * devCParBH.eps * (pow3s2d2  * (2.0f * pow3s2d2 - 1.0f)) / tmp;
									flj -= devCParBH.ljCutOffForce;
									fx +=  dx * flj;
									fy +=  dy * flj;
									fz +=  dz * flj;

									/*if (fx != fx || fy != fy || fz != fz) {	//nan
										printf("NAN in particle iteraction (After LJ)[%d] and meta:[%d]. Flj = %f\n", i, n, flj);
										printf("x = %f, y = %f, z = %f,\n", px, py, pz);
										printf("x = %f, y = %f, z = %f,\n", xd[n], xd[n], xd[n]);
										printf("fx = %f, fy = %f, fz = %f,\n", fx, fy, fz);
										printf("hx = %f, hy = %f, hz = %f,\n", hx, hy, hz);
										if (fx != fx) fx = 0.0f;
										if (fy != fy) fy = 0.0f;
										if (fz != fz) fz = 0.0f;
									}*/
								}
							}
						} else {
							// push cell onto stack
							depth++;
							if (sbase == threadIdx.x) {
								node[depth] = n;
								pos[depth] = 0;
							}
						}
					} else {
						depth = max(j, depth - 1);	// early out because all remaining children are also zero
					}
				}
				depth--;	// done with this level
			}

		//	fz += 0.7f;

			if (devCParBH.cf == ConstParams::BARREL) {
				//making center (0;0)
				px = px - devCParBH.barrelR;
				py = py - devCParBH.barrelR;
				float xr, yr;
				//calculating nearest point of the circle to current particle
				if (fabs(px) < 1E-6) {
					xr = 0;
					if (py > 0) yr = devCParBH.barrelR;
					else yr = -devCParBH.barrelR;
				} else {
					xr = sqrtf(SQ(devCParBH.barrelR) / (1.0f + SQ(py/px))/*, 0.5f*/);

					//there are 2 roots, we should define + or -
					if (px < 0.0f) xr = -xr;

					yr = xr * py / px;
				}

				dx = px - xr;
				dy = py - yr;
				tmp = (dx * dx  + dy * dy );

#define walls 15.0f

				pow3s2d2 = __powf(SQ(devCParBH.sigmaWall) / tmp, 3);
				flj = 24.0f * devCParBH.epsWall * (pow3s2d2  * (2.0f * pow3s2d2 - 1.0f)) / tmp;
				fx +=  dx * flj;	//maybe should be -dx...
				fy +=  dy * flj;
				fz += 24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/pz, 12) - __powf(devCParBH.sigmaWall/pz, 6)) / pz;
				fz += -24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/(devCParBH.lz - pz), 12) - __powf(devCParBH.sigmaWall/(devCParBH.lz - pz), 6)) / (devCParBH.lz - pz);


				/*flj = walls / (expf(walls * sgn(-sqrtf(px * px + py * py) + devCParBH.barrelR) * sqrtf(fabs(tmp))) + 1.0f);
				fx +=  -(dx) * flj;
				fy +=  -(dy) * flj;

				flj = walls / (expf(walls * sgn(pz) * sqrtf(fabs(pz))) + 1.0f);
				if (flj != flj) printf("1FLJ is NAN: pz = %f\n", pz);
				fz += flj;
				flj = - walls / (expf(walls * sgn(devCParBH.lz - pz) * sqrtf(fabs(devCParBH.lz - pz))) + 1.0f);
				if (flj != flj) printf("2FLJ is NAN: pz = %f\n", pz);
				fz += flj;*/
			} else {

			/*	if (fx < 0) fx = fx / (expf(-walls * px) + 1.0f);
				else fx = fx / (expf(-walls * (devCParBH.lx - px)) + 1.0f);

				if (fy < 0) fy = fy / (expf(-walls * py) + 1.0f);
				else fy = fy / (expf(-walls * (devCParBH.ly - py)) + 1.0f);

				if (fz < 0) fz = fz / (expf(-walls * pz) + 1.0f);
				else fz = fz / (expf(-walls * (devCParBH.lz - pz)) + 1.0f);*/

				fx += 24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/px, 12) - __powf(devCParBH.sigmaWall/px, 6)) / px;
				fy += 24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/py, 12) - __powf(devCParBH.sigmaWall/py, 6)) / py;
				fz += 24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/pz, 12) - __powf(devCParBH.sigmaWall/pz, 6)) / pz;

				fx += -24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/(devCParBH.lx - px), 12) - __powf(devCParBH.sigmaWall/(devCParBH.lx - px), 6)) / (devCParBH.lx - px);
				fy += -24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/(devCParBH.ly - py), 12) - __powf(devCParBH.sigmaWall/(devCParBH.ly - py), 6)) / (devCParBH.ly - py);
				fz += -24.0f * devCParBH.epsWall * (2.0f * __powf(devCParBH.sigmaWall/(devCParBH.lz - pz), 12) - __powf(devCParBH.sigmaWall/(devCParBH.lz - pz), 6)) / (devCParBH.lz - pz);



		/*		flj = walls / (expf(walls * px) + 1.0f);
				if (flj != flj) printf("1FLJ is NAN: px = %f\n", px);
				fx += flj;
				flj = - walls / (expf(walls * (devCParBH.lx - px)) + 1.0f);
				if (flj != flj) printf("2FLJ is NAN: px = %f\n", px);
				fx += flj;

				flj = walls / (expf(walls * py) + 1.0f);
				if (flj != flj) printf("1FLJ is NAN: py = %f\n", py);
				fy += flj;
				flj = - walls / (expf(walls * (devCParBH.ly - py)) + 1.0f);
				if (flj != flj) printf("2FLJ is NAN: py = %f\n", py);
				fy += flj;

				flj = walls / (expf(walls * pz) + 1.0f);
				if (flj != flj) printf("1FLJ is NAN: pz = %f\n", pz);
				fz += flj;
				flj = - walls / (expf(walls * (devCParBH.lz - pz)) + 1.0f);
				if (flj != flj) printf("2FLJ is NAN: pz = %f\n", pz);
				fz += flj;*/

			}

			hx += devChParBH.hExtX;
			hy += devChParBH.hExtY;
			hz += devChParBH.hExtZ;

#define	NX (uy * hz - uz * hy)
#define	NY (uz * hx - ux * hz)
#define	NZ (ux * hy - uy * hx)

			deltaPhyd[i] = phy1d[i] * devChParBH.dTimeCurrent;
			deltaThetad[i] = theta1d[i] * devChParBH.dTimeCurrent;

			deltaXd[i] = x1d[i] * devChParBH.dTimeCurrent;
			deltaYd[i] = y1d[i] * devChParBH.dTimeCurrent;
			deltaZd[i] = z1d[i] * devChParBH.dTimeCurrent;
#define MAX_DELTA 0.1f
			if (fabs(deltaXd[i]) > MAX_DELTA) {
				printf("Too big dx[%d] = %f\n", i, deltaXd[i]);
//				deltaXd[i] = MAX_DELTA * deltaXd[i] /  fabs(deltaXd[i]);
				deltaXd[i] = 1E10;
			}
			if (fabs(deltaYd[i]) > MAX_DELTA) {
				printf("Too big dy[%d] = %f\n", i, deltaYd[i]);
//				deltaYd[i] = MAX_DELTA * deltaYd[i] /  fabs(deltaYd[i]);
				deltaYd[i] = 1E10;
			}
			if (fabs(deltaZd[i]) > MAX_DELTA) {
				printf("Too big dz[%d] = %f\n", i, deltaZd[i]);
//				deltaZd[i] = MAX_DELTA * deltaZd[i] /  fabs(deltaZd[i]);
				deltaZd[i] = 1E10;
			}

/*			if (devCParBH.gravitation == true) {
				fz += - devCParBH.r * g * (devCParBH.roParticles - devCParBH.roEnvironment)
						/ (4.0f * M_PI * devCParBH.myu * devCParBH.myu * myu0 / 3.0f);
			}*/

			if (fx != fx || fy != fy || fz != fz ||
				hx != hx || hy != hy || hz != hz) {	//nan
				printf("Force Kernel: NAN in particle[%d]\n", i);
				printf("x = %f, y = %f, z = %f,\n", px, py, pz);
				printf("fx = %f, fy = %f, fz = %f,\n", fx, fy, fz);
				printf("hx = %f, hy = %f, hz = %f,\n", hx, hy, hz);
/*				if (fx != fx) fx = 0.0f;
				if (fy != fy) fy = 0.0f;
				if (fz != fz) fz = 0.0f;
				if (hx != hx) hx = 0.0f;
				if (hy != hy) hy = 0.0f;
				if (hz != hz) hz = 0.0f;*/
				fx = 1E10;
				fy = 1E10;
				fz = 1E10;
				hx = 0.0f;
				hy = 0.0f;
				hz = 0.0f;
			}

			if (devCParBH.thermalBath == true) {
				localState = rndStatesd[i];

				phy1d[i] += 2.5f * ((NZ - phy1d[i] * devCParBH.nyu) * devChParBH.dTimeCurrent
						+ curand_normal(&localState) /* 3.0f * devCParBH.nyu*/ * devChParBH.sqrtdTime * devCParBH.qr);

				theta1d[i] += 2.5f * ((- NX * __sinf(phyd[i]) + NY * __cosf(phyd[i])
					- theta1d[i] * devCParBH.nyu) * devChParBH.dTimeCurrent
					+ curand_normal(&localState) /* 3.0f * devCParBH.nyu*/ * devChParBH.sqrtdTime * devCParBH.qr);

				x1d[i] += (fx - x1d[i] * devCParBH.eta) * devChParBH.dTimeCurrent
						+ curand_normal(&localState) /* 2.0f * devCParBH.eta*/ * devChParBH.sqrtdTime * devCParBH.qt;
				y1d[i] += (fy - y1d[i] * devCParBH.eta) * devChParBH.dTimeCurrent
					+ curand_normal(&localState) /* 2.0f * devCParBH.eta*/ * devChParBH.sqrtdTime * devCParBH.qt;
				z1d[i] += (fz - z1d[i] * devCParBH.eta) * devChParBH.dTimeCurrent
					+ curand_normal(&localState) /* 2.0f * devCParBH.eta*/ * devChParBH.sqrtdTime * devCParBH.qt;

				rndStatesd[i] = localState;


				//if(i % 100 == 0) printf("x1[%d] = %f, y1[%d] = %f, z1[%d] = %f\n", i, x1d[i], i, y1d[i], i, z1d[i]);
				//printf("%d rnd %f\n", i, curand_normal(&localState));
			} else {
				phy1d[i] += 2.5f * (NZ - phy1d[i] * devCParBH.nyu) * devChParBH.dTimeCurrent;
				theta1d[i] += 2.5f * (- NX * __sinf(phyd[i]) + NY * __cosf(phyd[i])
					- theta1d[i] * devCParBH.nyu) * devChParBH.dTimeCurrent;

				x1d[i] += (fx - x1d[i] * devCParBH.eta) * devChParBH.dTimeCurrent;
				y1d[i] += (fy - y1d[i] * devCParBH.eta) * devChParBH.dTimeCurrent;
				z1d[i] += (fz - z1d[i] * devCParBH.eta) * devChParBH.dTimeCurrent;
				//if(i % 100 == 0) printf("x1[%d] = %f, y1[%d] = %f, z1[%d] = %f\n", i, x1d[i], i, y1d[i], i, z1d[i]);

			}
		}
	}
	//rndStatesd[IND] = localState;
}


/******************************************************************************/
/*** advance bodies ***********************************************************/
/******************************************************************************/

__device__
bool checkLocation(float x, float y, float z){
	if (devCParBH.cf == ConstParams::BARREL) {
		if (SQ(x - devCParBH.barrelR) + SQ(y - devCParBH.barrelR) < SQ(devCParBH.barrelR)
				&& z > 0.0f && z < devCParBH.lz) return true;
		else return false;
	} else {
		if (x > 0.0f && x < devCParBH.lx &&
			y > 0.0f && y < devCParBH.ly &&
			z > 0.0f && z < devCParBH.lz) return true;
		else return false;
	}
}

__device__
bool checkLocationFull(float x, float y, float z, int n){
	float xi, yi, zi;
	if (checkLocation(x, y, z) == false) return false;

	for (int i = 0; i < nbodiesd; i++){
		xi = xd[i];
		yi = yd[i];
		zi = zd[i];

		if (fabs(x - xi) < 2.0f &&
			fabs(y - yi) < 2.0f &&
			fabs(z - zi) < 2.0f &&
			i != n) {
			if (SQ(x - xi) + SQ(y - yi) + SQ(z - zi) < 4.0f) return false;	//4 = 2R * 2R = 4R = 4
		}
	}
	return true;
}

__global__
__launch_bounds__(THREADS6, FACTOR6)
void integrationKernel() {
	register int i, inc;
	register float deltaPhy, deltaTheta, deltaX, deltaY, deltaZ;
	register float x, y, z;
	// iterate over all bodies assigned to thread
	inc = blockDim.x * gridDim.x;
	for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc) {
		// integrate
		deltaPhy = deltaPhyd[i];
		deltaTheta = deltaThetad[i];
		deltaX = deltaXd[i];
		deltaY = deltaYd[i];
		deltaZ = deltaZd[i];

		if(deltaPhy != deltaPhy || deltaTheta != deltaTheta ||
				deltaX != deltaX || deltaY != deltaY || deltaZ != deltaZ) {
			printf("Integration kernel: NAN in particle[%d]\n", i);
			printf("x = %f, y = %f, z = %f,\n", xd[i], yd[i], zd[i]);
			printf("deltaX = %f, deltaY = %f, deltaZ = %f,\n", deltaX, deltaY, deltaZ);
			printf("deltaPhy = %f, deltaTheta = %f\n", deltaPhy, deltaTheta);
			if (deltaPhy != deltaPhy) deltaPhy = 0;
			if (deltaTheta != deltaTheta) deltaTheta = 0;
			if (deltaX != deltaX) deltaX = 1E10;	//move particle fare from current place
			if (deltaY != deltaY) deltaY = 1E10; //then it will be automatically located in a random place
			if (deltaZ != deltaZ) deltaZ = 1E10;
		}

		phyd[i] += deltaPhy;
		thetad[i] += deltaTheta;
		/*xd[i] += deltaX;
		yd[i] += deltaY;
		zd[i] += deltaZ;*/

		x = xd[i];
		y = yd[i];
		z = zd[i];

		x += deltaX;
		y += deltaY;
		z += deltaZ;

		if (checkLocation(x, y, z) == false){
			curandState localState;
			localState = rndStatesd[IND];

			do {
				x = 1.1f + curand_uniform(&localState) * (devCParBH.lx - 2.2f);
				y = 1.1f + curand_uniform(&localState) * (devCParBH.ly - 2.2f);
				z = 1.1f + curand_uniform(&localState) * (devCParBH.lz - 2.2f);
			} while(checkLocationFull(x, y, z, i) == false);

			printf("Particle[%d] out of borders: moved.\n", i);
			x1d[i] = 0.0f;
			y1d[i] = 0.0f;
			z1d[i] = 0.0f;
			rndStatesd[IND] = localState;
		}


		xd[i] = x;
		yd[i] = y;
		zd[i] = z;
	}
}

__global__
__launch_bounds__(THREADS0, FACTOR0)
void getCurrentPhysKernel(float* ringPhy) {
	register int i, inc;
	float x, y, phy;
	// iterate over all nodes assigned to thread
	inc = blockDim.x * gridDim.x;
	for (i = threadIdx.x + blockIdx.x * blockDim.x; i < nbodiesd; i += inc) {
		x = xd[i] - devCParBH.barrelR;
		y = yd[i] - devCParBH.barrelR;
		phy = atan2f(y, x);
		ringPhy[i] = phy;
	}
}

__global__
__launch_bounds__(THREADS0, FACTOR0)
void getRingDistrKernel(float* oldPhy, float* newPhy, int nRings, float* dphy, int* partPos) {
	float dPhyCached, r;
	int interval;
	int inc = blockDim.x * gridDim.x;
	//int ind = IND;
	float x, y;
	__shared__ float dPhyAddayShared[SHARED_ARRAY_BH];	//nRings should not be bigger than SHARED_ARRAY;
	__shared__ int sharedCounter[SHARED_ARRAY_BH];

	if (threadIdx.x < nRings) {
		dPhyAddayShared[threadIdx.x] = 0.0f;
		sharedCounter[threadIdx.x] = 0;
	}
	if (IND < nRings) {
		dphy[IND] = 0.0f;
		partPos[IND] = 0;
	}
	__syncthreads();

	for (int i = IND; i < nbodiesd; i += inc) {
		x = xd[i];
		y = yd[i];

		dPhyCached = newPhy[i] - oldPhy[i];
		if (dPhyCached > M_PI) dPhyCached -= 2.0f * M_PI;
		if (dPhyCached < -M_PI) dPhyCached += 2.0f * M_PI;

		r = sqrtf(SQ(x - devCParBH.barrelR) + SQ(y - devCParBH.barrelR)/*, 0.5f*/);
		interval = (int) floorf((float)nRings * r / devCParBH.barrelR);	//ly - is R of barrel
		atomicAdd(&(dPhyAddayShared[interval]), dPhyCached);
		atomicAdd(&(sharedCounter[interval]), 1);
	}
	__syncthreads();

	if (threadIdx.x == 0) {
		for (int i = 0; i < nRings; i++) {
			if (sharedCounter[i] > 0) atomicAdd(&(dphy[i]), dPhyAddayShared[i] / (sharedCounter[i] * gridDim.x));
			atomicAdd(&(partPos[i]), sharedCounter[i]);
		}
	}
}



static void CudaTest(char *msg) {
	cudaError_t e;

	cudaThreadSynchronize();
	if (cudaSuccess != (e = cudaGetLastError())) {
		fprintf(stderr, "%s: %d\n", msg, e);
		fprintf(stderr, "%s\n", cudaGetErrorString(e));
		exit(-1);
	}
}

void init(int blocks, curandState* devRndStates) {
	initializationKernel<<<blocks * FACTOR5, THREADS5>>>(devRndStates, time(NULL));
	CudaTest("init kernel launch failed");

	cudaFuncSetCacheConfig(averageU, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(calcUKernel, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(boundingBoxKernel, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(boundingBoxKernel, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(treeBuildingKernel, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(summarizationKernel, cudaFuncCachePreferShared);
	cudaFuncSetCacheConfig(sortKernel, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(forceCalculationKernel, cudaFuncCachePreferL1);
	cudaFuncSetCacheConfig(integrationKernel, cudaFuncCachePreferL1);

	cudaGetLastError();	// reset error value
}

void calcU(int blocks) {
	calcUKernel<<<blocks * FACTOR0, THREADS0>>>();
	CudaTest("kernel 0 launch failed");
	cudaThreadSynchronize();
}

void buildBox(int blocks) {
	boundingBoxKernel<<<blocks * FACTOR1, THREADS1>>>();
	CudaTest("kernel 1 launch failed");
	cudaThreadSynchronize();
}

void buildTree(int blocks) {
	treeBuildingKernel<<<blocks * FACTOR2, THREADS2>>>();
	CudaTest("kernel 2 launch failed");
	cudaThreadSynchronize();
}

void summarize(int blocks) {
	summarizationKernel<<<blocks * FACTOR3, THREADS3>>>();
	CudaTest("kernel 3 launch failed");
	cudaThreadSynchronize();
}

void sort(int blocks) {
	sortKernel<<<blocks * FACTOR4, THREADS4>>>();
	CudaTest("kernel 4 launch failed");
	cudaThreadSynchronize();
}

void force(int blocks) {
	forceCalculationKernel<<<blocks * FACTOR5, THREADS5>>>();
	CudaTest("kernel 5 launch failed");
	cudaThreadSynchronize();
}

void integrate(int blocks) {
	integrationKernel<<<blocks * FACTOR6, THREADS6>>>();
	CudaTest("kernel 6 launch failed");
	cudaThreadSynchronize();
}

void fillConstantPointers(int nbodies, int nnodes, float* mass, BHArrays arrl, Box boxl, PartParams devMatrixes) {
	float epssq = 0.05f * 0.05f;
	float itolsq = 1.0f / (0.5f * 0.5f);

	if (cudaSuccess != cudaMemcpyToSymbol(nnodesd, &nnodes, sizeof(int), 0, cudaMemcpyHostToDevice))
		throw DeviceMemCpyToSymbolException("copying of nnodes to device failed\n");	//CudaTest("nnode copy to device failed");
	if (cudaSuccess != cudaMemcpyToSymbol(nbodiesd, &nbodies, sizeof(int), 0, cudaMemcpyHostToDevice))
		throw DeviceMemCpyToSymbolException("copying of nnodes to device failed\n");	//CudaTest("nnode copy to device failed");

	if (cudaSuccess != cudaMemcpyToSymbol(epssqd, &epssq, sizeof(float), 0, cudaMemcpyHostToDevice))
		throw DeviceMemCpyToSymbolException("copying of epssq to device failed\n");	//CudaTest("epssq copy to device failed");
	if (cudaSuccess != cudaMemcpyToSymbol(itolsqd, &itolsq, sizeof(float), 0, cudaMemcpyHostToDevice))
		throw DeviceMemCpyToSymbolException("copying of itolsq to device failed\n");	//CudaTest("itolsq copy to device failed");

	if (cudaSuccess != cudaMemcpyToSymbol(errd, &(arrl.err), sizeof(void*), 0, cudaMemcpyHostToDevice))
				throw DeviceMemCpyToSymbolException("copying of arrl.err to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(sortd, &(arrl.sort), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of arrl.err to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(childd, &(arrl.child), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of arrl.err to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(countd, &(arrl.count), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of arrl.err to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(startd, &(arrl.start), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of arrl.err to device failed\n");

	if (cudaSuccess != cudaMemcpyToSymbol(deltaPhyd, &(devMatrixes.deltaPhy), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaPhy to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(deltaThetad, &(devMatrixes.deltaTheta), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaTheta to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(phyd, &(devMatrixes.phy), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaPhy to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(thetad, &(devMatrixes.theta), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.Theta to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(phy1d, &(devMatrixes.phy1), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.Phy1 to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(theta1d, &(devMatrixes.theta1), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.Theta1 to device failed\n");

	if (cudaSuccess != cudaMemcpyToSymbol(deltaXd, &(devMatrixes.deltaX), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaX to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(deltaYd, &(devMatrixes.deltaY), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaY to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(deltaZd, &(devMatrixes.deltaZ), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaZ to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(xd, &(devMatrixes.x), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.x to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(yd, &(devMatrixes.y), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.y to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(zd, &(devMatrixes.z), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.z to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(x1d, &(devMatrixes.x1), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.x1 to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(y1d, &(devMatrixes.y1), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.y1 to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(z1d, &(devMatrixes.z1), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.z1 to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(uxd, &(devMatrixes.ux), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.ux to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(uyd, &(devMatrixes.uy), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.uy to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(uzd, &(devMatrixes.uz), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.uz to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(rndStatesd, &(devMatrixes.rndStates), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaPhy to device failed\n");

	if (cudaSuccess != cudaMemcpyToSymbol(massd, &(mass), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of devMatrixes.deltaPhy to device failed\n");

	if (cudaSuccess != cudaMemcpyToSymbol(maxxd, &(boxl.maxx), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of boxl.maxx to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(maxyd, &(boxl.maxy), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of boxl.maxx to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(maxzd, &(boxl.maxz), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of boxl.maxx to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(minxd, &(boxl.minx), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of boxl.maxx to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(minyd, &(boxl.miny), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of boxl.maxx to device failed\n");
	if (cudaSuccess != cudaMemcpyToSymbol(minzd, &(boxl.minz), sizeof(void*), 0, cudaMemcpyHostToDevice))
			throw DeviceMemCpyToSymbolException("copying of boxl.maxx to device failed\n");
}


void fillGloabalChangableBH(ChangableParams* chPar) {
	cudaMemcpyToSymbol(devChParBH, chPar, sizeof(ChangableParams), 0, cudaMemcpyHostToDevice);
}

void fillGloabalConstantBH(ConstParams* cPar) {
	cudaMemcpyToSymbol(devCParBH, cPar, sizeof(ConstParams), 0, cudaMemcpyHostToDevice);
}

float3 averageU(int blocks) {
	float3 u;
	float3* devU;
	if (cudaSuccess != cudaMalloc((void**)&devU, sizeof(float3)))
		throw DeviceMemoryAllocationException ("devU Allocation exception");
	averageUKernel<<<blocks * FACTOR0, THREADS0>>>(devU);
	CudaTest("AverageU launch failed");
	if (cudaSuccess != cudaMemcpy(&u, devU, sizeof(float3), cudaMemcpyDeviceToHost))
		throw DeviceMemoryCopyException ("devU copy exception");

	cudaFree(devU);
	return u;
}

void getCurrentPhysBH(int blocks, float* devOldPhy) {
	getCurrentPhysKernel<<<blocks * FACTOR0, THREADS0>>> (devOldPhy);
}

void getRingStatBH(int blocks, int nRings, int nPart, float** devOldPhy, float* dphy, int* pd) {
	float* devNewPhy;
	float* devDPhy;
	int* devPD;
	float* tmp;

	if (cudaSuccess != cudaMalloc((void**)&devNewPhy, sizeof(float) * nPart))
		throw DeviceMemoryAllocationException("Error allocation of devNewPhy");
	if (cudaSuccess != cudaMalloc((void**)&devDPhy, sizeof(float) * nRings))
		throw DeviceMemoryAllocationException("Error allocation of devDPhy");
	if (cudaSuccess != cudaMalloc((void**)&devPD, sizeof(int) * nRings))
		throw DeviceMemoryAllocationException("Error allocation of devPD");

	getCurrentPhysBH(blocks, devNewPhy);
	getRingDistrKernel<<<blocks * FACTOR0, THREADS0>>>(*devOldPhy, devNewPhy, nRings, devDPhy, devPD);

	cudaMemcpy(dphy, devDPhy, sizeof(float) * nRings, cudaMemcpyDeviceToHost);
	cudaMemcpy(pd, devPD, sizeof(int) * nRings, cudaMemcpyDeviceToHost);

	tmp = *devOldPhy;
	*devOldPhy = devNewPhy;

	if (cudaSuccess != cudaFree(tmp))
		throw DeviceMemoryException(std::string("Error free of tmp") +
				std::string(cudaGetErrorString(cudaGetLastError())));
	if (cudaSuccess != cudaFree(devDPhy))
		throw DeviceMemoryException(std::string("Error free of devDPhy") +
				std::string(cudaGetErrorString(cudaGetLastError())));
	if (cudaSuccess != cudaFree(devPD))
		throw DeviceMemoryException(std::string("Error free of devPD") +
				std::string(cudaGetErrorString(cudaGetLastError())));
}
