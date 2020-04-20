/*
 * This file is part of the BSGS distribution (https://github.com/JeanLucPons/BSGS).
 * Copyright (c) 2020 Jean Luc PONS.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include "HashTable.h"
#include <stdio.h>
#include <math.h>
#ifndef WIN64
#include <string.h>
#endif

#define safe_free(x) if(x) {free(x);x=NULL;}


HashTable::HashTable() {

  memset(E,0,sizeof(E));
  memset(M,0,sizeof(M));
  memset(C,0,sizeof(C));

}

void HashTable::Reset() {

  memset(M,0,sizeof(M));
  memset(C,0,sizeof(C));
  for (uint32_t h = 0; h < HASH_SIZE; h++) {
    safe_free(E[h]);
  }

}

uint64_t HashTable::GetNbItem() {

  uint64_t totalItem = 0;
  for(uint64_t h = 0; h < HASH_SIZE; h++) 
    if(E[h])
      totalItem += (uint64_t)C[h];

  return totalItem;
}

void HashTable::Add(Int *x,uint64_t p) {

  uint64_t h = (x->bits64[3] & HASH_MASK);

  if (C[h]>=M[h]) {
    // (re)allocate
    M[h] += 4;
    ENTRY *old = E[h];
    E[h] = (ENTRY *)malloc(sizeof(ENTRY) * M[h]);
    if(old) {
      memcpy(E[h],old,C[h]*sizeof(ENTRY));
      safe_free(old);
    }
  }

  E[h][C[h]].b0 = (uint32_t)(x->bits64[0] & B0MASK) | (uint32_t)((p & 0xF00000000ULL) >> 32);
  E[h][C[h]].p  = (uint32_t)(p & 0xFFFFFFFFULL);
  C[h]++;

}

bool HashTable::Get(Int *x,std::vector<uint64_t> &p) {

  p.clear();
  uint64_t h = (x->bits64[3] & HASH_MASK);

  // Find the leftmost element
  uint16_t L = 0;
  uint16_t R = C[h];
  uint16_t m;
  uint32_t T = (uint32_t)(x->bits64[0] & B0MASK);

  while (L < R) {
    m = (L + R) / 2;
    if ( (E[h][m].b0 & B0MASK) < T ) {
      L = m + 1;
    } else {
      R = m;
    }
  }

  if( L<C[h] && (E[h][L].b0 & B0MASK)==T ) {
    p.push_back((uint64_t)E[h][L].p + (((uint64_t)E[h][L].b0 & ~B0MASK)<<32));
    L++;
    while( L<C[h] && ((E[h][L].b0 & B0MASK) == T) ) {
      p.push_back((uint64_t)E[h][L].p + (((uint64_t)E[h][L].b0 & ~B0MASK)<<32));
      L++;
    }
  }

  return p.size()>0;

}
 
void HashTable::Sort(uint32_t h) {

  if(C[h]>1) {

    ENTRY *A = E[h];
    bool endOfSort = false;
    uint16_t s = C[h]-1;
    while(!endOfSort) {
      endOfSort = true;
      for(uint16_t j=0;j<s;j++) {
        if( (A[j].b0 & B0MASK) > (A[j+1].b0 & B0MASK) ) {
          ENTRY T = A[j];
          A[j] = A[j+1];
          A[j+1] = T;
          endOfSort = false;
        }          
      }
      s--;
    }

  }

}

double HashTable::GetSizeMB() {

  uint64_t byte = sizeof(E) + sizeof(C) + sizeof(M);

  for (int h = 0; h < HASH_SIZE; h++) {
    byte += sizeof(ENTRY) * M[h];
  }

  return (double)byte / (1024.0*1024.0);

}

void HashTable::PrintInfo() {

  uint16_t max = 0;
  uint32_t maxH = 0;
  uint16_t min = 65535;
  uint32_t minH = 0;
  double std = 0;
  double avg = (double)GetNbItem() / (double)HASH_SIZE;

  for(uint32_t h=0;h<HASH_SIZE;h++) {
    if(C[h]>max) {
      max=C[h];
      maxH = h;
    }
    if(C[h]<min) {
      min=C[h];
      minH = h;
    }
    std += (avg - (double)C[h])*(avg - (double)C[h]);
  }
  std /= (double)HASH_SIZE;
  std = sqrt(std);

  ::printf("Size: %f MB\n",GetSizeMB());
  ::printf("Item: 2^%.2f \n",log2((double)GetNbItem()));
  ::printf("Max : %d [@ %06X]\n",max,maxH);
  ::printf("Min : %d [@ %06X]\n",min,minH);
  ::printf("Avg : %.2f \n",avg);
  ::printf("SDev: %.2f \n",std);

}
