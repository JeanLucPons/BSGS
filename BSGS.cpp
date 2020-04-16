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

#include "BSGS.h"
#include <fstream>
#include "SECPK1/IntGroup.h"
#include "Timer.h"
#include <string.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <algorithm>
#ifndef WIN64
#include <pthread.h>
#endif

#ifdef WIN64
#define LOCK(mutex) WaitForSingleObject(mutex,INFINITE);
#define UNLOCK(mutex) ReleaseMutex(mutex);
#else
#define LOCK(mutex) pthread_mutex_lock(&(mutex));
#define UNLOCK(mutex) pthread_mutex_unlock(&(mutex));
#endif

using namespace std;

// Baby steps
Point Gn[CPU_GRP_SIZE / 2];
Point _2Gn;

// Giant steps
Point GSn[CPU_GRP_SIZE / 2];
Point _2GSn;

// ----------------------------------------------------------------------------

BSGS::BSGS(Secp256K1 *secp) {

  this->secp = secp;

  // Compute Generator table G[n] = (n+1)*G (Baby steps)
  Point g = secp->G;
  Gn[0] = g;
  g = secp->DoubleDirect(g);
  Gn[1] = g;
  for(int i = 2; i < CPU_GRP_SIZE / 2; i++) {
    g = secp->AddDirect(g,secp->G);
    Gn[i] = g;
  }
  // _2Gn = CPU_GRP_SIZE*G
  _2Gn = secp->DoubleDirect(Gn[CPU_GRP_SIZE / 2 - 1]);

  // Init mutex
#ifdef WIN64
  ghMutex = CreateMutex(NULL,FALSE,NULL);
#else
  ghMutex = PTHREAD_MUTEX_INITIALIZER;
#endif

}

// ----------------------------------------------------------------------------

bool BSGS::ParseConfigFile(std::string fileName) {

  // Check file
  FILE *fp = fopen(fileName.c_str(),"rb");
  if(fp == NULL) {
    ::printf("Error: Cannot open %s %s\n",fileName.c_str(),strerror(errno));
    return false;
  }
  fclose(fp);

  // Get lines
  vector<string> lines;
  int nbLine = 0;
  string line;
  ifstream inFile(fileName);
  while(getline(inFile,line)) {

    // Remove ending \r\n
    int l = (int)line.length() - 1;
    while(l >= 0 && isspace(line.at(l))) {
      line.pop_back();
      l--;
    }

    if(line.length() > 0) {
      lines.push_back(line);
      nbLine++;
    }

  }

  if(lines.size()<4) {
    ::printf("Error: %s not enough arguments\n",fileName.c_str());
    return false;
  }

  bsSize = ::stoll(lines[0],NULL,16);
  rangeStart.SetBase16((char *)lines[1].c_str());
  rangeEnd.SetBase16((char *)lines[2].c_str());
  for(int i=3;i<(int)lines.size();i++) {
    
    Point p;
    bool isCompressed;
    if( !secp->ParsePublicKeyHex(lines[i],p,isCompressed) ) {
      ::printf("%s, error line %d: %s\n",fileName.c_str(),i,lines[i].c_str());
    }
    keysToSearch.push_back(p);

  }

#ifdef WIN64
  ::printf("BabyStep:0x%016I64X\n",bsSize);
#else
  ::printf("BabyStep:0x%" PRIx64 "\n",bsSize);
#endif

  ::printf("Start:%s\n",rangeStart.GetBase16().c_str());
  ::printf("Stop :%s\n",rangeEnd.GetBase16().c_str());
  ::printf("Keys :%d\n",(int)keysToSearch.size());

  return true;

}

// ----------------------------------------------------------------------------
#ifdef WIN64

THREAD_HANDLE BSGS::LaunchThread(LPTHREAD_START_ROUTINE func, TH_PARAM *p) {
  p->obj = this;
  return CreateThread(NULL, 0, func, (void*)(p), 0, NULL);
}
void  BSGS::JoinThreads(THREAD_HANDLE *handles,int nbThread) {
  WaitForMultipleObjects(nbThread, handles, TRUE, INFINITE);
}
void  BSGS::FreeHandles(THREAD_HANDLE *handles, int nbThread) {
  for (int i = 0; i < nbThread; i++)
    CloseHandle(handles[i]);
}
#else

THREAD_HANDLE BSGS::LaunchThread(void *(*func) (void *), TH_PARAM *p) {
  THREAD_HANDLE h;
  p->obj = this;
  pthread_create(&h, NULL, func, (void*)(p));
  return h;
}
void  BSGS::JoinThreads(THREAD_HANDLE *handles, int nbThread) {
  for (int i = 0; i < nbThread; i++)
    pthread_join(handles[i], NULL);
}
void  BSGS::FreeHandles(THREAD_HANDLE *handles, int nbThread) {
}
#endif


// ----------------------------------------------------------------------------

void BSGS::FillBabySteps(TH_PARAM *ph) {

  // Global init
  int thId = ph->threadId;
  counters[thId] = 0;
  uint64_t nbStep = kPerThread / CPU_GRP_SIZE;

  // CPU Thread
  IntGroup *grp = new IntGroup(CPU_GRP_SIZE / 2 + 1);

  // Group Init
  Point startP;

  Int dx[CPU_GRP_SIZE / 2 + 1];
  Point pts[CPU_GRP_SIZE];

  Int dy;
  Int dyn;
  Int _s;
  Int _p;
  Point pp;
  Point pn;
  grp->Set(dx);

  Int km(&ph->startKey);
  km.Add((uint64_t)CPU_GRP_SIZE / 2);
  startP = secp->ComputePublicKey(&km);

  ph->hasStarted = true;

#ifdef WIN64
    ::printf("BabyStep Thread %d: 0x%016I64X -> 0x%016I64X\n",ph->threadId,ph->startKey.bits64[0],ph->startKey.bits64[0]+kPerThread-1);
#else
    ::printf("BabyStep Thread %d: 0x%" PRIx64 " -> 0x%" PRIx64 "\n",ph->threadId,ph->startKey.bits64[0],ph->startKey.bits64[0]+kPerThread-1);
#endif

  // Baby Step Hashtable contains G,2G,3G,.....,(bsSize).G

  for(uint64_t s=0;s<nbStep;s++) {

    // Fill group
    int i;
    int hLength = (CPU_GRP_SIZE / 2 - 1);

    for(i = 0; i < hLength; i++) {
      dx[i].ModSub(&Gn[i].x,&startP.x);
    }
    dx[i].ModSub(&Gn[i].x,&startP.x);  // For the first point
    dx[i + 1].ModSub(&_2Gn.x,&startP.x); // For the next center point

    // Grouped ModInv
    grp->ModInv();

    // We use the fact that P + i*G and P - i*G has the same deltax, so the same inverse
    // We compute key in the positive and negative way from the center of the group

    // center point
    pts[CPU_GRP_SIZE / 2] = startP;

    for(i = 0; i<hLength; i++) {

      pp = startP;
      pn = startP;

      // P = startP + i*G
      dy.ModSub(&Gn[i].y,&pp.y);

      _s.ModMulK1(&dy,&dx[i]);       // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pp.x.ModNeg();
      pp.x.ModAdd(&_p);
      pp.x.ModSub(&Gn[i].x);           // rx = pow2(s) - p1.x - p2.x;

#if 0
      pp.y.ModSub(&Gn[i].x,&pp.x);
      pp.y.ModMulK1(&_s);
      pp.y.ModSub(&Gn[i].y);           // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      // P = startP - i*G  , if (x,y) = i*G then (x,-y) = -i*G
      dyn.Set(&Gn[i].y);
      dyn.ModNeg();
      dyn.ModSub(&pn.y);

      _s.ModMulK1(&dyn,&dx[i]);      // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pn.x.ModNeg();
      pn.x.ModAdd(&_p);
      pn.x.ModSub(&Gn[i].x);          // rx = pow2(s) - p1.x - p2.x;

#if 0
      pn.y.ModSub(&Gn[i].x,&pn.x);
      pn.y.ModMulK1(&_s);
      pn.y.ModAdd(&Gn[i].y);          // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      pts[CPU_GRP_SIZE / 2 + (i + 1)] = pp;
      pts[CPU_GRP_SIZE / 2 - (i + 1)] = pn;

    }

    // First point (startP - (GRP_SZIE/2)*G)
    pn = startP;
    dyn.Set(&Gn[i].y);
    dyn.ModNeg();
    dyn.ModSub(&pn.y);

    _s.ModMulK1(&dyn,&dx[i]);
    _p.ModSquareK1(&_s);

    pn.x.ModNeg();
    pn.x.ModAdd(&_p);
    pn.x.ModSub(&Gn[i].x);

#if 0
    pn.y.ModSub(&Gn[i].x,&pn.x);
    pn.y.ModMulK1(&_s);
    pn.y.ModAdd(&Gn[i].y);
#endif

    pts[0] = pn;

    // Add to table
    LOCK(ghMutex);
    for(uint64_t i=0;i<CPU_GRP_SIZE;i++)
      hashTable.Add(&pts[i].x, (uint64_t)( ph->startKey.bits64[0] + s*CPU_GRP_SIZE + i) );
    UNLOCK(ghMutex);
    counters[thId] += CPU_GRP_SIZE;

    // Next start point (startP + GRP_SIZE*G)
    pp = startP;
    dy.ModSub(&_2Gn.y,&pp.y);

    _s.ModMulK1(&dy,&dx[i + 1]);
    _p.ModSquareK1(&_s);

    pp.x.ModNeg();
    pp.x.ModAdd(&_p);
    pp.x.ModSub(&_2Gn.x);

    pp.y.ModSub(&_2Gn.x,&pp.x);
    pp.y.ModMulK1(&_s);
    pp.y.ModSub(&_2Gn.y);
    startP = pp;

  }

  delete grp;
  ph->isRunning = false;

}

// ----------------------------------------------------------------------------

void BSGS::SolveKey(TH_PARAM *ph) {

  // Global init
  int thId = ph->threadId;
  counters[thId] = 0;
  vector<uint64_t> off;

  // CPU Thread
  IntGroup *grp = new IntGroup(CPU_GRP_SIZE / 2 + 1);

  // Group Init
  Point startP;

  Int dx[CPU_GRP_SIZE / 2 + 1];
  Point pts[CPU_GRP_SIZE];

  Int dy;
  Int dyn;
  Int _s;
  Int _p;
  Point pp;
  Point pn;
  grp->Set(dx);

  // Substart startRange to the point to solve
  Int km(&ph->startKey);
  km.Neg();
  km.Add(&secp->order);
  km.Sub((uint64_t)(CPU_GRP_SIZE/2)*bsSize);
  startP = secp->ComputePublicKey(&km);
  startP = secp->AddDirect(keyToSearch,startP);

  ph->hasStarted = true;

  if(keyIdx==0)
    ::printf("GiantStep Thread %d: %s\n",ph->threadId,ph->startKey.GetBase16().c_str());

  // Substart ((s*CPU_GRP_SIZE+i)*bsSize).G to the point to solve and look for a match into the hashtable

  for(uint64_t s = 0; s<ph->nbStep && !endOfSearch; s++) {

    // Fill group
    int i;
    int hLength = (CPU_GRP_SIZE / 2 - 1);

    for(i = 0; i < hLength; i++) {
      dx[i].ModSub(&GSn[i].x,&startP.x);
    }
    dx[i].ModSub(&GSn[i].x,&startP.x);  // For the first point
    dx[i + 1].ModSub(&_2GSn.x,&startP.x); // For the next center point

    // Grouped ModInv
    grp->ModInv();

    // We use the fact that P + i*G and P - i*G has the same deltax, so the same inverse
    // We compute key in the positive and negative way from the center of the group

    // center point
    pts[CPU_GRP_SIZE / 2] = startP;

    for(i = 0; i<hLength; i++) {

      pp = startP;
      pn = startP;

      // P = startP + i*G
      dy.ModSub(&GSn[i].y,&pp.y);

      _s.ModMulK1(&dy,&dx[i]);       // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pp.x.ModNeg();
      pp.x.ModAdd(&_p);
      pp.x.ModSub(&GSn[i].x);           // rx = pow2(s) - p1.x - p2.x;

#if 0
      pp.y.ModSub(&GSn[i].x,&pp.x);
      pp.y.ModMulK1(&_s);
      pp.y.ModSub(&GSn[i].y);           // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      // P = startP - i*G  , if (x,y) = i*G then (x,-y) = -i*G
      dyn.Set(&GSn[i].y);
      dyn.ModNeg();
      dyn.ModSub(&pn.y);

      _s.ModMulK1(&dyn,&dx[i]);      // s = (p2.y-p1.y)*inverse(p2.x-p1.x);
      _p.ModSquareK1(&_s);            // _p = pow2(s)

      pn.x.ModNeg();
      pn.x.ModAdd(&_p);
      pn.x.ModSub(&GSn[i].x);          // rx = pow2(s) - p1.x - p2.x;

#if 0
      pn.y.ModSub(&GSn[i].x,&pn.x);
      pn.y.ModMulK1(&_s);
      pn.y.ModAdd(&GSn[i].y);          // ry = - p2.y - s*(ret.x-p2.x);  
#endif

      pts[CPU_GRP_SIZE / 2 + (i + 1)] = pp;
      pts[CPU_GRP_SIZE / 2 - (i + 1)] = pn;

    }

    // First point (startP - (GRP_SZIE/2)*G)
    pn = startP;
    dyn.Set(&GSn[i].y);
    dyn.ModNeg();
    dyn.ModSub(&pn.y);

    _s.ModMulK1(&dyn,&dx[i]);
    _p.ModSquareK1(&_s);

    pn.x.ModNeg();
    pn.x.ModAdd(&_p);
    pn.x.ModSub(&GSn[i].x);

#if 0
    pn.y.ModSub(&GSn[i].x,&pn.x);
    pn.y.ModMulK1(&_s);
    pn.y.ModAdd(&GSn[i].y);
#endif

    pts[0] = pn;

    // Check key
    for(uint64_t i = 0; i<CPU_GRP_SIZE; i++) {
      if( hashTable.Get(&pts[i].x,off) ) {
        for(int o=0;o<off.size();o++) {
         Int pk(bsSize);
         Int bigS(s*CPU_GRP_SIZE+i);
         pk.Mult(&bigS);
         Int bigO(off[o]);
         pk.Add(&bigO);
         pk.Add(&ph->startKey);
         // Check
         Point p = secp->ComputePublicKey(&pk);
         if(p.equals(keyToSearch)) {
           // Key solved
           ::printf("\nKey#%2d Pub:  0x%s \n",keyIdx,secp->GetPublicKeyHex(true,p).c_str());
           ::printf("       Priv: 0x%s \n",pk.GetBase16().c_str());
           endOfSearch = true;
         }
        }
      }
    }
    counters[thId] += CPU_GRP_SIZE;

    // Next start point (startP + GRP_SIZE*G)
    pp = startP;
    dy.ModSub(&_2GSn.y,&pp.y);

    _s.ModMulK1(&dy,&dx[i + 1]);
    _p.ModSquareK1(&_s);

    pp.x.ModNeg();
    pp.x.ModAdd(&_p);
    pp.x.ModSub(&_2GSn.x);

    pp.y.ModSub(&_2GSn.x,&pp.x);
    pp.y.ModMulK1(&_s);
    pp.y.ModSub(&_2GSn.y);
    startP = pp;

  }

  delete grp;
  ph->isRunning = false;

}

// ----------------------------------------------------------------------------

#ifdef WIN64
DWORD WINAPI _FillBS(LPVOID lpParam) {
#else
void *_FillBS(void *lpParam) {
#endif
  TH_PARAM *p = (TH_PARAM *)lpParam;
  p->obj->FillBabySteps(p);
  return 0;
}

// ----------------------------------------------------------------------------

#ifdef WIN64
DWORD WINAPI _SolveKey(LPVOID lpParam) {
#else
void *_SolveKey(void *lpParam) {
#endif
  TH_PARAM *p = (TH_PARAM *)lpParam;
  p->obj->SolveKey(p);
  return 0;
}

// ----------------------------------------------------------------------------

bool BSGS::isAlive(TH_PARAM *p) {

  bool isAlive = false;
  int total = nbCPUThread;
  for(int i=0;i<total;i++)
    isAlive = isAlive || p[i].isRunning;

  return isAlive;

}

// ----------------------------------------------------------------------------

bool BSGS::hasStarted(TH_PARAM *p) {

  bool hasStarted = true;
  int total = nbCPUThread;
  for (int i = 0; i < total; i++)
    hasStarted = hasStarted && p[i].hasStarted;

  return hasStarted;

}

// ----------------------------------------------------------------------------

bool BSGS::isWaiting(TH_PARAM *p) {

  bool isWaiting = true;
  int total = nbCPUThread;
  for (int i = 0; i < total; i++)
    isWaiting = isWaiting && p[i].isWaiting;

  return isWaiting;

}

// ----------------------------------------------------------------------------

uint64_t BSGS::getCPUCount() {

  uint64_t count = 0;
  for(int i=0;i<nbCPUThread;i++)
    count += counters[i];
  return count;

}

// ----------------------------------------------------------------------------

string BSGS::GetTimeStr(double dTime) {

  char tmp[256];

  double nbDay = dTime / 86400.0;
  if (nbDay >= 1) {

    double nbYear = nbDay / 365.0;
    if (nbYear > 1) {
      if (nbYear < 5)
        sprintf(tmp, "%.1fy", nbYear);
      else
        sprintf(tmp, "%gy", nbYear);
    } else {
      sprintf(tmp, "%.1fd", nbDay);
    }

  } else {

    int iTime = (int)dTime;
    int nbHour = (int)((iTime % 86400) / 3600);
    int nbMin = (int)(((iTime % 86400) % 3600) / 60);
    int nbSec = (int)(iTime % 60);

    if (nbHour == 0) {
      if (nbMin == 0) {
        sprintf(tmp, "%02ds", nbSec);
      } else {
        sprintf(tmp, "%02d:%02d", nbMin, nbSec);
      }
    } else {
      sprintf(tmp, "%02d:%02d:%02d", nbHour, nbMin, nbSec);
    }

  }

  return string(tmp);

}

// ----------------------------------------------------------------------------

void BSGS::Run(int nbThread) {

  double t0 = Timer::get_tick();

  nbCPUThread = nbThread;
  endOfSearch = false;

  TH_PARAM *params = (TH_PARAM *)malloc(nbCPUThread * sizeof(TH_PARAM));
  THREAD_HANDLE *thHandles = (THREAD_HANDLE *)malloc(nbCPUThread * sizeof(THREAD_HANDLE));

  memset(params, 0, nbCPUThread * sizeof(TH_PARAM));
  memset(counters, 0, sizeof(counters));
  printf("Number of CPU thread: %d\n", nbCPUThread);

  // Check input parameters
  uint64_t k = 1;
  kPerThread = bsSize/nbThread;
  if( bsSize%nbThread != 0 ) {
    int gSize = CPU_GRP_SIZE;
    ::printf("Warning, BS size is not a multiple of nbThread\n");
  }
  if(kPerThread % CPU_GRP_SIZE != 0) {
    int gSize = CPU_GRP_SIZE;
    ::printf("Warning, BSSize/nbThread is not a multiple of %d\n",gSize);
  }

  // Launch Baby Step threads
  for(int i = 0; i < nbCPUThread; i++) {
    params[i].threadId = i;
    params[i].isRunning = true;
    params[i].startKey.bits64[0] = k;
    thHandles[i] = LaunchThread(_FillBS,params + i);
    k += kPerThread;
  }

  // Wait for end of baby step calculation
  Process(params);
  JoinThreads(thHandles,nbCPUThread);
  FreeHandles(thHandles,nbCPUThread);

  // Compute range per thread
  Int bs(bsSize);
  Int nbTh;
  Int r;
  nbTh.SetInt32(nbThread);
  Int rgPerTh(&rangeEnd);
  rgPerTh.Sub(&rangeStart);
  rgPerTh.AddOne();
  rgPerTh.Div(&nbTh,&r);
  if(!r.IsZero()) {
    ::printf("Warning, range is not a multiple of nbThread\n");
  }
  Int stepPerThred(&rgPerTh);
  Int grpSize;
  grpSize.SetInt32(CPU_GRP_SIZE);
  stepPerThred.Div(&bs,&r);
  if(!r.IsZero()) stepPerThred.AddOne();
  stepPerThred.Div(&grpSize,&r);
  if(!r.IsZero()) {
    ::printf("Warning, range is not a multiple of nbThread*%d\n",CPU_GRP_SIZE);
  }

  // Compute Giant steps adding table GSn[n] = -(n+1)*BS
  bs.Neg();
  bs.Add(&secp->order);

  Point bsP = secp->ComputePublicKey(&bs);
  Point g = bsP;
  GSn[0] = g;
  g = secp->DoubleDirect(g);
  GSn[1] = g;
  for(int i = 2; i < CPU_GRP_SIZE / 2; i++) {
    g = secp->AddDirect(g,bsP);
    GSn[i] = g;
  }
  // _2Gn = -CPU_GRP_SIZE*BS
  _2GSn = secp->DoubleDirect(GSn[CPU_GRP_SIZE / 2 - 1]);

  for(keyIdx =0; keyIdx<keysToSearch.size(); keyIdx++) {

    keyToSearch = keysToSearch[keyIdx];

    // Lanch Giant Step threads
    endOfSearch = false;
    Int sk(&rangeStart);
    for(int i = 0; i < nbCPUThread; i++) {
      params[i].threadId = i;
      params[i].isRunning = true;
      params[i].startKey.Set(&sk);
      params[i].nbStep=stepPerThred.bits64[0];
      thHandles[i] = LaunchThread(_SolveKey,params + i);
      sk.Add(&rgPerTh);
    }

    // Wait for end
    Process(params);
    JoinThreads(thHandles,nbCPUThread);
    FreeHandles(thHandles,nbCPUThread);

  }

  double t1 = Timer::get_tick();

  ::printf("\nDone: Total time %s \n" , GetTimeStr(t1-t0).c_str());

}

// Wait for end of threads and display stats
void BSGS::Process(TH_PARAM *params) {

  double t0;
  double t1;

#ifndef WIN64
  setvbuf(stdout,NULL,_IONBF,0);
#endif

  uint64_t lastCount = 0;
  double avgKeyRate;
  uint64_t count;

  // Key rate smoothing filter
#define FILTER_SIZE 8
  double lastkeyRate[FILTER_SIZE];
  uint32_t filterPos = 0;

  double keyRate = 0.0;

  memset(lastkeyRate,0,sizeof(lastkeyRate));

  // Wait that all threads have started
  while(!hasStarted(params))
    Timer::SleepMillis(50);

  t0 = Timer::get_tick();
  startTime = t0;

  while(isAlive(params)) {

    int delay = 2000;
    while(isAlive(params) && delay > 0) {
      Timer::SleepMillis(50);
      delay -= 50;
    }

    count = getCPUCount();

    t1 = Timer::get_tick();
    keyRate = (double)(count - lastCount) / (t1 - t0);
    lastkeyRate[filterPos%FILTER_SIZE] = keyRate;
    filterPos++;

    // KeyRate smoothing
    avgKeyRate = 0.0;
    uint32_t nbSample;
    for(nbSample = 0; (nbSample < FILTER_SIZE) && (nbSample < filterPos); nbSample++) {
      avgKeyRate += lastkeyRate[nbSample];
    }
    avgKeyRate /= (double)(nbSample);

    if(isAlive(params) && !endOfSearch) {

      ::printf("\r[%.2f MKs][Cnt 2^%.2f][%s][%.1fMB]  ",
        avgKeyRate / 1000000.0,
        log2((double)count),
        GetTimeStr(t1 - startTime).c_str(),
        hashTable.GetSizeMB()
        );

    }

    lastCount = count;
    t0 = t1;

  }

  count = hashTable.GetNbItem();
  
  if( !endOfSearch ) {
    ::printf("\r[%.2f MKs][Cnt 2^%.2f][%s][%.1fMB]  \n",
      avgKeyRate / 1000000.0,
      log2((double)count),
      GetTimeStr(t1 - startTime).c_str(),
      hashTable.GetSizeMB()
    );
  }

}
