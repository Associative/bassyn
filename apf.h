//
// S A M P L E   G E N E R A T I O N   &   P R O C E S S I N G
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// 
// Various waveform generators (sine, saw, triangle) and some
// filters (hp, lp, bp, shelve etc.)
//
#ifndef _APF_
#define _APF_

#include <stdio.h>
#include <math.h>
#include "SDL.h"

// synthesizer module types
enum moduleTypes {
  // basic functions
  MT_Value,
  MT_Sum,
  MT_Gain,

  // generator modules
  MT_G_Sine,
  MT_G_Triangle,
  MT_G_Square,

  // second order filter modules
  MT_F2_HighPass,
  MT_F2_LowPass,
  MT_F2_BandPass,
  MT_F2_HighShelve,
  MT_F2_LowShelve,
  MT_F2_Peak,
  MODULE_TYPE_COUNT
};

// basic module representation
struct module {
  int    type;            // module type
  int    level;           // level in processing hierarchy (determines calculation order of the modules)
  double value;           // output value
  struct module **param;  // array of pointers to parameter values (addresses of other modules)
  void   *data;           // module dependent data
  struct module *link;    // link in a linked list for sequential processing
};

// second order filter memory structure
struct filterMemory2 {
  double x[3];
  double y[3];
};

// module description structure
struct moduleDesc {
  const char *name;
  const int  paramCount;
  const char **paramStrings;
  const int  dataSize;
  void (*process)(struct module*);
};

// global variables
extern unsigned long samplerate;
extern unsigned long long samplecount;

// general functions
int calcLevel(struct module *m);
double processModules(struct module *ml);
struct module* sortModules(struct module *ml);
struct module* createModule(int t, double v, struct module *l);
struct module* freeModule(struct module *m);
const char* getModuleName(struct module *m);

double clipEx(double x, double max);
double clipIn(double x, double max);
double getPeriodSample(double f, double ph);



//  M O D U L E   D E F I N I T I O N S 
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// module type descriptions
extern const struct moduleDesc moduleInfo[];

// module parameter strings
extern const char *valueStrings[];
extern const char *sumStrings[];
extern const char *gainStrings[];
extern const char *gSineStrings[];
extern const char *gTriangleStrings[];
extern const char *gSquareStrings[];
extern const char *f2HighPassStrings[];
extern const char *f2LowPassStrings[];
extern const char *f2BandPassStrings[];
extern const char *f2HighShelveStrings[];
extern const char *f2LowShelveStrings[];
extern const char *f2PeakStrings[];

// module processing functions
void valueProcess(struct module *m);
void sumProcess(struct module *m);
void gainProcess(struct module *m);
void gSineProcess(struct module *m);
void gTriangleProcess(struct module *m);
void gSquareProcess(struct module *m);
void f2HighPassProcess(struct module *m);
void f2LowPassProcess(struct module *m);
void f2BandPassProcess(struct module *m);
void f2HighShelveProcess(struct module *m);
void f2LowShelveProcess(struct module *m);
void f2PeakProcess(struct module *m);



/*
// Value module
enum MP_Value {
  MP_Value_COUNT=0
};

// Sum module
enum MP_Sum {
  MP_Sum_InputA=0,
  MP_Sum_InputB,
  MP_Sum_COUNT
};

// Gain module
enum MP_Gain {
  MP_Gain_Input=0,
  MP_Gain_Factor,
  MP_Gain_COUNT
};

// Sine wave generator module
enum MP_G_Sine {
  MP_G_Sine_Freq=0,
  MP_G_Sine_Phase,
  MP_G_Sine_Amplitude,
  MP_G_Sine_COUNT
};

// Triangle wave generator module
enum MP_G_Triangle {
  MP_G_Triangle_Freq=0,
  MP_G_Triangle_Phase,
  MP_G_Triangle_Amplitude,
  MP_G_Triangle_Skew,
  MP_G_Triangle_COUNT
};

// Square wave generator module
enum MP_G_Square {
  MP_G_Square_Freq=0,
  MP_G_Square_Phase,
  MP_G_Square_Duty,
  MP_G_Square_COUNT
};

// Second order high-pass filter module
enum MP_F2_HighPass {
  MP_F2_HighPass_Bypass=0,
  MP_F2_HighPass_Freq,
  MP_F2_HighPass_COUNT
}; 

// Second order low-pass filter module
enum MP_F2_LowPass {
  MP_F2_LowPass_Bypass=0,
  MP_F2_LowPass_Freq,
  MP_F2_LowPass_COUNT
}; 

// second order band-pass filter module
enum MP_F2_BandPass {
  MP_F2_BandPass_Bypass=0,
  MP_F2_BandPass_Freq,
  MP_F2_BandPass_Q,
  MP_F2_BandPass_COUNT
}; 

// second order high-frequency shelving filter module
enum MP_F2_HighShelve {
  MP_F2_HighShelve_Bypass=0,
  MP_F2_HighShelve_Freq,
  MP_F2_HighShelve_Boost,
  MP_F2_HighShelve_COUNT
}; 

// second order low-frequency shelving filter module
enum MP_F2_LowShelve {
  MP_F2_LowShelve_Bypass=0,
  MP_F2_LowShelve_Freq,
  MP_F2_LowShelve_Boost,
  MP_F2_LowShelve_COUNT
}; 

// second order peak filter module
enum MP_F2_Peak {
  MP_F2_Peak_Bypass=0,
  MP_F2_Peak_Freq,
  MP_F2_Peak_Q,
  MP_F2_Peak_COUNT
}; 
*/






/*
// basic waveform types (other waveforms can be derived from these)
enum generatorTypes {
  GT_Sine,
  GT_Triangle,
  GT_Square
};

// filter types
enum filterTypes {
  FT_HighPass,
  FT_LowPass,
  FT_BandPass,
  FT_HighShelve,
  FT_LowShelve,
  FT_Peak
};



// waveform description
struct generatorParams {
  char enable;
  int type;
  double A;
  double f;
  double ph;
  double s;
};

// waveform generator functions (clipped at A>1)
double generator(struct generatorParams *p, Uint64 s, Uint32 ts);
//double genSinus(float f, float A, float ph, Uint64 s, Uint32 fs);
//double genTriangle(float f, float A, float ph, float d, Uint64 s, Uint32 fs);


// filter parameters
struct filterParams {
  char bypass;
  int type;
  int order;
  double fc;
  double *a;
  double *b;
};

// filter memory
struct filterMem {
  int size;
  double *x;
  double *y;
};

// filter function
double filter(struct filterParams *p, struct filterMem *m, double s);

// clear filter memory
void filter_clear(struct filterMem *m);

// shift filter memory
void filter_shift(struct filterMem *m);

// parameter calculations
void filter_setParams(struct filterParams *p, Uint32 fs);
*/


//  W a v e   g e n e r a t o r   f u n c t i o n s
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
// generic waveform generator
double generator(struct generatorParams *p, Uint64 s, Uint32 fs) {
  double x, y=0;
  double A;
  
  A = p->A;
  if (p->A<1.0f && p->A>0.0f) A = 1.0f;
  if (p->A<0.0f && p->A>-1.0f) A = -1.0f;

  x = fmod((2*M_PI*p->f*s)/fs + p->ph, 2*M_PI);

  switch (p->type) {
    case GT_Sine:
      y = A*sin(x);
      break;

    case GT_Triangle:
      if (x < 2*M_PI*p->s) y = -A + (x*A)/(M_PI*p->s);
      if (x >= 2*M_PI*p->s) y = -A - A*(x-2*M_PI)/(M_PI*(1-p->s));
      break;  
  }
  
  if (y > 1.0f) y = 1.0f;
  if (y < -1.0f) y = -1.0f;

  return y;
}
*/
/*
// sine wave generator
double genSinus(float f, float A, float ph, Uint64 s, Uint32 fs) {
  double y, x;

  if (A<1.0f && A>0.0f) A = 1.0f;
  if (A<0.0f && A>-1.0f) A = -1.0f;
  x = 2*M_PI*f*s/fs + ph;
  y = A*sin(x);
  if (y > 1.0f) y = 1.0f;
  if (y < -1.0f) y = -1.0f;

  return y;
} 

// triangle/saw wave generator
double genTriangle(float f, float A, float ph, float d, Uint64 s, Uint32 fs) {
  double x,y=0;

  if (A<1.0f && A>0.0f) A = 1.0f;
  if (A<0.0f && A>-1.0f) A = -1.0f;
  x = fmod((2*M_PI*f*s)/fs + ph, 2*M_PI);
  if (x < 2*d*M_PI) y = -A + (x*A)/(M_PI*d);
  if (x >= 2*d*M_PI) y = -A - A*(x-2*M_PI)/(M_PI*(1-d));
  if (y > 1.0f) y = 1.0f;
  if (y < -1.0f) y = -1.0f;  

  return y;
} 
*/


//  F i l t e r   f u n c t i o n s
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
/*
// The filter function
double filter(struct filterParams *p, struct filterMem *m, double s) {
  int i; 

  m->x[0] = s;
  m->y[0] = p->a[0] * m->x[0];

  for (i=1; i<=p->order; i++) {
    m->y[0] += (p->a[i] * m->x[i]) - (p->b[i] * m->y[i]);
  }
  
  filter_shift(m);
  return m->y[0];
}

// clear filter memory
void filter_clear(struct filterMem *m) {
  for (int i=0; i<m->size; i++) {
    m->x[i] = 0;
    m->y[i] = 0;
  } 
}

// shift filter memory
void filter_shift(struct filterMem *m) {
  for (int i=m->size-1; i>0; i--) {
    m->x[i] = m->x[i-1];
    m->y[i] = m->y[i-1];
  }
}

// parameter calculation
void filter_setParams(struct filterParams *p, Uint32 fs) {
  double K, denom;

  K = tan(M_PI*p->fc/fs);

  switch (p->type) {
  case FT_HighPass:
    break;

  case FT_LowPass:
    denom = 1 + sqrt(2)*K + K*K;
    p->a[0] = (K*K) / denom;
    p->a[1] = 2*p->a[0];
    p->a[2] = p->a[0];
    p->b[0] = 0;
    p->b[1] = p->a[1] - 2/denom;
    p->b[2] = (1-sqrt(2)*K+K*K) / denom;
    break;

  case FT_BandPass:
    break;

  case FT_HighShelve:
    break;

  case FT_LowShelve:
    break;

  case FT_Peak:
    break;
  }
}
*/

/*
// Second order low pass filter
void filterLP2(float f, Uint32 fs, Sint16* stream, int len) {
  static double x1, x2, y1, y2;
  double a1, a2, a3, b1, b2;
  double K, K2, denom;
  double x0, y0;

  K = tan(M_PI*f/fs);
  K2 = K*K;
  denom = 1 + sqrt(2)*K + K2;

  a1 = K2/denom;
  a2 = 2*a1;
  a3 = a1;
  b1 = a2-(2/denom);
  b2 = (1-sqrt(2)*K+K2)/denom;

  for (int i=0; i<(len/2); i++) {
    x0 = stream[i];
    y0 = a1*x0 + a2*x1 + a3*x2 - b1*y1 - b2*y2;
    x2 = x1;
    x1 = x0;
    y2 = y1;
    y1 = y0;
    if (y0>32767) y0 = 32767;
    if (y0<-32767) y0 = -32767;
    stream[i] = y0;
  }   
}
*/

#endif
