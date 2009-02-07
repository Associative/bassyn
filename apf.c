//
//  B A S S Y N   -   B A S I C   S Y N T H E S I Z E R
//
//  Modular synthesizer with various generators and filters
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  apf.c
//

#include <math.h>
#include <stdlib.h>
#include "apf.h"



// global variables
unsigned long samplerate;
unsigned long long samplecount;

// module parameter strings
const char *valueStrings[]        = {NULL}; 
const char *sumStrings[]          = {"Input A", "Input B"};
const char *gainStrings[]         = {"Input", "Gain Factor"};
const char *gSineStrings[]        = {"Frequency", "Phaseshift", "Amplitude"};
const char *gTriangleStrings[]    = {"Frequency", "Phaseshift", "Amplitude", "Skew ratio"};
const char *gSquareStrings[]      = {"Frequency", "Phaseshift", "Duty cycle"};
const char *f2HighPassStrings[]   = {"Input", "Bypass", "Cutoff frequency"};
const char *f2LowPassStrings[]    = {"Input", "Bypass", "Cutoff frequency"};
const char *f2BandPassStrings[]   = {"Input", "Bypass", "Center frequency", "Q-factor"};
const char *f2HighShelveStrings[] = {"Input", "Bypass", "Corner frequency", "Boost factor"};
const char *f2LowShelveStrings[]  = {"Input", "Bypass", "Corner frequency", "Boost factor"};
const char *f2PeakStrings[]       = {"Input", "Bypass", "Center frequency", "Q-factor"};

// module descriptors
const struct moduleDesc moduleInfo[] = {

  // MT_Value
  { "Constant",
    0, valueStrings,
    0, valueProcess
  },

  // MT_Sum
  { "Sum",
    2, sumStrings,
    0, sumProcess
  },

  // MT_Gain
  { "Gain",
    2, gainStrings,
    0, gainProcess
  },

  // MT_G_Sine
  { "Sine generator",
    3, gSineStrings,
    0, gSineProcess
  },

  // MT_G_Triangle
  { "Triangle generator",
    4, gTriangleStrings,
    0, gTriangleProcess
  },

  // MT_G_Square
  { "Square generator",
    3, gSquareStrings,
    0, gSquareProcess
  },

  // MT_G_F2_HighPass
  { "High-pass filter 2",
    3, f2HighPassStrings,
    sizeof(struct filterMemory2), f2HighPassProcess
  },

  // MT_F2_LowPass
  { "Low-pass filter 2",
    3, f2LowPassStrings,
    sizeof(struct filterMemory2), f2LowPassProcess
  },

  // MT_F2_BandPass
  { "Band-pass filter 2",
    4, f2BandPassStrings,
    sizeof(struct filterMemory2), f2BandPassProcess
  },

  // MT_F2_HighShelve
  { "High shelving filter 2",
    4, f2HighShelveStrings,
    sizeof(struct filterMemory2), f2HighShelveProcess
  },

  // MT_F2_LowShelve
  { "Low shelving filter 2",
    4, f2LowShelveStrings,
    sizeof(struct filterMemory2), f2LowShelveProcess
  },

  // MT_F2_Peak
  { "Peaking filter 2",
    4, f2PeakStrings,
    sizeof(struct filterMemory2), f2PeakProcess
  }
  
};

//  G E N E R A L   F U N C T I O N S
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// calculate the level in the processing hierarchy
int calcLevel(struct module *m) {
  struct module *n;
  int l=0;

  for (int i=0; i<m->info->paramCount; i++) {
    n = m->param[i];
    if ((n != NULL) && (n->level+1 > l)) l = n->level+1;
  }
 
  return l;
}

// sort the system according to level
struct module* sortModules(struct module *ml) {
  struct module *m, *mr;

  if (ml->link == NULL) return ml;

  mr = sortModules(ml->link);
  if (ml->level < mr->level) {
    ml->link = mr;
    return ml;
  }

  m = mr;
  while ((m->link != NULL) && (m->level > m->link->level)) {
    m = m->link; 
  }

  ml->link = m->link;
  m->link = ml;

  return mr;
}

// process the system
double processModules(struct module *ml) {
  struct module *m = ml;

  while (m->link != NULL) {
    m->info->process(m);
    m = m->link;
  }
  m->info->process(m);

  return m->value;
}

// create a new module
struct module* createModule(int t, double v, struct module *l) {
  struct module *m;

  m = (struct module *) malloc(sizeof(struct module));
  m->info = &moduleInfo[t];
  m->level = 0;
  m->value = v;
  m->param = (struct module **) malloc(moduleInfo[t].paramCount * sizeof(struct module *));
  m->data = malloc(moduleInfo[t].dataSize);
  m->link = l;

  return m;
}

// delete a module (returns the next module in the chain)
struct module* freeModule(struct module *m) {
  struct module *l = m->link; 

  free(m->param);
  free(m->data);
  free(m);

  return l;
}



//  M O D U L E   S P E C I F I C   F U N C T I O N S
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

void valueProcess(struct module *m) {
  ;
}

void sumProcess(struct module *m) {
  m->value = m->param[0]->value + m->param[1]->value;
}

void gainProcess(struct module *m) {
  m->value = m->param[0]->value * m->param[1]->value;
}

void gSineProcess(struct module *m) {
  double x, y;
  x = getPeriodSample(m->param[0]->value, m->param[1]->value);
  y = clipIn(m->param[2]->value, 1.0f) * sin(x);
  m->value = clipEx(y, 1.0f);
}

void gTriangleProcess(struct module *m) {
  m->value = 0;
}

void gSquareProcess(struct module *m) {
  m->value = 0;
}

void f2HighPassProcess(struct module *m) {
  m->value = 0;
}

void f2LowPassProcess(struct module *m) {
  double K;
  m->value = 0;
}

void f2BandPassProcess(struct module *m) {
  m->value = 0;
}

void f2HighShelveProcess(struct module *m) {
  m->value = 0;
}

void f2LowShelveProcess(struct module *m) {
  m->value = 0;
}

void f2PeakProcess(struct module *m) {
  m->value = 0;
}


//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

// min <= x <= max
double clipEx(double x, double max) {
  if (x>max) return max;
  if (x<-max) return -max;
  return x;
}

// x <= -max || max >= x
double clipIn(double x, double max) {
  if (x>0 && x<max) return max;
  if (x>-max && x<0) return -max;
  return x;
}

// get the samplenumber within a certain period
double getPeriodSample(double f, double ph) {
  return fmod((2*M_PI*f*samplecount)/samplerate + ph, 2*M_PI);
}
