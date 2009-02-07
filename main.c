//
//  B A S S Y N   -   B A S I C   S Y N T H E S I Z E R
//
//  Modular synthesizer with various generators and filters
//
//  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
//
//  main.c
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef __APPLE__
  #include "glut.h"
#else
  #include "GL/glut.h"
#endif

#include "SDL.h"
#include "SDL_audio.h"
#include "apf.h"

Uint32 running;
double wave[2048];
float win_w, win_h;

// synthesizer entry and endpoint pointers
struct module *synthRoot;
struct module *synthOutput;

// program specific functions
void initGL(void);
void handleEvents(void);
void procAudio(void *userdata, Sint16 *stream, int len);
void drawFrame(float aDiv, float tDiv, float trig, char pf); 


// program entry
int main(int argc, char *argv[])
{
  struct module *tp;
  SDL_AudioSpec *desired, *obtained;

  win_w = 600;
  win_h = 400;

  if (SDL_Init(SDL_INIT_VIDEO | SDL_INIT_AUDIO) < 0) {
    fprintf(stderr, "Couldn't initialize SDL: %s\n", SDL_GetError());
    exit(1);
  }
  SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER, 1);
  SDL_SetVideoMode(win_w, win_h, 0, SDL_OPENGL);
  atexit(SDL_Quit);

  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

  samplecount = 0;
  samplerate = 44100; 
  
  desired = malloc(sizeof(SDL_AudioSpec));
  obtained = malloc(sizeof(SDL_AudioSpec));

  desired->freq = samplerate;
  desired->format = AUDIO_S16SYS;
  desired->channels = 1;
  desired->samples = 2048;
  desired->callback = procAudio;
  desired->userdata = wave;

  if (SDL_OpenAudio(desired, obtained) < 0) {
    fprintf(stderr, "Couldn't open audio: %s\n", SDL_GetError());
    exit(1);
  }

  memset(wave, 0, 2048*sizeof(double));

////////////////////////////
  synthRoot = createModule(MT_Value, 1e3f, NULL);
  synthRoot->level = calcLevel(synthRoot); 
  tp = synthRoot;

  tp->link = createModule(MT_Value, 0, NULL);
  tp = tp->link;
  tp->level = calcLevel(tp);

  tp->link = createModule(MT_Value, 1.0f, NULL);
  tp = tp->link;
  tp->level = calcLevel(tp);

  tp->link = createModule(MT_G_Sine, 0, NULL);
  tp = tp->link;
  tp->param[0] = synthRoot;
  tp->param[1] = synthRoot->link;
  tp->param[2] = synthRoot->link->link;
  tp->level = calcLevel(tp);

  synthRoot = sortModules(synthRoot);

  tp = synthRoot;
  while (tp != NULL) {
    printf("%2d %s value: %f\n", tp->level, tp->info->name, tp->value);
    tp = tp->link;
  }  
  //printf("result: %f\n", processModules(synthRoot));
////////////////////////////

  // Play!
  SDL_PauseAudio(0);
  
  running = 1;
  while (running) {
    handleEvents();
    drawFrame(1.0f, 5e-3, 0.0f, 0);
  }

  return(0);
}


// OpenGL initialization
void initGL(void) {
  glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  glClearDepth(1.0f);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glViewport(0,0,500,300);
  glMatrixMode(GL_MODELVIEW);
}


// eventhandler for keypresses and the like
void handleEvents(void) {
  SDL_Event e;

  while (SDL_PollEvent(&e)) {
    if (e.type == SDL_QUIT)
      running--;

    if (e.type == SDL_KEYDOWN) {
      switch (e.key.keysym.sym) {
        case SDLK_q:
        case SDLK_ESCAPE:
          running--;
          break;
        default:
          break;
      }
    }
  }
} 


// draw video frame
void drawFrame(float aDiv, float tDiv, float trig, char pf) {
  float x, dx, y0, y1;
  int i, j;

  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glLoadIdentity();

  glLineWidth(1.0f);

  glScalef(1.0f, 0.75f, 1.0f);

  glBegin(GL_LINES);
    glColor3f(0.0f, 1.0f, 0.0f);
    glVertex2f(-1.0,0.0f);
    glVertex2f(1.0,0.0f);
    glVertex2f(0.0f,-0.5f);
    glVertex2f(0.0f,0.5f);
    glVertex2f(-1.0f, 1.0f);
    glVertex2f(1.0f, 1.0f);
    glVertex2f(-1.0f, -1.0f);
    glVertex2f(1.0f, -1.0f);
  glEnd();

  glColor3f(1.0f,0.0f,0.0f);

  glBegin(GL_LINE_STRIP);
  dx = samplerate*tDiv;
  y1 = wave[0];
  for (i=0, j=-1; i<2048 && j<dx; i++) {
    y0 = wave[i];
    if ((j<0) && (i>dx || pf ? (y1<trig && y0>=trig) : (y1>trig && y0<=trig))) j=0;
    if (j>-1) {
      x = j*2.0f/dx - 1.0f;
      glVertex2f(x,y0*aDiv);
      j++;
    }
    y1 = y0;
  }
  glEnd();

  SDL_GL_SwapBuffers();
}

// audio processing callback
void procAudio(void *userdata, Sint16 *stream, int len) {
  double *wave = userdata;

  for (int i=0; i<(len/2); i++) {
    wave[i] = processModules(synthRoot);
    stream[i] = wave[i]*32767;
    ++samplecount;
  }
}
