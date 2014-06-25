#pragma once

//2q9 it is
static const unsigned int stencilSize = 9;

//stencil indices; keep the ordering intact!!
#define C  0
#define N  1
#define NE 2
#define E  3
#define SE 4
#define S  5
#define SW 6
#define W  7
#define NW 8

//stencil neighbours
static int stenNbs[2][stencilSize] = {
// X: C  N  NE E  SE  S  SW   W  NW
     {0, 0, 1, 1, 1,  0, -1, -1,-1},
// Y: C  N  NE E  SE  S  SW   W  NW
     {0, 1, 1, 0,-1, -1, -1,  0, 1}};

static const double weigths[stencilSize] = 
//   C      N     NE    E      SE      S      SW      W      NW
  {4./9.,1./9.,1./36.,1./9., 1./36., 1./9., 1./36., 1./9., 1./36.};

//Directions for velocity-vector
#define X 0
#define Y 1
