/*
 *  robust.h
 *  robust
 *
 *  Created by Kjell Konis on 24/04/2006.
 *  Copyright 2006. All rights reserved.
 *
 */

#ifdef USING_R
  typedef double Sfloat;
  typedef int Sint;
  #define SINT_MAX INT_MAX
  #define SINT_MIN INT_MIN
#else
  typedef double Sfloat;
  typedef long Sint;
  #define SINT_MAX LONG_MAX
  #define SINT_MIN LONG_MIN
#endif

