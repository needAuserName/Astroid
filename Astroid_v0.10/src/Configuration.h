/*    
 *    Configuration.h		
 *
 *    Copyright (C) 2014 University of Kentucky and
 *                       Yan Huang
 *
 *    Authors: Yan Huang
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef CONFIGURATION

#define CONFIGURATION

#define UNIX

#ifdef UNIX

#include <fstream>
#include <cstring>
#include <iostream>
#include <cmath>
#include <cstdlib> 
#include <ctime>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <queue>
#include <stack>
#include <sstream>
#include <assert.h>

#else

#include <stdio.h>
#include <fstream>
#include <string>
#include <iostream>
#include <math.h>
#include <stdlib.h> 
#include <time.h>
#include <vector>
#include <algorithm>
#include <queue>
#include <stack>
#include <sstream>
#include <assert.h>

#endif

using namespace std;

/************************************************************************/
/* GENERAL SETTINGS                                                     */
/************************************************************************/
//#define UNIX

//////////////////////////////////////////////////////////////////////////
//*******FOR DATA TYPE AND SCALE*******

//#define JUNCTIONONLY

#ifndef JUNCTIONONLY
#define DO_ESTIMATION
const int SUPPORT_VECTOR_SIZE = 1;//8;
#else
const int SUPPORT_VECTOR_SIZE = 1;
#endif

#define NORMALIZE_COVERAGE
#define TRANSCRIPTION_DIRECTION
#define FILTER_JUNCTION
#define FILTER_FRAGMENTS
const double JSD_RELIABILITY_MIN_COVERAGE = 10;
const int JSD_RELIABILITY_THRESHOLD = 5;
const double EXPRESSION_RELIABILITY_MIN_COVERAGE = 10;
//#define DO_MERGEALTERSITES
#define COVERAGE
#define COUNT_MAJOR_PATHS_ONLY
//#define COUTSTEPS


const long MAX_NUMBER = 2000000000;
const long default_dataset_num = 100;

#endif