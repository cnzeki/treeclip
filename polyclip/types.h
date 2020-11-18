#pragma once
#include <cmath>
#include <ctime>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cfloat>
#include <exception>
#include "clipper/clipper.hpp"
#include "common/polygon.h"
#include "common/utilities.h"

using namespace std;
typedef std::vector< Point > Poly;
typedef std::vector< Poly > Polys;

#define and    &&
#define or     ||
#define not    !
#define equal  ==
