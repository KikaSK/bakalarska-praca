#ifndef MESH_H
#define MESH_H

#include "assertm.h"
#include "edge.h"
#include "point.h"
#include "triangle.h"
#include "vector.h"
#include <ginac/ginac.h>
#include <iostream>
#include <vector>

using namespace GiNaC;
using std::endl, std::vector;

// co chcem robit
// pridavat trojuholniky, hrany, body
// vyberat a pridavat do active_edges
// mat normaly trojuholnikov ku ktorym viem pristupovat
// mat hrany, vediet povedat, ze tuto hranu maju tieto 1-2 trojuholniky
// mat body, vediet povedat ktore hrany vychadzaju z bodu
// vediet vypocitat vektor kolmy na hranu v rovine susedneho trojuholnika

class Mesh {
private:
public:
  Mesh(Triangle T);
  Mesh()
};

#endif
