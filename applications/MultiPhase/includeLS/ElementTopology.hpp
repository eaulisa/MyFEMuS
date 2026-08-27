#pragma once

#include <cstddef>
#include <stdexcept>

namespace ex40_topo {

enum SideShape : unsigned { Quad, Tri, Line, Point };

struct ConnView {
  const unsigned *p;
  std::size_t n;

  const unsigned *begin() const { return p; }
  const unsigned *end() const { return p + n; }
  const unsigned &operator[](std::size_t i) const { return p[i]; }
};

struct SideInfo {
  SideShape shape;
  ConnView nodes;    // tutti i nodi della side
  ConnView vertices; // solo vertici della side
};

struct EdgeInfo {
  ConnView nodes;    // nodi sull’edge (incl. mid)
  ConnView vertices; // i 2 vertici dell’edge
};

struct ElementTopology {
  unsigned typeId;
  unsigned dim;
  unsigned nNodes;
  unsigned nVertices;

  unsigned nSides; // facce (3D) / edges (2D) / estremi (1D)
  const SideInfo *sides;

  unsigned nEdges; // edges geometrici (utile soprattutto 3D)
  const EdgeInfo *edges;
};

// ------------------------------------------------------------
// Tabelle (internal linkage: static const in header = ok)
// ------------------------------------------------------------

// Hex27 (type id = 0)
static const unsigned HEX_FACE_NODES[6][9] = {
    {0, 1, 5, 4, 8, 17, 12, 16, 20},  {1, 2, 6, 5, 9, 18, 13, 17, 21},
    {2, 3, 7, 6, 10, 19, 14, 18, 22}, {3, 0, 4, 7, 11, 16, 15, 19, 23},
    {0, 3, 2, 1, 11, 10, 9, 8, 24},   {4, 5, 6, 7, 12, 13, 14, 15, 25}};

static const unsigned HEX_FACE_VERTS[6][4] = {{0, 1, 5, 4}, {1, 2, 6, 5},
                                              {2, 3, 7, 6}, {3, 0, 4, 7},
                                              {0, 3, 2, 1}, {4, 5, 6, 7}};

static const unsigned HEX_EDGE_NODES[12][3] = {
    {0, 1, 8},  {1, 2, 9},  {2, 3, 10}, {3, 0, 11}, {4, 5, 12}, {5, 6, 13},
    {6, 7, 14}, {7, 4, 15}, {0, 4, 16}, {1, 5, 17}, {2, 6, 18}, {3, 7, 19}};

static const unsigned HEX_EDGE_VERTS[12][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0},
                                               {4, 5}, {5, 6}, {6, 7}, {7, 4},
                                               {0, 4}, {1, 5}, {2, 6}, {3, 7}};

static const SideInfo HEX_SIDES[6] = {
    {Quad, {HEX_FACE_NODES[0], 9}, {HEX_FACE_VERTS[0], 4}},
    {Quad, {HEX_FACE_NODES[1], 9}, {HEX_FACE_VERTS[1], 4}},
    {Quad, {HEX_FACE_NODES[2], 9}, {HEX_FACE_VERTS[2], 4}},
    {Quad, {HEX_FACE_NODES[3], 9}, {HEX_FACE_VERTS[3], 4}},
    {Quad, {HEX_FACE_NODES[4], 9}, {HEX_FACE_VERTS[4], 4}},
    {Quad, {HEX_FACE_NODES[5], 9}, {HEX_FACE_VERTS[5], 4}}};

static const EdgeInfo HEX_EDGES[12] = {
    {{HEX_EDGE_NODES[0], 3}, {HEX_EDGE_VERTS[0], 2}},
    {{HEX_EDGE_NODES[1], 3}, {HEX_EDGE_VERTS[1], 2}},
    {{HEX_EDGE_NODES[2], 3}, {HEX_EDGE_VERTS[2], 2}},
    {{HEX_EDGE_NODES[3], 3}, {HEX_EDGE_VERTS[3], 2}},
    {{HEX_EDGE_NODES[4], 3}, {HEX_EDGE_VERTS[4], 2}},
    {{HEX_EDGE_NODES[5], 3}, {HEX_EDGE_VERTS[5], 2}},
    {{HEX_EDGE_NODES[6], 3}, {HEX_EDGE_VERTS[6], 2}},
    {{HEX_EDGE_NODES[7], 3}, {HEX_EDGE_VERTS[7], 2}},
    {{HEX_EDGE_NODES[8], 3}, {HEX_EDGE_VERTS[8], 2}},
    {{HEX_EDGE_NODES[9], 3}, {HEX_EDGE_VERTS[9], 2}},
    {{HEX_EDGE_NODES[10], 3}, {HEX_EDGE_VERTS[10], 2}},
    {{HEX_EDGE_NODES[11], 3}, {HEX_EDGE_VERTS[11], 2}}};

// Tet15 (type id = 1)
static const unsigned TET_FACE_NODES[4][7] = {{0, 2, 1, 6, 5, 4, 10},
                                              {0, 1, 3, 4, 8, 7, 11},
                                              {1, 2, 3, 5, 9, 8, 12},
                                              {2, 0, 3, 6, 7, 9, 13}};

static const unsigned TET_FACE_VERTS[4][3] = {
    {0, 2, 1}, {0, 1, 3}, {1, 2, 3}, {2, 0, 3}};

static const unsigned TET_EDGE_NODES[6][3] = {{0, 1, 4}, {1, 2, 5}, {0, 2, 6},
                                              {0, 3, 7}, {1, 3, 8}, {2, 3, 9}};

static const unsigned TET_EDGE_VERTS[6][2] = {{0, 1}, {1, 2}, {0, 2},
                                              {0, 3}, {1, 3}, {2, 3}};

static const SideInfo TET_SIDES[4] = {
    {Tri, {TET_FACE_NODES[0], 7}, {TET_FACE_VERTS[0], 3}},
    {Tri, {TET_FACE_NODES[1], 7}, {TET_FACE_VERTS[1], 3}},
    {Tri, {TET_FACE_NODES[2], 7}, {TET_FACE_VERTS[2], 3}},
    {Tri, {TET_FACE_NODES[3], 7}, {TET_FACE_VERTS[3], 3}}};

static const EdgeInfo TET_EDGES[6] = {
    {{TET_EDGE_NODES[0], 3}, {TET_EDGE_VERTS[0], 2}},
    {{TET_EDGE_NODES[1], 3}, {TET_EDGE_VERTS[1], 2}},
    {{TET_EDGE_NODES[2], 3}, {TET_EDGE_VERTS[2], 2}},
    {{TET_EDGE_NODES[3], 3}, {TET_EDGE_VERTS[3], 2}},
    {{TET_EDGE_NODES[4], 3}, {TET_EDGE_VERTS[4], 2}},
    {{TET_EDGE_NODES[5], 3}, {TET_EDGE_VERTS[5], 2}}};

// Wedge21 (type id = 2)
static const unsigned WDG_QUAD_FACE_NODES[3][9] = {
    {0, 1, 4, 3, 6, 13, 9, 12, 15},
    {1, 2, 5, 4, 7, 14, 10, 13, 16},
    {2, 0, 3, 5, 8, 12, 11, 14, 17}};
static const unsigned WDG_QUAD_FACE_VERTS[3][4] = {
    {0, 1, 4, 3}, {1, 2, 5, 4}, {2, 0, 3, 5}};

static const unsigned WDG_TRI_FACE_NODES[2][7] = {{0, 2, 1, 8, 7, 6, 18},
                                                  {3, 4, 5, 9, 10, 11, 19}};
static const unsigned WDG_TRI_FACE_VERTS[2][3] = {{0, 2, 1}, {3, 4, 5}};

static const unsigned WDG_EDGE_NODES[9][3] = {
    {0, 1, 6},  {1, 2, 7},  {2, 0, 8},  {3, 4, 9}, {4, 5, 10},
    {5, 3, 11}, {0, 3, 12}, {1, 4, 13}, {2, 5, 14}};
static const unsigned WDG_EDGE_VERTS[9][2] = {
    {0, 1}, {1, 2}, {2, 0}, {3, 4}, {4, 5}, {5, 3}, {0, 3}, {1, 4}, {2, 5}};

static const SideInfo WDG_SIDES[5] = {
    {Quad, {WDG_QUAD_FACE_NODES[0], 9}, {WDG_QUAD_FACE_VERTS[0], 4}},
    {Quad, {WDG_QUAD_FACE_NODES[1], 9}, {WDG_QUAD_FACE_VERTS[1], 4}},
    {Quad, {WDG_QUAD_FACE_NODES[2], 9}, {WDG_QUAD_FACE_VERTS[2], 4}},
    {Tri, {WDG_TRI_FACE_NODES[0], 7}, {WDG_TRI_FACE_VERTS[0], 3}},
    {Tri, {WDG_TRI_FACE_NODES[1], 7}, {WDG_TRI_FACE_VERTS[1], 3}}};

static const EdgeInfo WDG_EDGES[9] = {
    {{WDG_EDGE_NODES[0], 3}, {WDG_EDGE_VERTS[0], 2}},
    {{WDG_EDGE_NODES[1], 3}, {WDG_EDGE_VERTS[1], 2}},
    {{WDG_EDGE_NODES[2], 3}, {WDG_EDGE_VERTS[2], 2}},
    {{WDG_EDGE_NODES[3], 3}, {WDG_EDGE_VERTS[3], 2}},
    {{WDG_EDGE_NODES[4], 3}, {WDG_EDGE_VERTS[4], 2}},
    {{WDG_EDGE_NODES[5], 3}, {WDG_EDGE_VERTS[5], 2}},
    {{WDG_EDGE_NODES[6], 3}, {WDG_EDGE_VERTS[6], 2}},
    {{WDG_EDGE_NODES[7], 3}, {WDG_EDGE_VERTS[7], 2}},
    {{WDG_EDGE_NODES[8], 3}, {WDG_EDGE_VERTS[8], 2}}};

// Quad9 (type id = 3)
static const unsigned Q9_EDGE_NODES[4][3] = {
    {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};
static const unsigned Q9_EDGE_VERTS[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
static const SideInfo Q9_SIDES[4] = {
    {Line, {Q9_EDGE_NODES[0], 3}, {Q9_EDGE_VERTS[0], 2}},
    {Line, {Q9_EDGE_NODES[1], 3}, {Q9_EDGE_VERTS[1], 2}},
    {Line, {Q9_EDGE_NODES[2], 3}, {Q9_EDGE_VERTS[2], 2}},
    {Line, {Q9_EDGE_NODES[3], 3}, {Q9_EDGE_VERTS[3], 2}}};
static const EdgeInfo Q9_EDGES[4] = {
    {{Q9_EDGE_NODES[0], 3}, {Q9_EDGE_VERTS[0], 2}},
    {{Q9_EDGE_NODES[1], 3}, {Q9_EDGE_VERTS[1], 2}},
    {{Q9_EDGE_NODES[2], 3}, {Q9_EDGE_VERTS[2], 2}},
    {{Q9_EDGE_NODES[3], 3}, {Q9_EDGE_VERTS[3], 2}}};

// Tri7 (type id = 4)
static const unsigned T7_EDGE_NODES[3][3] = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
static const unsigned T7_EDGE_VERTS[3][2] = {{0, 1}, {1, 2}, {2, 0}};
static const SideInfo T7_SIDES[3] = {
    {Line, {T7_EDGE_NODES[0], 3}, {T7_EDGE_VERTS[0], 2}},
    {Line, {T7_EDGE_NODES[1], 3}, {T7_EDGE_VERTS[1], 2}},
    {Line, {T7_EDGE_NODES[2], 3}, {T7_EDGE_VERTS[2], 2}}};
static const EdgeInfo T7_EDGES[3] = {
    {{T7_EDGE_NODES[0], 3}, {T7_EDGE_VERTS[0], 2}},
    {{T7_EDGE_NODES[1], 3}, {T7_EDGE_VERTS[1], 2}},
    {{T7_EDGE_NODES[2], 3}, {T7_EDGE_VERTS[2], 2}}};

// Line3 (type id = 5)
static const unsigned L3_SIDE_NODES[2][1] = {{0}, {1}};
static const SideInfo L3_SIDES[2] = {
    {Point, {L3_SIDE_NODES[0], 1}, {L3_SIDE_NODES[0], 1}},
    {Point, {L3_SIDE_NODES[1], 1}, {L3_SIDE_NODES[1], 1}}};

static const unsigned L3_EDGE_NODES[1][3] = {{0, 1, 2}};
static const unsigned L3_EDGE_VERTS[1][2] = {{0, 1}};
static const EdgeInfo L3_EDGES[1] = {
    {{L3_EDGE_NODES[0], 3}, {L3_EDGE_VERTS[0], 2}}};

// ------------------------------------------------------------
// Public API
// ------------------------------------------------------------
inline const ElementTopology &topologyFromTypeId(unsigned t) {
  static const ElementTopology HEX = {0, 3, 27, 8, 6, HEX_SIDES, 12, HEX_EDGES};
  static const ElementTopology TET = {1, 3, 15, 4, 4, TET_SIDES, 6, TET_EDGES};
  static const ElementTopology WDG = {2, 3, 21, 6, 5, WDG_SIDES, 9, WDG_EDGES};
  static const ElementTopology Q9 = {3, 2, 9, 4, 4, Q9_SIDES, 4, Q9_EDGES};
  static const ElementTopology T7 = {4, 2, 7, 3, 3, T7_SIDES, 3, T7_EDGES};
  static const ElementTopology L3 = {5, 1, 3, 2, 2, L3_SIDES, 1, L3_EDGES};

  switch (t) {
  case 0:
    return HEX;
  case 1:
    return TET;
  case 2:
    return WDG;
  case 3:
    return Q9;
  case 4:
    return T7;
  case 5:
    return L3;
  default:
    throw std::runtime_error(
        "topologyFromTypeId: unknown element type id (expected 0..5)");
  }
}

} // namespace ex40_topo
