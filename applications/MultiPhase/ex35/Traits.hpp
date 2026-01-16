#pragma once
#include <array>
#include <stdexcept>
#include <limits>

constexpr unsigned UMAX = std::numeric_limits<unsigned>::max(); // no static needed
// ============================================================
// Element type ids
// 0 Hex27, 1 Tet15, 2 Wedge21, 3 Quad9, 4 Tri7, 5 Line3
// ============================================================
enum ElType : unsigned {
  Hex27   = 0,
  Tet15   = 1,
  Wedge21 = 2,
  Quad9   = 3,
  Tri7    = 4,
  Line3   = 5
};

struct Range {
  unsigned start = 0;
  unsigned count = 0;
};

struct ElemTraits {
  ElType type{};
  unsigned dim      = 0;
  unsigned connSize = 0;

  unsigned nVert   = 0;
  unsigned nEdge   = 0;
  unsigned nFace   = 0;
  unsigned nCenter = 0;

  Range vert{};
  Range edge{};
  Range faceCenters{};
  Range center{};

  unsigned nSides = 0;

  static constexpr unsigned MaxSides     = 6;
  static constexpr unsigned MaxSideNodes = 9;

  std::array<unsigned, MaxSides> sideNodeCount{};
  std::array<std::array<unsigned, MaxSideNodes>, MaxSides> sideNodes{};
};

inline const ElemTraits& getTraits(unsigned elType) {
  static const std::array<ElemTraits, 6> T = []{
    std::array<ElemTraits, 6> a{};

    // ============================================================
    // HEX27
    //
    //         7------14-------6
    //        /|              /|
    //       / |             / |
    //     15  |   25      13  |
    //     /  19      22   /  18
    //    /    |          /    |
    //   4------12-------5     |
    //   | 23  |   26    | 21  |
    //   |     3------10-|-----2
    //   |    /          |    /
    //  16   /  20      17   /
    //   | 11      24    |  9
    //   | /             | /
    //   |/              |/
    //   0-------8-------1
    // ============================================================
    {
      ElemTraits t{};
      t.type = Hex27; t.dim = 3; t.connSize = 27;
      t.nVert=8; t.nEdge=12; t.nFace=6; t.nCenter=1;

      t.vert={0,8};
      t.edge={8,12};
      t.faceCenters={20,6};
      t.center={26,1};
      t.nSides=6;

      // f=0..3 sides, f=4 bottom, f=5 top
      t.sideNodeCount[0]=9; t.sideNodes[0]={ {0,1,5,4, 8,17,12,16,20} };
      t.sideNodeCount[1]=9; t.sideNodes[1]={ {1,2,6,5, 9,18,13,17,21} };
      t.sideNodeCount[2]=9; t.sideNodes[2]={ {2,3,7,6,10,19,14,18,22} };
      t.sideNodeCount[3]=9; t.sideNodes[3]={ {3,0,4,7,11,16,15,19,23} };
      t.sideNodeCount[4]=9; t.sideNodes[4]={ {0,3,2,1,11,10,9,8,24} };   // bottom (outward)
      t.sideNodeCount[5]=9; t.sideNodes[5]={ {4,5,6,7,12,13,14,15,25} }; // top

      a[Hex27]=t;
    }

    // ============================================================
    // TET15
    //
    //            3
    //           /|\
    //          / | \
    //         /  |  \
    //        9   |12 8
    //       /    |    \
    //      /   14|  11 \
    //     / 13   7      \
    //    2-------|5------1
    //     \      |      /
    //      \     |     /
    //       \    |10  /
    //        6   |   4
    //         \  |  /
    //          \ | /
    //           \|/
    //            0
    // ============================================================
    {
      ElemTraits t{};
      t.type = Tet15; t.dim = 3; t.connSize = 15;
      t.nVert=4; t.nEdge=6; t.nFace=4; t.nCenter=1;

      t.vert={0,4};
      t.edge={4,6};
      t.faceCenters={10,4};
      t.center={14,1};
      t.nSides=4;

      t.sideNodeCount[0]=7; t.sideNodes[0]={ {0,2,1, 6,5,4,10,0,0} };
      t.sideNodeCount[1]=7; t.sideNodes[1]={ {0,1,3, 4,8,7,11,0,0} };
      t.sideNodeCount[2]=7; t.sideNodes[2]={ {1,2,3, 5,9,8,12,0,0} };
      t.sideNodeCount[3]=7; t.sideNodes[3]={ {0,3,2, 7,9,6,13,0,0} };

      a[Tet15]=t;
    }

    // ============================================================
    // WEDGE21
    //
    //           5
    //          /|\
    //         / | \
    //        /  |  \
    //      11  19   10
    //      /   14    \
    //     /     |     \
    //    /      |      \
    //   3-------9-------4
    //   |  17  20  16   |
    //   |       2       |
    //   |      / \      |
    //   |     /   \     |
    //  12    / 15  \   13
    //   |   8       7   |
    //   |  /   18    \  |
    //   | /           \ |
    //   |/             \|
    //   0-------6-------1
    // ============================================================
    {
      ElemTraits t{};
      t.type = Wedge21; t.dim = 3; t.connSize = 21;
      t.nVert=6; t.nEdge=9; t.nFace=5; t.nCenter=1;

      t.vert={0,6};
      t.edge={6,9};
      t.faceCenters={15,5};
      t.center={20,1};
      t.nSides=5;

      t.sideNodeCount[0]=9; t.sideNodes[0]={ {0,1,4,3,6,13,9,12,15} };
      t.sideNodeCount[1]=9; t.sideNodes[1]={ {1,2,5,4,7,14,10,13,16} };
      t.sideNodeCount[2]=9; t.sideNodes[2]={ {2,0,3,5,8,12,11,14,17} };
      t.sideNodeCount[3]=7; t.sideNodes[3]={ {0,2,1,8,7,6,18,0,0} }; // bottom
      t.sideNodeCount[4]=7; t.sideNodes[4]={ {3,4,5,9,10,11,19,0,0} }; // top

      a[Wedge21]=t;
    }

    // ============================================================
    // QUAD9
    //
    //      3-----6-----2
    //      |           |
    //      |           |
    //      7     8     5
    //      |           |
    //      |           |
    //      0-----4-----1
    // ============================================================
    {
      ElemTraits t{};
      t.type = Quad9; t.dim = 2; t.connSize = 9;
      t.nVert=4; t.nEdge=4; t.nFace=0; t.nCenter=1;

      t.vert={0,4};
      t.edge={4,4};
      t.center={8,1};
      t.nSides=4;

      t.sideNodeCount[0]=3; t.sideNodes[0]={ {0,1,4,0,0,0,0,0,0} };
      t.sideNodeCount[1]=3; t.sideNodes[1]={ {1,2,5,0,0,0,0,0,0} };
      t.sideNodeCount[2]=3; t.sideNodes[2]={ {2,3,6,0,0,0,0,0,0} };
      t.sideNodeCount[3]=3; t.sideNodes[3]={ {3,0,7,0,0,0,0,0,0} };

      a[Quad9]=t;
    }

    // ============================================================
    // TRI7
    //
    //      2
    //      | \
    //      |   \
    //      5     4
    //      |   6   \
    //      |         \
    //      0-----3----1
    // ============================================================
    {
      ElemTraits t{};
      t.type = Tri7; t.dim = 2; t.connSize = 7;
      t.nVert=3; t.nEdge=3; t.nFace=0; t.nCenter=1;

      t.vert={0,3};
      t.edge={3,3};
      t.center={6,1};
      t.nSides=3;

      t.sideNodeCount[0]=3; t.sideNodes[0]={ {0,1,3,0,0,0,0,0,0} };
      t.sideNodeCount[1]=3; t.sideNodes[1]={ {1,2,4,0,0,0,0,0,0} };
      t.sideNodeCount[2]=3; t.sideNodes[2]={ {2,0,5,0,0,0,0,0,0} };

      a[Tri7]=t;
    }

    // ============================================================
    // LINE3
    //
    //	0-----2-----1
    // ============================================================
    {
      ElemTraits t{};
      t.type = Line3; t.dim = 1; t.connSize = 3;
      t.nVert=2; t.nEdge=0; t.nFace=0; t.nCenter=1;

      t.vert={0,2};
      t.center={2,1};
      t.nSides=2;

      t.sideNodeCount[0]=1; t.sideNodes[0]={ {0,0,0,0,0,0,0,0,0} };
      t.sideNodeCount[1]=1; t.sideNodes[1]={ {1,0,0,0,0,0,0,0,0} };

      a[Line3]=t;
    }

    return a;
  }();

  if (elType >= T.size())
    throw std::runtime_error("getTraits: invalid elType");

  return T[elType];
}

inline const ElemTraits& getTraits(ElType t) { return getTraits(static_cast<unsigned>(t)); }

