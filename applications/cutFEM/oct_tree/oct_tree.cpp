#include <iostream>
#include <vector>
#include <cstdlib>
#include <ctime>

// Define a 3D point
struct Point3D {
    float x, y, z;

    Point3D(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
};

// Define an Octree Node
class OctreeNode {
public:
    Point3D minBounds, maxBounds;
    bool isLeaf;
    std::vector<OctreeNode*> children;

    OctreeNode(const Point3D& _minBounds, const Point3D& _maxBounds)
        : minBounds(_minBounds), maxBounds(_maxBounds), isLeaf(true) {}

    // Randomly subdivide the node into 8 equal child nodes up to a specified depth
    void subdivide(int maxDepth, int currentDepth = 0) {
        if (currentDepth >= maxDepth || !isLeaf) {
            return;
        }

        isLeaf = false;
        float midX = (minBounds.x + maxBounds.x) / 2.0f;
        float midY = (minBounds.y + maxBounds.y) / 2.0f;
        float midZ = (minBounds.z + maxBounds.z) / 2.0f;

        children.push_back(new OctreeNode(minBounds, {midX, midY, midZ}));
        children.push_back(new OctreeNode({midX, minBounds.y, minBounds.z}, {maxBounds.x, midY, midZ}));
        children.push_back(new OctreeNode({minBounds.x, midY, minBounds.z}, {midX, maxBounds.y, midZ}));
        children.push_back(new OctreeNode({midX, midY, minBounds.z}, {maxBounds.x, maxBounds.y, midZ}));
        children.push_back(new OctreeNode({minBounds.x, minBounds.y, midZ}, {midX, midY, maxBounds.z}));
        children.push_back(new OctreeNode({midX, minBounds.y, midZ}, {maxBounds.x, midY, maxBounds.z}));
        children.push_back(new OctreeNode({minBounds.x, midY, midZ}, {midX, maxBounds.y, maxBounds.z}));
        children.push_back(new OctreeNode({midX, midY, midZ}, maxBounds));
/*
        for (OctreeNode* child : children) {
            child->subdivide(maxDepth, currentDepth + 1);
        }*/
    }

    // Recursive random subdivision with a maximum depth
    void randomSubdivide(int maxDepth, int currentDepth = 0) {
        if (currentDepth >= maxDepth || !isLeaf) {
            return;
        }

        // Randomly decide whether to subdivide
        if (rand() % 2 == 1) {
            subdivide(maxDepth, currentDepth);
            for (OctreeNode* child : children) {
                child->randomSubdivide(maxDepth, currentDepth + 1);
            }
        }
    }

    // Recursive search for the smallest sub-cube containing a given point
    OctreeNode* search(const Point3D& point) {
        if (isLeaf) {
            return this;
        }

        for (OctreeNode* child : children) {
            if (point.x >= child->minBounds.x && point.x <= child->maxBounds.x &&
                point.y >= child->minBounds.y && point.y <= child->maxBounds.y &&
                point.z >= child->minBounds.z && point.z <= child->maxBounds.z) {
                return child->search(point);
            }
        }

        return nullptr; // Should not reach here under normal circumstances
    }
};

// Function to print the octree structure
void printOctreeStructure(OctreeNode* node, int depth = 0) {
    for (int i = 0; i < depth; ++i) {
        std::cout << "  ";
    }

    std::cout << "Node Bounds: ";
    std::cout << "(" << node->minBounds.x << ", " << node->minBounds.y << ", " << node->minBounds.z << ") to ";
    std::cout << "(" << node->maxBounds.x << ", " << node->maxBounds.y << ", " << node->maxBounds.z << ")";

    if (node->isLeaf) {
        std::cout << " [Leaf]\n";
    } else {
        std::cout << " [Non-Leaf]\n";
        for (OctreeNode* child : node->children) {
            printOctreeStructure(child, depth + 1);
        }
    }
}

int main() {
    // Create the root node representing the unit cube
    OctreeNode root({0.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f});

    // Specify the maximum depth of the Octree
    int maxDepth = 4;

    // Randomly subdivide the root node and its descendants
    srand(time(0));  // Seed for random number generation
    root.randomSubdivide(maxDepth);

    // Print the octree structure
    std::cout << "Octree Structure:\n";
    printOctreeStructure(&root);

    // Example: search for the smallest sub-cube containing a given point
    Point3D searchPoint(0.7f, 0.3f, 0.5f);
    OctreeNode* result = root.search(searchPoint);
    // Display the result
    std::cout << "\nSearch Point: (" << searchPoint.x << ", " << searchPoint.y << ", " << searchPoint.z << ")\n";
    std::cout << "Smallest Sub-cube Bounds: ";
    std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
    std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";

    return 0;
}















// #include <iostream>
// #include <vector>
// #include <cstdlib>
// #include <ctime>
//
// // Define a 3D point
// struct Point3D {
//     float x, y, z;
//
//     Point3D(float _x, float _y, float _z) : x(_x), y(_y), z(_z) {}
// };
//
// // Define an Octree Node
// class OctreeNode {
// public:
//     Point3D minBounds, maxBounds;
//     bool isLeaf;
//     std::vector<OctreeNode*> children;
//
//     OctreeNode(const Point3D& _minBounds, const Point3D& _maxBounds)
//         : minBounds(_minBounds), maxBounds(_maxBounds), isLeaf(true) {}
//
//     // Randomly subdivide the node into 8 equal child nodes
//     void subdivide() {
//         isLeaf = false;
//         float midX = (minBounds.x + maxBounds.x) / 2.0f;
//         float midY = (minBounds.y + maxBounds.y) / 2.0f;
//         float midZ = (minBounds.z + maxBounds.z) / 2.0f;
//
//         children.push_back(new OctreeNode(minBounds, {midX, midY, midZ}));
//         children.push_back(new OctreeNode({midX, minBounds.y, minBounds.z}, {maxBounds.x, midY, midZ}));
//         children.push_back(new OctreeNode({minBounds.x, midY, minBounds.z}, {midX, maxBounds.y, midZ}));
//         children.push_back(new OctreeNode({midX, midY, minBounds.z}, {maxBounds.x, maxBounds.y, midZ}));
//         children.push_back(new OctreeNode({minBounds.x, minBounds.y, midZ}, {midX, midY, maxBounds.z}));
//         children.push_back(new OctreeNode({midX, minBounds.y, midZ}, {maxBounds.x, midY, maxBounds.z}));
//         children.push_back(new OctreeNode({minBounds.x, midY, midZ}, {midX, maxBounds.y, maxBounds.z}));
//         children.push_back(new OctreeNode({midX, midY, midZ}, maxBounds));
//     }
//
//     // Recursive search for the smallest sub-cube containing a given point
//     OctreeNode* search(const Point3D& point) {
//         if (isLeaf) {
//             return this;
//         }
//
//         for (OctreeNode* child : children) {
//             if (point.x >= child->minBounds.x && point.x <= child->maxBounds.x &&
//                 point.y >= child->minBounds.y && point.y <= child->maxBounds.y &&
//                 point.z >= child->minBounds.z && point.z <= child->maxBounds.z) {
//                 return child->search(point);
//             }
//         }
//
//         return nullptr; // Should not reach here under normal circumstances
//     }
// };
//
// // Function to randomly subdivide the node and its descendants
// void randomSubdivide(OctreeNode* node) {
//     if (node->isLeaf) {
//         // Randomly decide whether to subdivide or not
//         if (rand() % 2 == 1) {
//             node->subdivide();
//             for (OctreeNode* child : node->children) {
//                 if (rand() % 2 == 1) {
//                    child->subdivide();
//                 }
//             }
//         }
//     }
// }
//
// // Function to print the octree structure
// void printOctreeStructure(OctreeNode* node, int depth = 0) {
//     for (int i = 0; i < depth; ++i) {
//         std::cout << "  ";
//     }
//
//     std::cout << "Node Bounds: ";
//     std::cout << "(" << node->minBounds.x << ", " << node->minBounds.y << ", " << node->minBounds.z << ") to ";
//     std::cout << "(" << node->maxBounds.x << ", " << node->maxBounds.y << ", " << node->maxBounds.z << ")";
//
//     if (node->isLeaf) {
//         std::cout << " [Leaf]\n";
//     } else {
//         std::cout << " [Non-Leaf]\n";
//         for (OctreeNode* child : node->children) {
//             printOctreeStructure(child, depth + 1);
//         }
//     }
// }
//
// int main() {
//     // Create the root node representing the unit cube
//     OctreeNode root({0.0f, 0.0f, 0.0f}, {1.0f, 1.0f, 1.0f});
//
//     // Randomly subdivide the root node and its descendants
//     srand(time(0));  // Seed for random number generation
//     for(int ii=0 ; ii<3 ; i++){
//       randomSubdivide(&root);
//     }
//
//     // Print the octree structure
//     std::cout << "Octree Structure:\n";
//     printOctreeStructure(&root);
//
//     // Example: search for the smallest sub-cube containing a given point
//     Point3D searchPoint(0.7f, 0.3f, 0.5f);
//     OctreeNode* result = root.search(searchPoint);
//
//     // Display the result
//     std::cout << "\nSearch Point: (" << searchPoint.x << ", " << searchPoint.y << ", " << searchPoint.z << ")\n";
//     std::cout << "Smallest Sub-cube Bounds: ";
//     std::cout << "(" << result->minBounds.x << ", " << result->minBounds.y << ", " << result->minBounds.z << ") to ";
//     std::cout << "(" << result->maxBounds.x << ", " << result->maxBounds.y << ", " << result->maxBounds.z << ")\n";
//
//     return 0;
// }
