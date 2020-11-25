//
// Created by Edge on 2020/10/17.
//

#ifndef _BVH_H
#define _BVH_H

#include "vectorbase.h"
#include "mesh.h"
#include <vector>

namespace Manta {
class Ray {
public:
    Ray() : origin(Vec3::Zero), velocity(Vec3::Zero) {}
    Ray(const Vec3& pos) : origin(pos), velocity(Vec3::Zero) {}
    Vec3 origin;
    Vec3 velocity;
};

class Extent {
public:
    Extent() {
        for (int i = 0; i < 7; ++i) {
            m_d[i][0] = 1e9;
            m_d[i][1] = -1e9;
        }
        m_triangleIndex = 0;
    }

    void extendBy(const Extent &extent);

    bool isHit(const Ray &ray, Real precomputedNumerator[], Real precomputeDenominator[], Real&tNear, Real&tFar,
               uint8_t &planeIndex);

    Real m_d[7][2];
    static const Vec3 m_planeSetNormals[7];
    int m_triangleIndex;
};

struct OctreeNode {
    OctreeNode *m_child[8] = {};
    std::vector<const Extent *> m_data;
    Extent m_nodeExtent;
    bool m_isLeaf = true;
    ~OctreeNode() {
        for(int i = 0;i < 8;++i) {
            if(m_child[i]) {
                delete m_child[i];
            }
        }
    }
};

class Octree {
public:
    Octree(const Extent &extent);

    void insert(const Extent *extent) {
        insert(m_root, extent, m_bound[0], m_bound[1], 0);
    }

    void build() {
        build(m_root, m_bound[0], m_bound[1]);
    }

    struct QueueElement {
        const OctreeNode *node;
        Real t;
        QueueElement(const OctreeNode *node, Real tHit) : node(node), t(tHit) {}
        friend bool operator < (const QueueElement &a, const QueueElement &b) {return a.t>b.t;};
    };
    ~Octree() {
        if(m_root) {
            delete m_root;
        }
    }

private:
    void insert(OctreeNode *node, const Extent *extent, const Vec3&boundMin, const Vec3&boundMax, int depth);

    void computeChildBound(int index, const Vec3&nodeCentroid, const Vec3&boundMin, const Vec3&boundMax, Vec3&pMin, Vec3&pMax);

    void build(OctreeNode *node, const Vec3&boundMin, const Vec3&boundMax);

public:
    Vec3 m_bound[2]; // for compute node centroid for later space partition
    OctreeNode *m_root = nullptr;
};

PYTHON() class MeshBVH : public Mesh {
public:
    PYTHON() MeshBVH(FluidSolver* parent);
    virtual ~MeshBVH() {
        delete m_octree;
    }
    virtual MeshBVH* clone();

    PYTHON() virtual void load (std::string name, bool append = false);

    void computeTriangleBound(const Vec3&planeNormal, const Triangle& triangle, Real &dNear, Real &dFar);

    bool isHitExact(Real tmin, int index, const Ray &ray) const;

    Vec3 getMinBound() const;

    Vec3 getMaxBound() const;

    PYTHON() bool isInsideMesh(const Vec3& pos, Real tmin=0.00001) const;

protected:
    std::vector<Extent> m_extent;

    Octree *m_octree;
};

}
#endif //_BVH_H
