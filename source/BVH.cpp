//
// Created by Edge on 2020/10/17.
//
#include "BVH.h"
#include <queue>
#include "randomstream.h"

static const int PC_DEBUGLEVEL = 1;

using namespace std;

namespace Manta {
    const Vec3 Extent::m_planeSetNormals[7] = {
            {1,                  0,                  0},
            {0,                  1,                  0},
            {0,                  0,                  1},
            {sqrtf(3.0) / 3.0f,  sqrtf(3.0) / 3.0f,  sqrtf(3.0) / 3.0f},
            {-sqrtf(3.0) / 3.0f, sqrtf(3.0) / 3.0f,  sqrtf(3.0) / 3.0f},
            {-sqrtf(3.0) / 3.0f, -sqrtf(3.0) / 3.0f, sqrtf(3.0) / 3.0f},
            {sqrtf(3.0) / 3.0f,  -sqrtf(3.0) / 3.0f, sqrtf(3.0) / 3.0f}
    };

    void Extent::extendBy(const Extent &extent) {
        for (int i = 0; i < 7; ++i) {
            m_d[i][0] = min(m_d[i][0], extent.m_d[i][0]);
            m_d[i][1] = max(m_d[i][1], extent.m_d[i][1]);
        }
    }

    bool
    Extent::isHit(const Ray &ray, Real precomputedNumerator[], Real precomputeDenominator[], Real&tNear, Real&tFar,
                  uint8_t &planeIndex) {

        for (int i = 0; i < 7; ++i) {
            Real tn = (m_d[i][0] - precomputedNumerator[i]) / precomputeDenominator[i];
            Real tf = (m_d[i][1] - precomputedNumerator[i]) / precomputeDenominator[i];
            if (precomputeDenominator[i] < 0) swap(tn, tf);
            if (tn > tNear) tNear = tn, planeIndex = i;
            if (tf < tFar) tFar = tf;
            if (tNear > tFar) return false;
        }
        return true;

    }

    Octree::Octree(const Extent &extent) {
        // constructor compute bound by first computing centroid and then add dimension to form bound
        Real xDiff = extent.m_d[0][1] - extent.m_d[0][0];
        Real yDiff = extent.m_d[1][1] - extent.m_d[1][0];
        Real zDiff = extent.m_d[2][1] - extent.m_d[2][0];
        Real dim = max(xDiff, max(yDiff, zDiff)) / 2.0;
        // technically, this is not true centroid, just use this for computational efficiency
        Vec3 centroid = Vec3(extent.m_d[0][0] + extent.m_d[0][1], extent.m_d[1][0] + extent.m_d[1][1],
                          extent.m_d[2][0] + extent.m_d[2][1]) / 2.0;
        m_bound[0] = (centroid - Vec3(dim, dim, dim));
        m_bound[1] = (centroid + Vec3(dim, dim, dim));
        m_root = new OctreeNode;
    }

    void Octree::insert(OctreeNode *node, const Extent *extent, const Vec3&boundMin,
                        const Vec3&boundMax, int depth) {
        if (node->m_isLeaf) {
            // insert extent into node or reallocate node extents' to child node
            if (node->m_data.empty() || depth == 16) {
                node->m_data.push_back(extent);
            } else {
                node->m_isLeaf = false;
                while (!node->m_data.empty()) {
                    // recursion and go to non leaf handling part
                    insert(node, node->m_data.back(), boundMin, boundMax, depth);
                    node->m_data.pop_back();
                }
                insert(node, extent, boundMin, boundMax, depth);
            }
        } else {
            // compute centroid and determine which child this extent should go
            Vec3 extentCentroid ((extent->m_d[0][0] + extent->m_d[0][1]) * 0.5f, (extent->m_d[1][0] + extent->m_d[1][1]) * 0.5f,
                                            (extent->m_d[2][0] + extent->m_d[2][1]) * 0.5f );
            Vec3 nodeCentroid = (boundMin + boundMax) * 0.5f;
            int childIndex = 0;
            if (extentCentroid.x > nodeCentroid.x) childIndex += 4;
            if (extentCentroid.y > nodeCentroid.y) childIndex += 2;
            if (extentCentroid.z > nodeCentroid.z) childIndex += 1;
            Vec3 childBoundMin, childBoundMax;
            // compute child bound and pass down
            computeChildBound(childIndex, nodeCentroid, boundMin, boundMax, childBoundMin, childBoundMax);
            if (node->m_child[childIndex] == nullptr) node->m_child[childIndex] = new OctreeNode;
            insert(node->m_child[childIndex], extent, childBoundMin, childBoundMax, depth + 1);
        }
    }

    void Octree::computeChildBound(int index, const Vec3&nodeCentroid, const Vec3&boundMin, const Vec3&boundMax,
        Vec3&pMin, Vec3&pMax) {
        pMin.x = (index & 4) ? nodeCentroid.x : boundMin.x;
        pMin.y = (index & 2) ? nodeCentroid.y : boundMin.y;
        pMin.z = (index & 1) ? nodeCentroid.z : boundMin.z;
        pMax.x = (index & 4) ? boundMax.x : nodeCentroid.x;
        pMax.y = (index & 2) ? boundMax.y : nodeCentroid.y;
        pMax.z = (index & 1) ? boundMax.z : nodeCentroid.z;
    }

    void Octree::build(OctreeNode *node, const Vec3&boundMin, const Vec3&boundMax) {
        if (node->m_isLeaf) {
            // update node extent
            for (auto extent:node->m_data) {
                node->m_nodeExtent.extendBy(*extent);
            }
        } else {
            for (int i = 0; i < 8; ++i) {
                if (node->m_child[i]) {
                    Vec3 nodeCentroid = (boundMin + boundMax) * 0.5f;
                    Vec3 childBoundMin, childBoundMax;
                    // compute child bound and pass down
                    computeChildBound(i, nodeCentroid, boundMin, boundMax, childBoundMin, childBoundMax);
                    build(node->m_child[i], childBoundMin, childBoundMax);
                    // collect child extent compute result
                    node->m_nodeExtent.extendBy(node->m_child[i]->m_nodeExtent);
                }
            }
        }
    }

    MeshBVH::MeshBVH(FluidSolver* parent) : Mesh(parent) {

    }

    MeshBVH* MeshBVH::clone() {
        // TODO
        return nullptr;
    }

    void MeshBVH::load(std::string name, bool append) {
        Mesh::load(name, append);

        Extent sceneExtent;
        m_extent.resize(numTris(), sceneExtent);

        for (int j = 0; j < numTris(); ++j) {
            for (int k = 0; k < 7; ++k) {
                computeTriangleBound(Extent::m_planeSetNormals[k], mTris[j], m_extent[j].m_d[k][0],
                                     m_extent[j].m_d[k][1]);
            }
            m_extent[j].m_triangleIndex = j;
            sceneExtent.extendBy(m_extent[j]);
        }
        m_octree = new Octree(sceneExtent);
        for (int j = 0; j < numTris(); ++j) {
            m_octree->insert(&m_extent[j]);
        }
        m_octree->build();
    }

    void
    MeshBVH::computeTriangleBound(const Vec3&planeNormal, const Triangle &triangle,
        Real &dNear,
        Real &dFar) {
        for (int i = 0; i < 3; ++i) {
            const Vec3 &triangleCoord = mNodes[triangle.c[i]].pos;
            auto d = dot(planeNormal, triangleCoord);
            if (d < dNear) dNear = d;
            if (d > dFar) dFar = d;
        }
    }

    bool MeshBVH::isHitExact(Real tmin, int index, const Ray &ray) const {
        const Real EPSILON = 1e-7;
        const Triangle& triangle = mTris[index];
        const Vec3(&point)[3] = { mNodes[triangle.c[0]].pos, mNodes[triangle.c[1]].pos,
                                   mNodes[triangle.c[2]].pos };
        Vec3 planeVector[2] = {point[1] - point[0],
                                   point[2] - point[0]};
        Vec3 h = cross(ray.velocity, planeVector[1]);
        Real a = dot(planeVector[0], h);
        if (std::abs(a) < EPSILON)
            return false;
        Real f = 1.0 / a;
        Vec3 s = ray.origin - point[0];
        Real u = f * dot(s, h);
        if (u < 0.0 || u > 1.0)
            return false;
        Vec3 q = cross(s, planeVector[0]);
        Real v = f * dot(ray.velocity, q);
        if (v < 0.0 || v + u > 1.0)
            return false;

        Real t = f * dot(planeVector[1], q);
        if (t > EPSILON && t > tmin) {
            // ray.origin + t * ray.velocity is intersection point
            return true;
        }
        return false;
    }

    bool MeshBVH::isInsideMesh(const Vec3 &pos, Real tmin) const {
        // We assume that mesh is closed, then basically we can use arbitrary ray to test what times it hit the mesh,
        // if the amounts is odd, then the ray's origin is inside the mesh, and vice versa.
        Ray testRay(pos);
        RandomStream mRand(9832);
        testRay.velocity = Vec3(1,1,1) - 2.0 * mRand.getVec3();// arbitrary ray
        Real precomputedNumerator[7], precomputedDenominator[7];
        for (int i = 0; i < 7; ++i) {
            precomputedNumerator[i] = dot(Extent::m_planeSetNormals[i], testRay.origin);
            precomputedDenominator[i] = dot(Extent::m_planeSetNormals[i], testRay.velocity);
        }
        Real tNear = -1e9, tFar = 1e9;
        uint8_t planeIndex;
        if (!m_octree->m_root->m_nodeExtent.isHit(testRay, precomputedNumerator, precomputedDenominator, tNear, tFar,
                                                  planeIndex)) {
            return false;
        }
//        debMsg("isInsideMesh test pass", PC_DEBUGLEVEL);
        priority_queue<Octree::QueueElement> pq;
        pq.push({Octree::QueueElement(m_octree->m_root, 0)});
        int hitAmount = 0;
        // FIXME refactor
        while (!pq.empty()) {
            const OctreeNode *node = pq.top().node;
            pq.pop();
            if (node->m_isLeaf) {
                for (int i = 0; i < node->m_data.size(); ++i) {
                    if (isHitExact(tmin, node->m_data[i]->m_triangleIndex, testRay)) {
                        hitAmount++;
                    }
                }
            } else {
                for (int i = 0; i < 8; ++i) {
                    if (node->m_child[i]) {
                        Real tNearChild = 0, tFarChild = tFar;
                        if (node->m_child[i]->m_nodeExtent.isHit(testRay, precomputedNumerator, precomputedDenominator,
                                                                 tNearChild, tFarChild, planeIndex)) {
                            Real t = (tNearChild < 0 && tFarChild >= 0) ? tFarChild : tNearChild;
                            pq.push(Octree::QueueElement(node->m_child[i], t));
                        }
                    }
                }
            }
        }
//        debMsg("isInsideMesh result=" << pos << " " << hitAmount, PC_DEBUGLEVEL);
        return hitAmount&1;
    }

    Vec3 MeshBVH::getMinBound() const {
        return m_octree->m_bound[0];
    }

    Vec3 MeshBVH::getMaxBound() const {
        return m_octree->m_bound[1];
    }
}