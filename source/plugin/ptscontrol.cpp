/******************************************************************************
 *
 * MantaFlow fluid solver framework
 * Copyright 2011 Tobias Pfaff, Nils Thuerey 
 *
 * This program is free software, distributed under the terms of the
 * Apache License, Version 2.0 
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Plugins for control particle
 *
 ******************************************************************************/
#include "vectorbase.h"
#include "grid.h"
#include "mesh.h"
#include "particle.h"
#include "kernel.h"
#include "shapes.h"
#include "BVH.h"

static const int PC_DEBUGLEVEL = 1;

using namespace std;
namespace Manta {

//! sample a shape with particles, use reset to clear the particle buffer,
//! and skipEmpty for a continuous inflow (in the latter case, only empty cells will
//! be re-filled once they empty when calling sampleShapeWithParticles during
//! the main loop).
PYTHON()
void sampleShapeSurfaceWithParticles(const Shape &shape, const FlagGrid &flags, BasicParticleSystem &parts,
                                     const int discretization, const Real randomness, const bool reset = false,
                                     const bool refillEmpty = false,
                                     const LevelsetGrid *exclude = NULL) {
    const bool is3D = flags.is3D();
    const Real jlen = randomness / discretization;
    const Vec3 disp(1.0 / discretization, 1.0 / discretization, 1.0 / discretization);
    RandomStream mRand(9832);

    if (reset) {
        parts.clear();
        parts.doCompress();
    }

    FOR_IJK_BND(flags, 0) {
        if (flags.isObstacle(i, j, k)) continue;
        if (refillEmpty && flags.isFluid(i, j, k)) continue;
        const Vec3 pos(i, j, k);
        for (int dk = 0; dk < (is3D ? discretization : 1); dk++)
        for (int dj = 0; dj < discretization; dj++)
        for (int di = 0; di < discretization; di++) {
            Vec3 subpos = pos + disp * Vec3(0.5 + di, 0.5 + dj, 0.5 + dk);
            subpos += jlen * (Vec3(1, 1, 1) - 2.0 * mRand.getVec3());
            if (!is3D) subpos[2] = 0.5;
            if (exclude && exclude->getInterpolated(subpos) <= 0.) continue;
            if (!shape.isNear(subpos)) continue;
            parts.addBuffered(subpos);
        }
    }

    parts.insertBufferedParticles();
}

//! sample a shape with particles, use reset to clear the particle buffer,
//! and skipEmpty for a continuous inflow (in the latter case, only empty cells will
//! be re-filled once they empty when calling sampleShapeWithParticles during
//! the main loop).
PYTHON()
void sampleShapeSurfaceWithOneParticle(const Shape &shape, const FlagGrid &flags, BasicParticleSystem &parts,
                                       const int discretization, const Real randomness, const bool reset = false,
                                       const bool refillEmpty = false,
                                       const LevelsetGrid *exclude = NULL) {
    const bool is3D = flags.is3D();
    const Real jlen = randomness / discretization;
    const Vec3 disp(1.0 / discretization, 1.0 / discretization, 1.0 / discretization);
    RandomStream mRand(9832);

    if (reset) {
        parts.clear();
        parts.doCompress();
    }

    FOR_IJK_BND(flags, 0) {
        if (flags.isObstacle(i, j, k)) continue;
        if (refillEmpty && flags.isFluid(i, j, k)) continue;
        const Vec3 pos(i, j, k);
        for (int dk = 0; dk < (is3D ? discretization : 1); dk++)
        for (int dj = 0; dj < discretization; dj++)
        for (int di = 0; di < discretization; di++) {
            Vec3 subpos = pos + disp * Vec3(0.5 + di, 0.5 + dj, 0.5 + dk);
            subpos += jlen * (Vec3(1, 1, 1) - 2.0 * mRand.getVec3());
            if (!is3D) subpos[2] = 0.5;
            if (exclude && exclude->getInterpolated(subpos) <= 0.) continue;
            if (!shape.isNear(subpos)) continue;
            parts.addBuffered(subpos);
            parts.insertBufferedParticles();
            return;
        }
    }
}

PYTHON()
void sampleMeshWithParticles(Mesh &mesh, const FlagGrid& flags, BasicParticleSystem &parts, const bool reset = false,
                             const int particleFlag = -1) {
    if (reset) {
        parts.clear();
        parts.doCompress();
    }

    for (IndexInt i = 0; i < mesh.numNodes(); ++i) {
        const Node node = mesh.nodes(i);
        const Vec3 &pos = node.pos;
        if (particleFlag < 0) {
            parts.addBuffered(pos);
        } else {
            parts.addBuffered(pos, particleFlag);
        }
    }
    // this may help
//    applyMeshToGrid()
    //    FOR_IJK_BND(flags, 0) {
//        if ( flags.isObstacle(i,j,k) ) continue;
//        if ( phi(i,j,k) < 1.733 ) {
//            const Vec3 pos (i,j,k);
//            for (int dk=0; dk<(is3D ? discretization : 1); dk++)
//                for (int dj=0; dj<discretization; dj++)
//                    for (int di=0; di<discretization; di++) {
//                        Vec3 subpos = pos + disp * Vec3(0.5+di, 0.5+dj, 0.5+dk);
//                        subpos += jlen * (Vec3(1,1,1) - 2.0 * mRand.getVec3());
//                        if(!is3D) subpos[2] = 0.5;
//                        if( phi.getInterpolated(subpos) > 0. ) continue;
//                        if(particleFlag < 0){
//                            parts.addBuffered(subpos);
//                        }
//                        else{
//                            parts.addBuffered(subpos, particleFlag);
//                        }
//                    }
//        }
//    }

    parts.insertBufferedParticles();
}

PYTHON()
void sampleMeshBVHWithParticles(const MeshBVH &mesh, const FlagGrid& flags, BasicParticleSystem &parts, const bool reset = false,
                             const int particleFlag = -1) {
    if (reset) {
        parts.clear();
        parts.doCompress();
    }

    RandomStream mRand(9832);
    int sampleAmount = 5;
    Vec3 minBound = mesh.getMinBound(), maxBound = mesh.getMaxBound();
    for(int k = minBound.z; k <= maxBound.z; ++k)
    for(int j = minBound.y; j <= maxBound.y; ++j)
    for(int i = minBound.x; i <= maxBound.x; ++i) {
        Vec3 pos(i, j, k);
        for(int di = 0; di < sampleAmount; ++di) {
            Vec3 subpos = pos + (Vec3(1,1,1) - 2.0 * mRand.getVec3());
            if(mesh.isInsideMesh(subpos)) {
                parts.addBuffered(subpos);
            }
        }
    }
    parts.insertBufferedParticles();
}

Real getLinearFallOff(Real d, Real h) {
    if (d <= h / 2.0) {
        return 1.0;
    } else if (d > h / 2.0 && d < h) {
        return 2 - 2 * d / h;
    } else {
        return 0;
    }
}

Real getNormalizedSpline(Real d, Real h) {
    if (d < h) {
        return 315.0 * pow(pow(h, 2.0) - pow(d, 2.0), 3.0) / (64.0 * acos(-1) * pow(h, 9.0));
    } else {
        return 0;
    }
}

KERNEL(pts)
void
KnGetParticleAttractionForceScaleFactor(const BasicParticleSystem &parts, const FlagGrid &flags, const Real radius,
                                        Real volume, ParticleDataImpl <Real> &target) {
    const Vec3 &pos = parts.getPos(idx);
    Real totalInfl = 0.0;
    int r = int(radius) + 1;
    int rZ = flags.is3D() ? r : 0;
    int times = 0;
    for (int zj = pos[2] - rZ; zj <= pos[2] + rZ; zj++)
    for (int yj = pos[1] - r; yj <= pos[1] + r; yj++)
    for (int xj = pos[0] - r; xj <= pos[0] + r; xj++) {
        times++;
        if (flags.isInBounds(Vec3i(xj, yj, zj)) && flags.isFluid(xj, yj, zj)) {
            totalInfl += volume * getNormalizedSpline(norm(Vec3(xj, yj, zj) - pos), radius);
//            debMsg("KnGetParticleAttractionForceScaleFactor linear fall off param=" << norm(Vec3(xj, yj, zj) - pos) << ", linear fall off=" << getNormalizedSpline(norm(Vec3(xj, yj, zj) - pos), radius) << ", totalInfl=" << totalInfl << ", times=" << times, PC_DEBUGLEVEL);
        }
    }
    target[idx] = 1.0 - min(Real(1.0), totalInfl);
//    debMsg("KnGetParticleAttractionForceScaleFactor idx=" << idx << ", totalInfluence=" << totalInfl, PC_DEBUGLEVEL);
//    if(target[idx] > 0.0) {
//        debMsg("KnGetParticleAttractionForceScaleFactor result=" << target[idx], PC_DEBUGLEVEL);
//    }

}

PYTHON()
void
getParticleAttractionForceScaleFactor(const BasicParticleSystem &parts, const FlagGrid &flags, const Real radius,
                                      Real volume, ParticleDataImpl <Real> &target) {
    KnGetParticleAttractionForceScaleFactor(parts, flags, radius, volume, target);
}

KERNEL()
void KnGetParticleAttractionForce(const BasicParticleSystem &parts, const Grid<int> &index,
                                  const ParticleIndexSystem &indexSys, const ParticleDataImpl <Real> &scaleFactor,
                                  const FlagGrid &flags, const Real radius, const Real attractionForceStrength,
                                  Grid <Vec3> &vel, bool scale) {

    if (!flags.isFluid(i, j, k)) return;

    Vec3 cellPos = Vec3(i, j, k);
    int r = int(1. * radius) + 1;
    int rZ = vel.is3D() ? r : 0;

    float gridScale = (scale) ? flags.getDx() : 1;

    for (int zj = k - rZ; zj <= k + rZ; zj++)
    for (int yj = j - r; yj <= j + r; yj++)
    for (int xj = i - r; xj <= i + r; xj++) {
        if (!vel.isInBounds(Vec3i(xj, yj, zj))) continue;

        IndexInt isysIdxS = index.index(xj, yj, zj);
        IndexInt pStart = index(isysIdxS), pEnd = 0;
        if (vel.isInBounds(isysIdxS + 1)) pEnd = index(isysIdxS + 1);
        else pEnd = indexSys.size();

        // now loop over particles in cell
        for (IndexInt p = pStart; p < pEnd; ++p) {
            const int psrc = indexSys[p].sourceIndex;
            const Vec3 pos = parts[psrc].pos;
            Vec3 direction = pos - cellPos;
            Real distance = norm(direction);
            Vec3 attractionForce = scaleFactor[psrc] * getNormalized(direction) *
                                   getNormalizedSpline(distance,
                                                       radius); //  * flags.getParent()->getDt() / gridScale

//            if(scaleFactor[psrc] * getNormalizedSpline(distance, radius) > 0.0) {
//                debMsg("KnGetParticleAttractionForce result=" << scaleFactor[psrc] * getNormalizedSpline(distance, radius) * flags.getParent()->getDt() / gridScale, PC_DEBUGLEVEL);
//            }
            vel(i, j, k) += attractionForceStrength * attractionForce;
        }
    }

}

PYTHON()
void addAttractionForce(
        const BasicParticleSystem &parts, const Grid<int> &index, const ParticleIndexSystem &indexSys,
        const ParticleDataImpl <Real> &scaleFactor, const FlagGrid &flags, const Real radius,
        const Real attractionForceStrength, Grid <Vec3> &vel, bool scale = true) {
    KnGetParticleAttractionForce(parts, index, indexSys, scaleFactor, flags, radius, attractionForceStrength, vel,
                                 scale);
}

KERNEL()
void KnGetParticleVelocityForce(const BasicParticleSystem &parts, const ParticleDataImpl <Vec3> &velParts,
                                const Grid<int> &index, const ParticleIndexSystem &indexSys,
                                const FlagGrid &flags, const Real radius, const Real velocityForceStrength,
                                const Grid <Vec3> &vel, Grid <Vec3> &force, bool scale) {
    // calculate near particle affect to current grid velocity force
    if (!flags.isFluid(i, j, k)) return;

    Vec3 cellVel = vel.get(i, j, k);
    int r = int(1. * radius) + 1;
    int rZ = vel.is3D() ? r : 0;

    float gridScale = (scale) ? flags.getDx() : 1;

    for (int zj = k - rZ; zj <= k + rZ; zj++)
    for (int yj = j - r; yj <= j + r; yj++)
    for (int xj = i - r; xj <= i + r; xj++) {
        if (!vel.isInBounds(Vec3i(xj, yj, zj))) continue;

        IndexInt isysIdxS = index.index(xj, yj, zj);
        IndexInt pStart = index(isysIdxS), pEnd = 0;
        if (vel.isInBounds(isysIdxS + 1)) pEnd = index(isysIdxS + 1);
        else pEnd = indexSys.size();

        // now loop over particles in cell
        for (IndexInt p = pStart; p < pEnd; ++p) {
            const int psrc = indexSys[p].sourceIndex;
            Vec3 direction = velParts[psrc] - cellVel;
            Vec3 velocityForce = direction * getNormalizedSpline(norm(direction), radius);
            force(i, j, k) += velocityForceStrength * velocityForce;
        }
    }
}

PYTHON()
void genVelocityForce(
        const BasicParticleSystem &parts, const ParticleDataImpl <Vec3> &velParts,
        const Grid<int> &index, const ParticleIndexSystem &indexSys,
        const FlagGrid &flags, const Real radius, const Real velocityForceStrength,
        const Grid <Vec3> &vel, Grid <Vec3> &force, bool scale = true) {
    KnGetParticleVelocityForce(parts, velParts, index, indexSys, flags, radius, velocityForceStrength, vel, force, scale);
}

KERNEL()
void KnApplyForceOnVel(const FlagGrid &flags, const Grid <Vec3> &force, Grid <Vec3> &vel) {
    if (!flags.isFluid(i, j, k)) return;
    vel(i, j, k) += force(i, j, k);
}

PYTHON()
void applyForceOnVel(const FlagGrid &flags, const Grid <Vec3> &force, Grid <Vec3> &vel) {
    KnApplyForceOnVel(flags, force, vel);
}

KERNEL(pts)
void KnCircleParticleVelocity(BasicParticleSystem &parts, ParticleDataImpl <Vec3> &velParts, float t) {
    float angle = t/acos(-1);
    velParts[idx] = {cos(angle), 0, sin(angle)};
    parts.setPos(idx, parts.getPos(idx) + velParts[idx]);
}

PYTHON()
void circleParticleVelocity(BasicParticleSystem &parts, ParticleDataImpl <Vec3> &velParts, float t) {
    KnCircleParticleVelocity(parts, velParts, t);
}

KERNEL(pts)
void KnComputeSmoothedFluidVelocityParticle(const BasicParticleSystem &parts, const FlagGrid &flags, const Grid <Vec3> &vel, const Real radius,
                                    ParticleDataImpl <Vec3> &target) {
    const Vec3 &pos = parts.getPos(idx);
    target[idx] = Vec3(0.0);
    int r = int(radius) + 1;
    int rZ = flags.is3D() ? r : 0;
    Real totalWeight = 0.0;
    for (int zj = pos[2] - rZ; zj <= pos[2] + rZ; zj++)
    for (int yj = pos[1] - r; yj <= pos[1] + r; yj++)
    for (int xj = pos[0] - r; xj <= pos[0] + r; xj++) {
        if (flags.isInBounds(Vec3i(xj, yj, zj)) && flags.isFluid(xj, yj, zj)) {
            Real weight = getNormalizedSpline(norm(Vec3(xj, yj, zj) - pos), radius);
            target[idx] += weight * vel.get(xj, yj, zj);
            totalWeight += weight;
        }
    }
    if (totalWeight >= 1e-6) {
        target[idx] /= totalWeight;
    }
}

PYTHON()
void computeSmoothedFluidVelocityParticle(const BasicParticleSystem& parts, const FlagGrid& flags, const Grid <Vec3>& vel, const Real radius,
    ParticleDataImpl <Vec3>& target) {
    KnComputeSmoothedFluidVelocityParticle(parts, flags, vel, radius, target);
}

KERNEL()
void KnComputeSmoothedFluidVelocity(const BasicParticleSystem &parts, const FlagGrid &flags, const Grid<int> &index, const ParticleIndexSystem &indexSys,
        const ParticleDataImpl <Vec3> &partVel, const Real radius, Grid <Vec3> &target) {
    if (!flags.isFluid(i, j, k)) return;
// FIXME change code
    Vec3 cellPos = Vec3(i, j, k);
    int r = int(1. * radius) + 1;
    int rZ = flags.is3D() ? r : 0;
    Real totalWeight = 0.0;

    for (int zj = k - rZ; zj <= k + rZ; zj++)
    for (int yj = j - r; yj <= j + r; yj++)
    for (int xj = i - r; xj <= i + r; xj++) {
        if (!flags.isInBounds(Vec3i(xj, yj, zj))) continue;

        IndexInt isysIdxS = index.index(xj, yj, zj);
        IndexInt pStart = index(isysIdxS), pEnd = 0;
        if (flags.isInBounds(isysIdxS + 1)) pEnd = index(isysIdxS + 1);
        else pEnd = indexSys.size();

        // now loop over particles in cell
        for (IndexInt p = pStart; p < pEnd; ++p) {
            const int psrc = indexSys[p].sourceIndex;
            const Vec3 pos = parts[psrc].pos;
            Real weight = getNormalizedSpline(norm(cellPos - pos), radius);
            totalWeight += weight;
            target(i, j, k) += partVel[psrc] * weight;
        }
    }
    if (totalWeight >= 1e-6) {
        target(i, j, k) /= totalWeight;
    }
}

PYTHON()
void computeSmoothedFluidVelocity(const BasicParticleSystem& parts, const FlagGrid& flags, const Grid<int>& index, const ParticleIndexSystem& indexSys,
    const ParticleDataImpl <Vec3>& partVel, const Real radius, Grid <Vec3>& target) {
    KnComputeSmoothedFluidVelocity(parts, flags, index, indexSys, partVel, radius, target);
}

} // end namespace
