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

const int PC_DEBUGLEVEL = 1;

using namespace std;
namespace Manta {

PYTHON() void sampleMeshWithParticles(Mesh& mesh, BasicParticleSystem& parts, const bool reset=false, const int particleFlag=-1) {
    if(reset) {
        parts.clear();
        parts.doCompress();
    }

    for(IndexInt i=0; i<mesh.numNodes();++i) {
        Node& node = mesh.nodes(i);
        const Vec3& pos = node.pos;
        if(particleFlag < 0){
            parts.addBuffered(pos);
        }
        else{
            parts.addBuffered(pos, particleFlag);
        }
    }

    parts.insertBufferedParticles();
}

Real getLinearFallOff(Real d, Real h) {
    if(d <= h/2.0) {
        return 1.0;
    } else if(d > h/2.0 && d < h) {
        return 2-2*d/h;
    } else {
        return 0;
    }
}

KERNEL(pts)
void KnGetParticleAttractionForceScaleFactor(const BasicParticleSystem& parts, const FlagGrid& flags, const Real radius, Real volume, ParticleDataImpl<Real>& target) {
    const Vec3& pos = parts.getPos(idx);
    Real totalInfl = 0.0;
    int r  = int(radius) + 1;
    int rZ = flags.is3D() ? r : 0;
    for(int zj=pos[2]-rZ; zj<= pos[2]+rZ; zj++)
    for(int yj= pos[1]-r ; yj<= pos[1]+r ; yj++)
    for(int xj= pos[0]-r ; xj<= pos[0]+r ; xj++) {
        if(flags.isFluid(xj, yj, zj)) {
            totalInfl += volume * getLinearFallOff(norm(Vec3(xj, yj, zj) - pos), radius);
//            debMsg("KnGetParticleAttractionForceScaleFactor linear fall off param=" << sqrt(normSquare(Vec3(xj, yj, zj) - pos)) << ", linear fall off=" << getLinearFallOff(sqrt(normSquare(Vec3(xj, yj, zj) - pos)), radius) << ", radius=" << radius << ", volume=" << volume.get(xj, yj, zj), PC_DEBUGLEVEL);
        }
    }
    target[idx] = 1.0-min(Real(1.0), totalInfl);
//    if(target[idx] > 0.0) {
//        debMsg("KnGetParticleAttractionForceScaleFactor result=" << target[idx], PC_DEBUGLEVEL);
//    }

}

PYTHON()
void getParticleAttractionForceScaleFactor(const BasicParticleSystem& parts, const FlagGrid& flags, const Real radius, Real volume, ParticleDataImpl<Real>& target)
{
    KnGetParticleAttractionForceScaleFactor(parts, flags, radius, volume, target);
}

Real getNormalizedSpline(Real d, Real h) {
    if(d < h) {
        return 315.0 * pow(pow(h, 2.0) - pow(d, 2.0), 3.0) / (64.0*acos(-1)*pow(h, 9.0));
    } else {
        return 0;
    }
}

KERNEL()
void KnGetParticleAttractionForce(const BasicParticleSystem& parts, const Grid<int>& index, const ParticleIndexSystem& indexSys, const ParticleDataImpl<Real>& scaleFactor, const FlagGrid& flags, const Real radius, const Real attractionForceStrength, Grid<Vec3>& vel, bool scale) {

    if(!flags.isFluid(i, j, k)) return;

    Vec3 cellPos = Vec3(i, j, k);
    int   r = int(1. * radius) + 1;
    int   rZ = vel.is3D() ? r : 0;

    float gridScale = (scale) ? flags.getDx() : 1;

    for(int zj=k-rZ; zj<=k+rZ; zj++)
    for(int yj=j-r ; yj<=j+r ; yj++)
    for(int xj=i-r ; xj<=i+r ; xj++) {
        if (!vel.isInBounds(Vec3i(xj,yj,zj))) continue;

        IndexInt isysIdxS = index.index(xj,yj,zj);
        IndexInt pStart = index(isysIdxS), pEnd=0;
        if(vel.isInBounds(isysIdxS+1)) pEnd = index(isysIdxS+1);
        else                             pEnd = indexSys.size();

        // now loop over particles in cell
        for(IndexInt p=pStart; p<pEnd; ++p) {
            const int psrc = indexSys[p].sourceIndex;
            const Vec3 pos = parts[psrc].pos;
            Vec3 direction = pos - cellPos;
            Real distance = norm(direction);
            Vec3 attractionForce = scaleFactor[psrc] * getNormalized(direction) * getNormalizedSpline(distance, radius) * flags.getParent()->getDt() / gridScale;

            if(scaleFactor[psrc] * getNormalizedSpline(distance, radius) > 0.0) {
                debMsg("KnGetParticleAttractionForce result=" << scaleFactor[psrc] * getNormalizedSpline(distance, radius) * flags.getParent()->getDt() / gridScale, PC_DEBUGLEVEL);
            }
            vel(i, j, k) += attractionForceStrength * attractionForce;
        }
    }

}

PYTHON()
void addAttractionForce(
        const BasicParticleSystem& parts, const Grid<int>& index, const ParticleIndexSystem& indexSys, const ParticleDataImpl<Real>& scaleFactor, const FlagGrid& flags, const Real radius, const Real attractionForceStrength, Grid<Vec3>& vel, bool scale=true)
{
    KnGetParticleAttractionForce(parts, index, indexSys, scaleFactor, flags, radius, attractionForceStrength, vel, scale);
}

KERNEL(pts)
void KnGetParticleVelocityForce(const BasicParticleSystem& parts, const ParticleDataImpl<Vec3>& velParts, const FlagGrid& flags, const MACGrid& vel, const Real radius, const Real velocityForceStrength, Grid<Vec3>& force) {
    const Vec3& particleVel = velParts.get(idx);
    const Vec3& pos = parts.getPos(idx);
    int r  = int(radius) + 1;
    int rZ = flags.is3D() ? r : 0;
    for (int zj = pos[2] - rZ; zj <= pos[2] + rZ; zj++)
    for (int yj = pos[1] - r; yj <= pos[1] + r; yj++)
    for (int xj = pos[0] - r; xj <= pos[0] + r; xj++) {
        if(flags.isFluid(xj, yj, zj)) {
            // get grid pos
            Vec3 direction = particleVel - vel.get(xj, yj, zj);
            Vec3 velocityForce = direction * getNormalizedSpline(sqrt(normSquare(pos - Vec3(xj, yj, zj))), radius);
            force(xj, yj, zj) += velocityForceStrength * velocityForce;
        }
    }
}

PYTHON()
void genVelocityForce(
        const BasicParticleSystem& parts, const ParticleDataImpl<Vec3>& velParts, const FlagGrid& flags, const MACGrid& vel, const Real radius, const Real velocityForceStrength, Grid<Vec3>& force)
{
    KnGetParticleVelocityForce(parts, velParts, flags, vel, radius, velocityForceStrength, force);
}

} // end namespace
