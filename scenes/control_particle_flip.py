#
# Simple fluid control based on flip with level set and basic resampling
#
from manta import *

# solver params
dim = 3
res = 64
# res = 128
gs = Vec3(res,res,res)
if (dim==2):
        gs.z=1
s = Solver(name='main', gridSize = gs, dim=dim)
s.timestep  = 0.8
minParticles = pow(2,dim)

# save particles for separate surface generation pass?
saveParts = False

# size of particles
radiusFactor = 1.0

# control params
# TODO fill param
mesh_sample_dist=1.0
attr_force_strength=3.0 # 2.0
vel_force_strength=0.2
volume=1.0

# prepare grids and particles
phi = s.create(LevelsetGrid)
flags = s.create(FlagGrid)

vel = s.create(MACGrid)
velOld = s.create(MACGrid)
pressure = s.create(RealGrid)
tmpVec3  = s.create(VecGrid)
tstGrid  = s.create(RealGrid)
tmpForce  = s.create(VecGrid)

pp = s.create(BasicParticleSystem)
pVel = pp.create(PdataVec3)
# test real value, not necessary for simulation
pTest    = pp.create(PdataReal)
mesh     = s.create(Mesh)

# acceleration data for particle nbs
pindex = s.create(ParticleIndexSystem)
gpi    = s.create(IntGrid)

ppControl = s.create(BasicParticleSystem)
ppCScale = ppControl.create(PdataReal)
ppCVel = ppControl.create(PdataVec3)

# show up weird
#mesh_obj = s.create(Mesh)

# scene setup
bWidth=1
flags.initDomain(boundaryWidth=bWidth)
basin = Box( parent=s, p0=gs*Vec3(0,0,0), p1=gs*Vec3(1,0.2,1))
pball = Sphere( parent=s , center=gs*Vec3(0.5,0.23,0.5), radius=0.02)
iball = Sphere( parent=s , center=gs*Vec3(0.4,0.7,0.5), radius=res*0.125)
ball = Sphere( parent=s , center=gs*Vec3(0.5,0.3,0.5), radius=res*0.125)
phi.setConst(1e10)
phi.join(basin.computeLevelset())
# phi.join(iball.computeLevelset())
# phi.join(box.computeLevelset())
flags.updateFromLevelset(phi)
sampleLevelsetWithParticles( phi=phi, flags=flags, parts=pp, discretization=2, randomness=0.05 )
# Setup control particle
ppControlTest = s.create(BasicParticleSystem)
testMesh     = s.create(MeshBVH)
testMesh.load(name="wolf.obj")
# sampleMeshWithParticles(mesh=testMesh, flags=flags, parts=ppControlTest, reset=True)
sampleMeshBVHWithParticles(mesh=testMesh, flags=flags, parts=ppControlTest, reset=True)
sampleMeshBVHWithParticles(mesh=testMesh, flags=flags, parts=ppControl, reset=True)
# sampleShapeSurfaceWithOneParticle(shape=pball, flags=flags, parts=ppControl, discretization=1, randomness=0.333, reset=True)
# sampleShapeSurfaceWithParticles(shape=ball, flags=flags, parts=ppControl, discretization=1, randomness=0.333, reset=True)
# sampleShapeWithParticles(shape=ball, flags=flags, parts=ppControl, discretization=1, randomness=0.333, reset=True)

# testing the real channel while resampling - original particles
# will have a value of 0.1, new particle will get a value from the tstGrid
testInitGridWithPos(tstGrid)
pTest.setConst( 0.1 )

if (GUI):
        gui = Gui()
        gui.show()
        gui.pause()

#main loop
for t in range(1000):
        mantaMsg('\nFrame %i, simulation time %f' % (s.frame, s.timeTotal))

        # FLIP
        pp.advectInGrid(flags=flags, vel=vel, integrationMode=IntRK4, deleteInObstacle=False )

        # make sure we have velocities throught liquid region
        mapPartsToMAC(vel=vel, flags=flags, velOld=velOld, parts=pp, partVel=pVel, weight=tmpVec3 )
        extrapolateMACFromWeight( vel=vel , distance=2, weight=tmpVec3 )  # note, tmpVec3 could be free'd now...
        markFluidCells( parts=pp, flags=flags )

        # create approximate surface level set, resample particles
        gridParticleIndex( parts=pp , flags=flags, indexSys=pindex, index=gpi )
        unionParticleLevelset( pp, pindex, flags, gpi, phi , radiusFactor )
        resetOutflow(flags=flags,parts=pp,index=gpi,indexSys=pindex)
        # extend levelset somewhat, needed by particle resampling in adjustNumber
        extrapolateLsSimple(phi=phi, distance=4, inside=True);

        # forces & pressure solve
        # no work
        if t%200<100:
                print("controled attraction force")
                gridParticleIndex( parts=ppControl , flags=flags, indexSys=pindex, index=gpi )
                getParticleAttractionForceScaleFactor(parts=ppControl, flags=flags, radius=2.5*mesh_sample_dist, volume=volume, target=ppCScale)
                addAttractionForce(parts=ppControl, index=gpi, indexSys=pindex, scaleFactor=ppCScale, flags=flags, radius=2.5*mesh_sample_dist, attractionForceStrength=attr_force_strength, vel=vel)
                setWallBcs(flags=flags, vel=vel)
                solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
                setWallBcs(flags=flags, vel=vel)

                if t%100>=50:
                        # TODO ppCVel need to be update
                        tmpForce.clear()
                        circleParticleVelocity(parts=ppControl, velParts=ppCVel, t=t)
                        genVelocityForce(parts=ppControl, velParts=ppCVel, index=gpi, indexSys=pindex, flags=flags, radius=2.5*mesh_sample_dist, velocityForceStrength=vel_force_strength, vel=vel, force=tmpForce)
                        applyForceOnVel(flags=flags, vel=vel, force=tmpForce)
                        setWallBcs(flags=flags, vel=vel)
                        solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
                        setWallBcs(flags=flags, vel=vel)

        addGravity(flags=flags, vel=vel, gravity=(0,-0.001,0))
        setWallBcs(flags=flags, vel=vel)
        solvePressure(flags=flags, vel=vel, pressure=pressure, phi=phi)
        setWallBcs(flags=flags, vel=vel)

        # set source grids for resampling, used in adjustNumber!
        pVel.setSource( vel, isMAC=True )
        pTest.setSource( tstGrid );
        adjustNumber( parts=pp, vel=vel, flags=flags, minParticles=1*minParticles, maxParticles=2*minParticles, phi=phi, radiusFactor=radiusFactor )

        # make sure we have proper velocities
        extrapolateMACSimple( flags=flags, vel=vel )

        flipVelocityUpdate(vel=vel, velOld=velOld, flags=flags, parts=pp, partVel=pVel, flipRatio=0.97 )

        if (dim==3):
                phi.createMesh(mesh)

        #s.printMemInfo()
        s.step()

        # generate data for flip03_gen.py surface generation scene
        if saveParts:
                pp.save( 'flipParts_%04d.uni' % t );

        # if 1 and (GUI):
        #         gui.screenshot( 'control_particle_%04d.png' % t );
        # gui.screenshot( 'result_vel/control_particle_%04d.png' % t );