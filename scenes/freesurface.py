#
# Simple example for free-surface simulation
# with MacCormack advection

from manta import *

# solver params
res = 64
gs = vec3(res,res,res)
s = Solver(name='main', gridSize = gs)
s.timestep = 0.25

# prepare grids and particles
flags = s.create(FlagGrid)
vel = s.create(MACGrid)
phi = s.create(LevelsetGrid)
pressure = s.create(RealGrid)
mesh = s.create(Mesh)

# scene setup
flags.initDomain()
drop = s.create(Sphere, center=gs*vec3(0.5,0.5,0.5), radius=res*0.15)
basin = s.create(Box, p0=gs*vec3(0,0,0), p1=gs*vec3(1,0.2,1))
fluidbox.computeLevelset(phi)
flags.updateFromLevelset(phi)
flip.adjustNumber(vel=vel, flags=flags, minParticles=8, maxParticles=30)
    
if (GUI):
    gui = Gui()
    gui.show()
    
#main loop
for t in range(200):
    
    # FLIP advect and writeback
    advectSemiLagrange(flags=flags, vel=vel, grid=vel, order=2)
    
    addGravity(flags=flags, vel=vel, gravity=(0,-0.002,0))
    
    # pressure solve
    setWallBcs(flags=flags, vel=vel)    
    setLiquidBcs(flags=flags, vel=vel)
    solvePressure(flags=flags, vel=vel, pressure=pressure)
    setLiquidBcs(flags=flags, vel=vel)    
    setWallBcs(flags=flags, vel=vel)
    
    # update and advect levelset
    phi.reinitMarching(flags=flags, ignoreWalls=False)
    advectSemiLagrange(flags=flags, vel=vel, grid=phi, order=2)
    flags.updateFromLevelset(phi)
    
    # note: these meshes are created by fast marching only, should smooth
    #       geometry and normals before rendering
    if (meshing):
        phi.createMesh(mesh)
        mesh.save('phi%04d.bobj.gz' % t)
    
    s.step()
