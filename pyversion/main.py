from .utils import *
from .cloth import *
# default input parameters
N = int(20)
separation = float(1.0)
mass = float(1.0)
fcon = float(10.0)
interact = int(2)
gravity = float(0.981)
ballsize = float(5.0)
offset = float(0.0)
dt = float(0.02)
update = int(2)
verbose = int(1)

#   read command line arguments
N, separation, mass, fcon, interact, gravity, ballsize, offset, dt, update, verbose = read_arg(sys.argv[1:])

if __name__ == "__main__":
    scene.autoscale = 0
    myball = ball(0, 0, 0, ballsize)
    
    #  create nodes in cloth
    nodes = []
    create_cloth(N, ballsize, nodes)
    PE = compute_force(interact, gravity, separation, fcon)

    iter = 0
    maxit = 400
    while iter < maxit:
        iter += 1
        if verbose:
            print("iteration and potential energy ", iter, PE)

        #   Update coordinates using same MD velocity verlet
        for node in nodes:
            node.pos += dt * (node.velocity + dt * node.force * 0.5)
            node.oldforce = node.force

        #   apply constraints (move nodes to surface of ball)
        for node in nodes:
            dist = node.pos - vector(myball.x, myball.y, myball.z)
            if dist.mag < myball.radius:
                fvector = dist / dist.mag * myball.radius
                node.velocity -= node.velocity.dot(fvector) * fvector / (fvector.mag ** 2)
                node.pos = vector(myball.x, myball.y, myball.z) + fvector

        if iter % update == 0:
            #   update the view if necessary
            for nx in range(N - 1):
                for ny in range(N - 1):
                    nodes[nx * N + ny].fill.v0.pos = nodes[nx * N + ny].pos
                    nodes[nx * N + ny].fill.v1.pos = nodes[(nx + 1) * N + ny].pos
                    nodes[nx * N + ny].fill.v2.pos = nodes[(nx + 1) * N + ny + 1].pos
                    nodes[nx * N + ny].fill.v3.pos = nodes[nx * N + ny + 1].pos

            sleep(0.05)
        PE = compute_force(interact, gravity, separation, fcon)

        #   Update velocity using same MD velocity verlet
        damp = 0.995
        for node in nodes:
            node.velocity += dt * (node.force + node.oldforce) * 0.5 * damp