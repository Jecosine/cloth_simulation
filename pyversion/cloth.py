def create_cloth(N, ballsize, nodes):
    """ create a cloth of N x N nodes

    Args:
        N (_type_): _description_
        ballsize (_type_): _description_
        nodes (_type_): _description_
    """
    #   create nodes in cloth
    for nx in range(N):
        x = nx * separation - (N - 1) * separation * 0.5 + offset
        for ny in range(N):
            y = ny * separation - (N - 1) * separation * 0.5 + offset
            # if you set radius to something you will see the nodes
            node = sphere(
                pos=vector(x, ballsize + 1.0, y), radius=0, color=vector(1, 1, 0)
            )
            node.force = vector(0.0, 0.0, 0.0)
            node.velocity = vector(0.0, 0.0, 0.0)
            node.oldforce = vector(0.0, 0.0, 0.0)
            nodes.append(node)
    # colour squares between nodes
    for nx in range(N - 1):
        for ny in range(N - 1):
            c = vector(random(), random(), random())
            pt1 = vertex(pos=nodes[nx * N + ny].pos, color=c)
            pt2 = vertex(pos=nodes[(nx + 1) * N + ny].pos, color=c)
            pt3 = vertex(pos=nodes[(nx + 1) * N + ny + 1].pos, color=c)
            pt4 = vertex(pos=nodes[nx * N + ny + 1].pos, color=c)
            nodes[nx * N + ny].fill = quad(v0=pt1, v1=pt2, v2=pt3, v3=pt4)


def compute_force(delta, gravity, separation, fcon):
    """The potential is compute by the Hooke's law:
    $$
    PE_{ij} = K * (R_{ij}-E_{ij})
    $$
    Then the force can be calculated:
    $$
    F_{x_{ij}} = K * \frac{(R_{ij}-E_{ij}) * (x_i - x_j)}{R_{ij}}
    $$

    where r is the distance between two nodes, r0 is the rest length of the spring

    Args:
        delta (float): calculate force using nodes within `delta` distance
        gravity (float): the gravity constant
        separation (float): the separation between nodes
        fcon (float): force constant determining how stiff the spring (or cloth) is

    Returns:
        float: Potential energy 
    """
    r12 = vector(0.0, 0.0, 0.0)
    PE = 0.0
    #   loop over nodes in x and y direction
    for nx in range(N):
        for ny in range(N):
            #   add gravitational force
            nodes[nx * N + ny].force = vector(0.0, -gravity, 0.0)
            #   for node (nx,ny) loop over surrounding nodes and eval force/PE
            for dx in range(max(nx - delta, 0), min(nx + delta + 1, N)):
                for dy in range(max(ny - delta, 0), min(ny + delta + 1, N)):
                    len = sqrt(float((nx - dx) ** 2 + (ny - dy) ** 2)) * separation
                    #   don't self interact
                    if nx != dx or ny != dy:
                        r12 = nodes[dx * N + dy].pos - nodes[nx * N + ny].pos
                        PE += fcon * (r12.mag - len) * (r12.mag - len)
                        nodes[nx * N + ny].force += fcon * r12.norm() * (r12.mag - len)
    return PE