class Ball(object):
    def __init__(self, x=0.0, y=0.0, z=0.0, r=1.0):

        # note we make radius of sphere slightly shorter than it should be,
        # In physics simulations, it's common to slightly reduce the size of 
        # objects to prevent them from sticking together due to numerical precision errors 
        self.x = x
        self.y = y
        self.z = z
        self.radius = r
        self.visible = sphere(
            pos=vector(x, y, z), radius=r * 0.95, color=vector(0, 1, 0)
        )
