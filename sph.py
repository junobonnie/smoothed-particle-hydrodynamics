from vectortools import *
from atom import *
import sys

L = 20 # smoothing length
K = 0.02*L**6 # isentropic exponent
alpha = 0.5
rho_0 = 0.01
mu = 0.0002*L**6

class Atom(Atom):
    def __init__(self, element, pos, vel = Vector(0, 0), rho = 0, p = 0, color = (0,0,0)):
        self.element = element
        self.pos = pos
        self.vel = vel
        self.rho = rho
        self.p = p   
        self.color = color

class Render(Render):
    def __init__(self, screen, width, height):
        pg.init()
        self.screen = screen
        self.width = width
        self.height = height
        self.render_vector = Vector(0, height)
        self.render_metric = Tensor(1, 0, 0, -1)
        self.origin_vector = Vector(width/2, height/2)

    def atom(self, atom):
        self.circle(atom.pos, atom.element.radius, atom.color)

class World(World):
    def __init__(self, t, atoms, walls, gravity = Vector(0, 0)):
        self.t = t
        self.atoms = atoms
        self.walls = walls
        self.gravity = gravity
            
class Simulator(Simulator):
    def __init__(self, dt, world, render):
        self.dt = dt
        self.world = world
        self.render = render
        self.count_screen = 0
        self.count_snapshot = 0
    
    def set_density(self):
        for atom in self.world.atoms:
            result = 0
            for other in self.world.atoms:
                r = atom.pos - other.pos
                if not atom == other and r.dot(r) < L**2:# and r.dot(r) > self.element.radius+other.element.radius:
                    result += (15*atom.element.mass/(3.14*L**6))*((L - abs(r))**3) #-other.element.mass/r.dot(r)*(r/abs(r))
            atom.rho = result
            #print('============')
            #print(atom.rho)

    def set_pressure(self):
        for atom in self.world.atoms:
            atom.p = max(K*(atom.rho-rho_0),0)
            #print(atom.rho)

    def pressure_acc(self, atom): 
        result = Vector(0, 0) 
        for other in self.world.atoms:
            r = atom.pos - other.pos
            if not atom == other and r.dot(r) < L**2: # and r.dot(r) > atom.element.radius+other.element.radius:
                #print(other.rho)
                result += (45/(3.14*L**6))*(r/(abs(r)+0.01))*((atom.p + other.p)/(2*other.rho+0.002))*((L-abs(r))**2) #-other.element.mass/r.dot(r)*(r/abs(r))
                #print(result)
        return result
    #Pᵢ = (− (45 M) / (π L⁶)) ∑ⱼ − (xⱼ − xᵢ) / (dᵢⱼ) (pⱼ + pᵢ) / (2 ρⱼ) (L − dᵢⱼ)² 

    def viscosity_acc(self, atom):
        result = Vector(0, 0)
        for other in self.world.atoms:
            r = atom.pos - other.pos
            if not atom == other and r.dot(r) < L**2:# and r.dot(r) > self.element.radius+other.element.radius:
                result += -alpha*r+(45*mu/(3.14*L**6))*((other.vel - atom.vel)/(other.rho+0.001))*(L-abs(r)) #-other.element.mass/r.dot(r)*(r/abs(r))
                #result += (15*mu/(2*3.14*L**3))*((other.vel - atom.vel)/(other.rho+0.01))*((-1.5*(r.dot(r)/L**3)+(2*abs(r)/(L**2))+(L/(2*r.dot(r))))) #-other.element.mass/r.dot(r)*(r/abs(r))
        return result
    #Vᵢ = (45 μ M) / (π L⁶) ∑ⱼ (uⱼ − uᵢ) / (ρⱼ) (L − dᵢⱼ)

    def main(self):
        x_ = []
        v_ = []
        self.set_density()
        self.set_pressure()
        for atom in self.world.atoms:
            new_v = atom.vel + self.pressure_acc(atom)*self.dt + self.viscosity_acc(atom)*self.dt + self.world.gravity*self.dt
            v_.append(new_v)
            x_.append(atom.pos + new_v*self.dt)
        
        count = 0
        for atom in self.world.atoms:
            atom.color = (0,0,max(255-int(5*255*atom.rho),0))
            #print(atom.rho)
            atom.pos = x_[count]
            atom.vel = v_[count]
            count = count + 1
            
if __name__ == '__main__':
    width = 1000
    height = 800

    screen = pg.display.set_mode((width, height))
    render = Render(screen, width, height)
    clock = pg.time.Clock()

    black = pg.Color('black')
    white = pg.Color('white')
    red = pg.Color('red')
    green = pg.Color('green')
    blue = pg.Color('blue')

    wall1 = Wall(1000, 50, 0, Vector(-500, -400), red)
    wall2 = Wall(50, 800, 0, Vector(-500, -400), blue)
    wall3 = Wall(50, 800, 0, Vector(450,-400), blue)
    wall4 = Wall(1000, 50, 0, Vector(-500, 350), blue)
    wall5 = Wall(1000, 50, m.pi/4, Vector(0, -400), blue)

    e1 = Element(name = 'Helium', mass = 100, radius = 5, color = red)
    e2 = Element(name = 'Uranium', mass = 100, radius = 5, color = red)
   
    # atom1 = Atom(e1, Vector(-200, 0), Vector(50, 0))
    # atom2 = Atom(e1, Vector(0, 0))
    # atom3 = Atom(e1, Vector(25, -10))
    # atom4 = Atom(e1, Vector(25, 10))
    # atom5 = Atom(e1, Vector(50, -20))
    # atom6 = Atom(e1, Vector(50, 0))
    # atom7 = Atom(e1, Vector(50, 20))

    walls = [wall1, wall2, wall3, wall4]#, wall5]
    atoms = [] # [atom1, atom2, atom3, atom4, atom5, atom6, atom7]
    
    import random as r
    import math as m
    for i in range(300):
        atoms.append(Atom(e1, Vector(500*r.random()+400-650, 100*r.random()+300-100)))
    
  
        

    gravity = Vector(0, -10)

    world = World(0, atoms, walls, gravity)

    simulator = Simulator(0.05, world, render)
    
    while True:
        t = simulator.clock()
        simulator.draw_background(white)
        simulator.draw_grid(100)
        simulator.draw_wall()
        simulator.main()
        simulator.atom_wall_collision()
        # simulator.atom_atom_collision()
        # simulator.atom_atom_fusion()
        simulator.draw_atom()

        # render.text('pos = (%.2f, %.2f)'%(atoms[0].pos.x, atoms[0].pos.y) , None, 30, Vector(atoms[0].pos.x -100, atoms[0].pos.y - 30), black)
        # render.text('vel = (%.2f, %.2f)'%(atoms[0].vel.x, atoms[0].vel.y) , None, 30, Vector(atoms[0].pos.x -100, atoms[0].pos.y - 50), black)

        # render.text('pos = (%.2f, %.2f)'%(atoms[50].pos.x, atoms[50].pos.y) , None, 30, Vector(atoms[50].pos.x -100, atoms[50].pos.y - 30), blue)
        # render.text('vel = (%.2f, %.2f)'%(atoms[50].vel.x, atoms[50].vel.y) , None, 30, Vector(atoms[50].pos.x -100, atoms[50].pos.y - 50), blue)
    
        for event in pg.event.get():
            if event.type == pg.QUIT:
                sys.exit()
        clock.tick(100)        #pg.display.update()
        
        # simulator.save_screen('images/gravity_model_demo_1')
        
        #if t > 20:
        #    break
        #print(t)
        #pg.display.flip()
        #view = pg.surfarray.array3d(screen)
        #  convert from (width, height, channel) to (height, width, channel)
        #view = view.transpose([1, 0, 2])

        #  convert from rgb to bgr
        #img_bgr = cv2.cvtColor(view, cv2.COLOR_RGB2BGR)

        #Display image, clear cell every 0.5 seconds
        #cv2_imshow(img_bgr)
        #time.sleep(0.5)
        #output.clear()

        pg.display.update()
        
        simulator.save_screen('images/raindrop')
        #simulator.save_snapshot('snapshots/pocket_ball_demo', 99)
