from collections import namedtuple
import csv
import os
#from numpy import sqrt, cos, sin, arctan, pi, array
#from numpy.linalg import norm
from math import sqrt, cos, sin, pi

array = lambda i: i
#import scipy
#import scipy.stats

Galaxy = namedtuple('Galaxy', 'id pos e1 e2')
Halo = namedtuple('Halo', 'x y')

XMAX = 4200
YMAX = 4200

halo_mass_dist_params = (1.794052966939786, -0.92368781464069372, 28.690645093008492)
mass_prior = dict()
with open('mass_prior.csv') as f:
    for line in f:
        w, l = line.split(',')
        mass_prior[int(w)] = float(l)




def galaxies_near(galaxies, center, size):
    n = 0
    for g in galaxies:
        r = sqrt((center.x - g.x)**2 + (center.y - g.y)**2)
        if r < size:
            n += 1

    return n

def load_sky(n, test=False):

    if not test:
        file_template = 'data/Train_Skies/Training_Sky%s.csv'
    else:
        file_template = 'data/Test_Skies/Test_Sky%s.csv'
    sky = []
    with open(file_template % n) as f:
        f.readline() # Throw away header
        for l in csv.reader(f):
            params = map(float, l[1:])
            sky.append(Galaxy(l[0], array(params[:2]), params[2], params[3]))
            #sky = [Galaxy(l[0], array(float(l[1]), float(l[2])), float(l[3]), float(l[4])) for l in csv.reader(f)]

    return sky

def load_halos(filename='data/Training_halos.csv'):
    sky = []

    with open(filename) as f:
        f.readline() # Throw away headers
        reader = csv.reader(f)
        for l in reader:
            halos = [array([float(l[4+2*i]), float(l[4+2*i+1])]) for i in range(int(l[1]))]
            sky.append(halos)

    return sky


def distance(p1, p2):
    return norm(p1 - p2)

def ellip_transform(e1, e2, phi):
    cos_phi = cos(2*phi)
    sin_phi = sin(2*phi)
    e_tan =  -(e1 * cos_phi + e2 * sin_phi)
    e_cross = -(e1 * -sin_phi + e2 * cos_phi)
    #e_cross =  -(e1 * cos(2*(phi + pi/ 4)) + e2 * sin(2*(phi + pi/4)))
    return e_tan, e_cross

def warp(galaxy, halo, halo_weight=90):
    delta = halo - galaxy.pos
    phi = arctan(delta[1] / delta[0])
    e_tan, e_cross = ellip_transform(galaxy.e1, galaxy.e2, phi)
    e_tan += halo_weight / distance(galaxy.pos, halo)
    e1, e2 = ellip_transform(e_tan, e_cross, -phi)
    galaxy_warp = Galaxy(galaxy.id, galaxy.pos, e1, e2)
    
    return galaxy_warp


def unwarp(galaxy, halo, halo_weight=90):
    return warp(galaxy, halo, -halo_weight)


def unwarp2(galaxy, halo, halo_weight=90):
    delta = [galaxy.pos[0] - halo[0], galaxy.pos[1] - halo[1]]
    r32 = (delta[0]**2 + delta[1]**2)**1.5
    e1 = galaxy.e1 - halo_weight * (delta[1]**2 - delta[0]**2) / r32
    e2 = galaxy.e2 + halo_weight * 2 * delta[0] * delta[1] / r32
    return e1, e2

def unwarp_sky(sky, halo, halo_weight=90):
    sky_u = []
    for g in sky:
        e1, e2 = unwarp2(g, halo, halo_weight)
        sky_u.append(Galaxy(g.id, g.pos, e1, e2))
    return sky_u

def tangential_ellipticity(galaxy, center):
    delta = center - galaxy.pos
    phi = arctan(delta[1] / delta[0])
    e = -(galaxy.e1 * cos(2*phi) + galaxy.e2 * sin(2*phi))
    return e
    

def out_of_bounds(halo_vec, max_pos=4200):
    for p in halo_vec:
        if not 0 < p < max_pos:
            return True
    return False

def sky_likelihood(sky, halo_vec, weights=None):
    ll = 0
    if weights is None:
        weights = [90,] * (len(halo_vec) / 2)

    for w in weights:
        ll += mass_prior[int(w)]

    for g in sky:
        for i in range(len(halo_vec) / 2):
            e1, e2 = unwarp2(g, halo_vec[2*i:2*(i+1)], weights[i])
            g = Galaxy(g.id, g.pos, e1, e2)
        ll -= e1**2 + e2**2
    return ll



