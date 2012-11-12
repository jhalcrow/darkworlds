import csv
from halomodel import *
import sys
import multiprocessing
import random
import itertools
from math import exp, sqrt


def best_n(lls, cutoff=300, n=10):

    if len(lls) < n:
        return lls

    bests = [lls[0]]
    for ll, c, w in lls:
        too_close = False
        if len(c) > 1:
            raise
        for ll_b, c_b, w_b in bests:
            if list_distance(c[0], c_b[0]) < cutoff:
                too_close = True
                break
        if not too_close:
            bests.append((ll, c, w))

        if len(bests) >= n:
            return bests

    return bests


def best_bipolar(sky, center, size=100, dx=25, dw=10, wmin=5, wmax=250, dist_cutoff=1500):
    sky_filt = [g for g in sky if list_distance(g.pos, center) < dist_cutoff]

    lls = brute_force(sky_filt, 2, dx=dx, dw=dw,
        wmin=wmin, wmax=wmax,
        xmin=center[0] - size, xmax=center[0] + size, 
        ymin=center[1] - size, ymax=center[1] + size)
    ll, pair, weights = lls[0]
    return HaloCluster(pair, weights)


def brute_force(sky, num_halos=1, dx=50, dw=30, xmin=25, xmax=XMAX, ymin=25, ymax=YMAX, wmin=10, wmax=200):

    x = [float(i*dx) + xmin for i in range(0, int((xmax-xmin)/dx))]
    y = [float(i*dx)  + ymin for i in range(0, int((ymax-ymin)/dx))]
    weights = [float(i*dw) + wmin for i in range(0, int((wmax-wmin)/dw))]
    coords = [(xx, yy) for xx in x for yy in y]

    ll_best = None
    c_best = None
    w_best = None

    prior = mass_prior if num_halos == 1 else semi_mass_prior
    lls = []
    for clist in itertools.combinations(coords, num_halos):
        penalty = 0
        c = list(clist)
        for w in itertools.combinations_with_replacement(weights, num_halos):
            w = list(w)
            ll = sky_likelihood(sky, c, w, prior=prior)
            lls.append((ll, c, w))
            if ll_best is None or ll > ll_best:
                ll_best = ll
                c_best = clist
                w_best = w
    
    lls.sort(key = lambda x: -x[0])
    return lls


def optimize_weight(sky, halos, w_min=10, w_max=200, dw=10):

    weights = range(w_min, w_max, dw)
    ll_best = None
    weight = float(w_min)
    w_best = weight
    num_halos = len(halos)

    for weight in itertools.combinations(weights, num_halos):
        ll = sky_likelihood(sky, halos, weight)
        if ll_best < ll or not ll_best:
            ll_best = ll
            w_best = weight
     
    print 'Best weight: ', w_best
    return w_best

def center_of_gravity(cluster):
    total_weight = sum(cluster.weights)
    centroid_x = sum([h[0] * w / total_weight for (h, w) in zip(cluster.halos, cluster.weights)])
    centroid_y = sum([h[1] * w / total_weight for (h, w) in zip(cluster.halos, cluster.weights)])
    return [centroid_x, centroid_y]

def paired_brute_force(sky, num_halos, dx, dw, cutoff=300, 
        xmin=0, xmax=XMAX, ymin=0, ymax=YMAX, wmin=5, wmax=250, bipolar=False, prior=None):
    lls = brute_force(sky, dx=dx, dw=dw,  xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, wmin=wmin, wmax=wmax)

    best = best_n(lls, cutoff, 40)
    best_ll = None
    best_h = None
    best_w = None

    if bipolar:
        best = [best_bipolar(sky, center=c[0], size=100, dx=50, dw=10, wmin=10, wmax=100) 
                    for ll, c, w in best]

    for guesses in itertools.combinations(best, num_halos):
        if bipolar:
            h =  [hh for cluster in guesses for hh in cluster.halos]
            w =  [ww for cluster in guesses for ww in cluster.weights]
        else:
            h =  [h for g in guesses for h in g[1]]
            w = [w for g in guesses for w in g[2]]

        ll = sky_likelihood(sky, h, w)
        if not best_ll or ll > best_ll:
            best_ll = ll
            if bipolar:
                best_h = [center_of_gravity(c) for c in guesses]
                best_w = [sum(c.weights) for c in guesses]
                print best_h, best_w
            else:
                best_h = h
                best_w = w

    return (best_ll, best_h, best_w)



def load_test_halos():
    with open('data/Test_haloCounts.csv') as f:
        f.readline()
        counts = [int(l.split(',')[1]) for l in f]
    return counts

def predict(outfile, sky_ids, test=False, datadir='/data'):
    halos = load_halos()

    with open(outfile, 'w') as f:
        out = csv.writer(f)
        for sky_id in sky_ids:
            sky = load_sky(sky_id, test)
            if test:
                num_halos= load_test_halos()[sky_id - 1]
                print 'Sky', sky_id, 'Num halos', num_halos
            else:
                true_halos = load_halos()[sky_id - 1]
                num_halos = len(true_halos)
                print 'Sky', sky_id, 'True value:', true_halos
                best_w = optimize_weight(sky, true_halos)
                print 'True LL', sky_likelihood(sky, true_halos, best_w)
                
            if num_halos == 1:
                lls = brute_force(sky, num_halos, dx=50, dw=10)
                ll, prediction, _ = lls[0]
            else:
                ll, prediction, w = paired_brute_force(sky, num_halos, dx=75, dw=10, bipolar=False)

            print 'Found: ', prediction, ll
            if not test:
                err = best_halo_dist(prediction, true_halos) / num_halos
                print 'Error:', err

            halo_pos = [0,0,0,0,0,0]#np.zeros((6))
            for i in range(len(prediction)):
                halo_pos[i] = prediction[i]
            out_line = ['Sky%s' % sky_id] + [hc for h in prediction for hc in h]
            out.writerow(out_line)
            f.flush()

def best_halo_dist(halos1, halos2):
    best_dist = None
    for hperm in itertools.permutations(halos1):
        dist = sum([list_distance(h1, h2) for h1, h2 in zip(hperm, halos2)])
        if not best_dist or dist < best_dist:
            best_dist = dist
    return best_dist




if __name__ == '__main__':
    #predict(sys.argv[1], range(1, 121), True)
    #predict(sys.argv[1], range(1, 3) + range(101, 103) + range(201, 203), False)
    predict(sys.argv[1], range(101, 105), False)

