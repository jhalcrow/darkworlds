import csv
from halomodel import sky_likelihood, load_sky, load_halos, Halo, out_of_bounds, unwarp_sky, XMAX, YMAX, distance
#import numpy as np
import sys
import multiprocessing
import random
import itertools
from math import exp, sqrt

def list_distance(l1, l2):
    d = 0
    for x, y in zip(l1, l2):
        d += (x-y)**2
    return sqrt(d)

def best_n(lls, cutoff=300, n=10):
    lls.sort(key = lambda x: -x[0])
    bests = [lls[0]]
    for ll, c, w in lls:
        too_close = False
        for ll_b, c_b, w_b in bests:
            if list_distance(c, c_b) < cutoff:
                too_close = True
                break
        if not too_close:
            bests.append((ll, c, w))

        if len(bests) >= n:
            return bests

    return bests

def metropolis_hastings(sky, num_halos, h=None, weights=None, runs=30000, sigma=800, sigma_w=20, verbose=True):

    #h = np.random.rand(2*num_halos) * 4200
    if not h:
        h = [random.random()*4200 for i in range(2*num_halos)]

    if not weights:
        weights = [random.random() * 60 + 20 for i in range(num_halos)]

    ll = sky_likelihood(sky, h)
    best_h = h
    best_w = weights
    best_ll = ll
    for i in range(runs):
        if verbose and i % 1000 == 0:
            print i, h, ll
        #hprime = np.random.randn(len(h)) * sigma + h
        hprime = [random.normalvariate(hh, sigma) for hh in h]

        while out_of_bounds(hprime):
            #hprime = np.random.randn(len(h)) * sigma + h
            hprime = [random.normalvariate(hh, sigma) for hh in h]

        weights_prime = [random.normalvariate(ww, sigma_w) for ww in weights]
        while min(weights_prime) < 0 and max(weights_prime) > 200:
            weights_prime = [random.normalvariate(ww, sigma_w) for ww in weights]

        ll_prime = sky_likelihood(sky, hprime, weights_prime)
        if random.random() < exp(ll_prime - ll):
            ll = ll_prime
            h = hprime
            weights = weights_prime

            if ll > best_ll:
                best_h = h
                best_ll = ll
                best_w = weights

    return best_h, best_ll


def brute_force_inner(h):
    return (h, sky_likelihood(SKY, h))

def fix_sky(sky):
    global SKY
    SKY = sky

def brute_force(sky, avoid=None, num_halos=1, dx=50, dw=30, xmin=25, xmax=XMAX, ymin=25, ymax=YMAX, wmin=10, wmax=200):

    x = [float(i*dx) + xmin for i in range(0, xmax/dx)]
    y = [float(i*dx)  + ymin for i in range(0, ymax/dx)]
    weights = [float(i*dw) + wmin for i in range(0, wmax/dw)]
    coords = [[xx, yy] for xx in x for yy in y]
    #weights = [90,]

    ll_best = None
    c_best = None
    w_best = None

    lls = []
    for clist in itertools.combinations(coords, num_halos):
        c = sum(clist, [])
        penalty = 0
        #if avoid is not None:
        #    min_dist = list_distance(avoid, c[:2])
        #    for i in range(2, num_halos):
        #        if

        for w in itertools.combinations_with_replacement(weights, num_halos):
            w = list(w)
            ll = sky_likelihood(sky, c, w)
            lls.append((ll, c, w))
            if ll_best is None or ll > ll_best:
                ll_best = ll
                c_best = c
                w_best = w
                #print 'New best: ', c_best, ll_best, w_best

    return lls
    return list(c_best), ll_best


def optimize_weight(sky, halos, w_min=10, w_max=200, dw=10):

    weights = range(w_min, w_max, dw)
    ll_best = None
    weight = float(w_min)
    w_best = weight
    num_halos = len(halos) / 2

    for weight in itertools.combinations(weights, num_halos):
        ll = sky_likelihood(sky, halos, weight)
        if ll_best < ll or not ll_best:
            ll_best = ll
            w_best = weight
     
    print 'Best weight: ', w_best
    return w_best

def paired_brute_force(sky, num_halos, cutoff=300):
    lls = brute_force(sky, dw=10)
    best = best_n(lls, cutoff, 40)
    best_ll = None
    best_h = None
    best_w = None

    for guesses in itertools.combinations(best, num_halos):
        h =  sum([g[1] for g in guesses], [])
        w = sum([g[2] for g in guesses], [])
        ll = sky_likelihood(sky, h, w)
        if not best_ll or ll > best_ll:
            best_ll = ll
            best_h = h
            best_w = w

    return (best_ll, best_h, best_w)

def alternating_brute_force(sky, num_halos, rounds=10, dx=50):
    halos = [[random.random()*4200, random.random()*4200] for i in range(num_halos)]
    weights = [90 for i in range(num_halos)]

    #halos = [prediction[2*i:2*(i+1)] for i in range(len(prediction)/2)]

    for round_num in range(rounds):
        max_delta = 0
        for h_id in range(num_halos):
            sky_u = sky
            toremove = set(range(num_halos))
            toremove.discard(h_id)
            for i in toremove:
                sky_u = unwarp_sky(sky_u, halos[i], weights[i])
            new_halo, ll = brute_force(sky_u, dx=dx)

            max_delta = max(max_delta, list_distance(new_halo, halos[h_id]))
            halos[h_id] = new_halo
            weights[h_id] = optimize_weight(sky_u, new_halo)

        if max_delta < 50:
            break
    return [hh for h in halos for hh in h], 1



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
                best_w = optimize_weight(sky, sum(true_halos, []))
                print 'True LL', sky_likelihood(sky, sum(true_halos, []), best_w)
                
            if num_halos == 1:
                lls = brute_force(sky, num_halos, dx=50, dw=10)
                lls.sort(key=lambda x: -x[0])
                ll, prediction, _ = lls[0]
            else:
                ll, prediction, w = paired_brute_force(sky, num_halos)
                #prediction, ll = metropolis_hastings(sky, num_halos, sigma=800, runs=100000)
            
            print 'Found: ', prediction, ll
            if not test:
                err = 0
                for n in range(num_halos):
                    err += list_distance(true_halos[n], prediction[2*n:2*(n+1)]) / num_halos

                print 'Error:', err

            halo_pos = [0,0,0,0,0,0]#np.zeros((6))
            for i in range(len(prediction)):
                halo_pos[i] = prediction[i]
            out_line = ['Sky%s' % sky_id, ] + list(halo_pos)
            out.writerow(out_line)
            f.flush()



if __name__ == '__main__':
    predict(sys.argv[1], range(1, 121), True)
    #predict(sys.argv[1], range(1, 11) + range(101, 111) + range(201, 211), False)
    #predict(sys.argv[1], range(101, 105), False)

