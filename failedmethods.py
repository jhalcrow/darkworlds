# Old models


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
