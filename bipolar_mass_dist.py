from predict import *

with open('masses_all.csv', 'w') as f:
    for sky_id in range(1, 300):
        print sky_id
        sky = load_sky(sky_id)
        halos = load_halos()[sky_id - 1]
        for halo in halos:
            bp = best_bipolar(sky, halo, wmax=150, dw=5, dx=50)
            print bp
            f.write('%s,%s,%s\n' % (list_distance(*bp.halos), bp.weights[0], bp.weights[1]))
            f.flush()

