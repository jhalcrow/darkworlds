#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <queue>
#include <algorithm>
#include <map>
#include <sstream>
#include <math.h>

#include "halomodel2.h"

using namespace std;

float dist(const Halo& h1, const Halo& h2) {
    return sqrt(pow((h1.pos.x - h2.pos.x), 2) + pow((h1.pos.y - h2.pos.y), 2));
}

float dist(const ScoredHalo& h1, const ScoredHalo& h2) {
    return dist(h1.halo, h2.halo);
}

vector<float> range(float min, float max, float dx) {
    int N = static_cast<int>((max - min) / dx);
    vector<float> vals(N);
    float cur = min;
    for(int i = 0; i < N; ++i) {
        vals.push_back(cur);
        cur += dx;
    }

    return vals;
}

Ellipticity unwarp(const Galaxy& g, const vector<Halo> halos) {
    Ellipticity e(g.e.e1, g.e.e2);

    for(auto h: halos) {
        double dx = g.pos.x - h.pos.x;
        double dy = g.pos.y - h.pos.y;
        double r3 = pow(dx*dx + dy*dy, 1.5);
        e.e1 += h.weight * (dx*dx - dy*dy) / r3;
        e.e2 += h.weight * 2 * dx * dy / r3;
    }

    return e;
}

Ellipticity unwarp(const Galaxy& g, const Halo& h) {
    Ellipticity e(g.e.e1, g.e.e2);

    double dx = g.pos.x - h.pos.x;
    double dy = g.pos.y - h.pos.y;
    double r3 = pow(dx*dx + dy*dy, 1.5);
    e.e1 += h.weight * (dx*dx - dy*dy) / r3;
    e.e2 += h.weight * 2 * dx * dy / r3;


    return e;
}

Sky::Sky(int id, bool test, const std::map<int,float>& prior_) : sky(), prior(prior_) {
    ifstream sky_file;
    string filename;
    stringstream ssfilename( "data/", stringstream::in | stringstream::out); 
    if(test) {
        ssfilename << "data/Test_Skies/Test_Sky" << id << ".csv";
    } else {
        ssfilename << "data/Train_Skies/Training_Sky" << id << ".csv";
    }
    ssfilename >> filename;
    sky_file.open(filename.c_str());

    if (sky_file) {
        string line;

        getline(sky_file, line); // Discard header

        while(getline(sky_file, line)) {

            stringstream ss( stringstream::out | stringstream::in );
            ss << line;
            
            char c = 100;
            //Eat the galaxy name
            int i = 0; 
            while(c != ',') {
                ss >> c;
            }
            double x,y,e1,e2;
            ss >> x >> c >> y >> c >> e1 >> c >> e2;
            Position pos(x, y);
            Ellipticity e(e1, e2);
            Galaxy g(pos, e);

            sky.push_back(g);
        }

    } else {
        cerr << "Error reading file " << filename << endl;
        exit(1);
    }
    sky_file.close();
}

double Sky::likelihood(const vector<Halo>& halos) const {
    double ll = 0.0;
    for (auto h: halos) {
        ll += this->prior.at(static_cast<int>(h.weight));
    }

    for (auto g : sky) {
        Ellipticity e = unwarp(g, halos);
        ll -= e.e1 * e.e1 + e.e2 * e.e2;
    }

    return ll;
}

double Sky::likelihood(const Halo& h) const {
    double ll = this->prior.at(static_cast<int>(h.weight));

    for (auto g : sky) {
        Ellipticity e = unwarp(g, h);
        ll -= e.e1 * e.e1 + e.e2 * e.e2;
    }

    return ll;
}

vector<ScoredHalo> scan_pos(const Sky& sky, const int n, const Range xrange, const Range yrange, const Range wrange) {
    const double dist_cutoff = 300.0;

    const int nx = static_cast<int>((xrange.max - xrange.min) / xrange.dx);
    const int ny = static_cast<int>((yrange.max - yrange.min) / yrange.dx);
    const int nw = static_cast<int>((wrange.max - wrange.min) / wrange.dx);

    auto comp = []( ScoredHalo h1, ScoredHalo h2 ) { return h1.score < h2.score; };
    priority_queue<ScoredHalo, vector<ScoredHalo>, decltype(comp)> scores( comp );

    Halo h(Position(xrange.min, yrange.min), wrange.min);

    for(int ix=0; ix < nx; ++ix) {
        h.pos.x += xrange.dx;
        h.pos.y = yrange.min;
        for(int iy=0; iy < ny; ++iy) {
            h.pos.y += yrange.dx;
            h.weight = wrange.min;
            for(int iw=0; iw < nw; ++iw) {
                h.weight += wrange.dx;
                double ll = sky.likelihood(h);
                ScoredHalo newh(Halo(Position(h.pos.x, h.pos.y), h.weight), ll);
                scores.push(newh);
            }
        }
    }

    vector<ScoredHalo> best;
    while(scores.size() > 0) {
        if(best.size() >= n) {
            break;
        }
        ScoredHalo hscore = scores.top();
        scores.pop();

        bool tooclose = false;

        for(auto h: best) {
            if (dist(h.halo, hscore.halo) < dist_cutoff) {
                tooclose = true;
                break;
            }
        }
        if(!tooclose) {
            best.push_back(hscore);
        }
    }
    auto comp2 = []( ScoredHalo h1, ScoredHalo h2 ) { return h1.score > h2.score; };
    sort(best.begin(), best.end(), comp2);

    return best;
}

map<int, float> load_mass_prior(string filename) {
    ifstream mass_file;
    mass_file.open(filename.c_str());
    map<int, float> prior;
    if (mass_file) {
        string line;
        while(getline(mass_file, line)) {

            stringstream ss( stringstream::out | stringstream::in );
            ss << line;
            
            char c = 100;
            int w;
            float pdf;

            ss >> w >> c >> pdf;

            prior[w] = pdf;
        }
    }
    return prior;
}

vector<int> load_halo_counts(const bool test) {
    string filename;
    vector<int> counts;

    if(test) {
        filename = "data/Test_haloCounts.csv";
    } else {
        filename = "data/Training_halos.csv";
    }
    ifstream count_file;
    count_file.open(filename.c_str());
    string line;
    getline(count_file, line); // Eat header

    while(getline(count_file, line)) {
        stringstream ss( stringstream::out | stringstream::in );
        ss << line;
        char c = 100;
        int ct;
        while(c != ',') {
            ss >> c;
        }
        ss >> ct;
        counts.push_back(ct);
    }
    return counts;
}

vector<Halo> best_n(const Sky& sky, const vector<ScoredHalo>& scores, const int n) {

    vector<Halo> best;
    double best_score=-1000000000000;

    if (n == 1) {
        best.push_back(scores[0].halo);
        return best;
        // This is awful 
    } else if (n == 2) {
        for(int i = 0; i < scores.size(); ++i) {
            for(int j = i + 1; j < scores.size(); ++j) {
                vector<Halo> halo_set;
                halo_set.push_back(scores[i].halo);
                halo_set.push_back(scores[j].halo);
                double score = sky.likelihood(halo_set);
                if (score > best_score) {
                    best_score = score;
                    best = halo_set;
                }

            }
        }
    }
    if (n == 3) {
        for(int i = 0; i < scores.size(); ++i) {
            for(int j = i + 1; j < scores.size(); ++j) {
                for(int k = j + 1; k < scores.size(); ++k) {
                    vector<Halo> halo_set;
                    halo_set.push_back(scores[i].halo);
                    halo_set.push_back(scores[j].halo);
                    halo_set.push_back(scores[k].halo);
                    double score = sky.likelihood(halo_set);
                    if (score > best_score) {
                        best_score = score;
                        best = halo_set;
                    }
                }
            }
        }
    }

    return best;

}


int main() {

    const bool test = false;
    map<int, float> prior = load_mass_prior("mass_prior.csv");
    vector<int> num_halos = load_halo_counts(test);

    Range xrange = {0, 4200, 50};
    Range yrange = {0, 4200, 50};
    Range wrange = {10, 240, 20};

    for(int i = 1; i < 2; ++i) {

        Sky sky(i, test, prior);
        int halos = num_halos[i];
        cout << "Sky" << i;

        vector<ScoredHalo> scores = scan_pos(sky, 40, xrange, yrange, wrange);
        
        vector<Halo> best_scores = best_n(sky, scores, halos);
        for(auto halo: best_scores) {
            cout << "," << halo.pos.x << "," << halo.pos.y;
        }
        
        while(halos < 3) {
            cout << ",0,0";
            halos++;
        } 
        cout << endl;
    }

}