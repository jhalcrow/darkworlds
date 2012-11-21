
struct Position {
    double x;
    double y;

    Position(double x_, double y_) : x(x_), y(y_) {};
};

struct Ellipticity {
    double e1;
    double e2;

    Ellipticity(double e1_, double e2_) : e1(e1_), e2(e2_) {};
};

struct Galaxy {
    Position pos;
    Ellipticity e;

    Galaxy(Position pos_, Ellipticity e_) : pos(pos_), e(e_) {};
};

struct Halo {
    Position pos;
    double weight;

    Halo(Position pos_, double weight_) : pos(pos_), weight(weight_) {};
};

struct ScoredHalo {
    Halo halo;
    double score;

    ScoredHalo(Halo halo_, double score_) : halo(halo_), score(score_) {};
};

struct Range {
    float min;
    float max;
    float dx;
};


class Sky {
private:


public:
    std::vector<Galaxy> sky;

    std::map<int, float> prior;

    Sky(int id, bool test, const std::map<int, float>& prior_);

    void unwarp_sky(std::vector<Halo> halos);

    double likelihood(const std::vector<Halo>& halos) const;
    double likelihood(const Halo& h) const;
};
