#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <set>

double EPSILON = 0.00001;
int NUM_RECTANGLES = 2000;
double LOWER_BOUND = 0;
double UPPER_BOUND = 100;

class Point {
public:
    Point() {
        std::random_device rd;
        std::mt19937 re(rd());
        std::uniform_real_distribution<double> unifDistrib(LOWER_BOUND, UPPER_BOUND);
        x = unifDistrib(re);
        y = unifDistrib(re);
    }
    Point(double _x, double _y): x(_x), y(_y) {}
    bool operator==(const Point &other) const {
        return x - other.x < EPSILON && y - other.y < EPSILON; //
    }
    void generate_again() {
        std::random_device rd;
        std::mt19937 re(rd());
        std::uniform_real_distribution<double> unifDistrib(LOWER_BOUND, UPPER_BOUND);
        x = unifDistrib(re);
        y = unifDistrib(re);
    }
    static bool pointsCollinear(const Point &a, const Point &b, const Point &c) {
        return std::fabs((a.x - c.x)*(b.y - c.y) - (b.x - c.x)*(a.y - c.y)) < EPSILON;
    }
    static bool rectangleDegenerate(const Point &a, const Point &b, const Point &c, const Point &d) {
        return !pointsCollinear(a, b, c) && !pointsCollinear(a, b, d) && !pointsCollinear(b, c, d);
    }
    double x;
    double y;
};

class Rectangle {
    public:
    Rectangle(Point a, Point b, Point c, Point d) {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        this->ab = Point((a.x + b.x)/2, (a.y+b.y)/2);
        this->bc = Point((b.x + c.x)/2, (b.y+c.y)/2);
        this->cd = Point((c.x + d.x)/2, (c.y+d.y)/2);
        this->da = Point((a.x + d.x)/2, (a.y+d.y)/2);
    }
    // todo: kako da se u fji osiguram da je vertices bar 16 mesta dug?
    void returnVerticeArray(double* vertices) const {
        vertices[0] = a.x;
        vertices[1] = a.y;
        vertices[2] = b.x;
        vertices[3] = b.y;
        vertices[4] = c.x;
        vertices[5] = c.y;
        vertices[6] = d.x;
        vertices[7] = d.y;
        vertices[8] = ab.x;
        vertices[9] = ab.y;
        vertices[10] = bc.x;
        vertices[11] = bc.y;
        vertices[12] = cd.x;
        vertices[13] = cd.y;
        vertices[14] = da.x;
        vertices[15] = da.y;
    }
    Point a;
    Point b;
    Point c;
    Point d;
    Point ab;
    Point bc;
    Point cd;
    Point da;
};

static void test_parallel(const Rectangle &rec, std::set<std::vector<int>> &hypothesis_set) {
    double vertices[16];
    rec.returnVerticeArray(vertices);

    for (auto it = hypothesis_set.begin(); it != hypothesis_set.end(); ) {
        // proveriti da li hipoteza vazi za ucitana temena
        // tj proveriti da li za ucitane vertices vazi
        // hipoteza[0] hipoteza[1] || hipoteza[2] hipoteza[3]
        // ako ne vazi, izbaciti hipotezu iz skupa
        int i = it->at(0);
        int j = it->at(1);
        int i1 = it->at(2);
        int j1 = it->at(3);
        // todo: ovde je EPSILON 1 i ne radi mi za manje EPSILON ??
        if (std::fabs((vertices[2*i] - vertices[2*j])*(vertices[2*i1+1] - vertices[2*j1+1]) -
                    (vertices[2*i1] - vertices[2*j1])*(vertices[2*i+1] - vertices[2*j+1])) < 1) {
            it++;
                    } else {
                        // ako je det veliko, nisu paralelne, pa izbacujemo hipotezu.
                        // tu hipotezu onda ne proveravam ni za naredne cetvorouglove.
                        it = hypothesis_set.erase(it);
                    }
    }
}

static void generate_hypothesis(std::set<std::vector<int>> &hypothesis_set, bool symmetric = false) {
    int num_dots = 8;
    if (symmetric) {
        for(int i = 0; i < num_dots; i++) {
            for(int j = i + 1; j < num_dots; j++) {
                for(int i1 = 0; i1 < num_dots; i1++) {
                    for(int j1 = i1 + 1; j1 < num_dots; j1++) {
                        // dodajem ((i, j), (i1, j1)) akko je prvi par leksikografski veci
                        if (i > i1) {
                            hypothesis_set.insert({i, j, i1, j1});
                        } else if (i == i1) {
                            if (j > j1) {
                                hypothesis_set.insert({i, j, i1, j1});
                            }
                        }
                    }
                }
            }
        }
        return;
    }
    // ako nije simetricno tvrdjenje, sve dodajemo
    for(int i = 0; i < num_dots; i++) {
        for(int j = i + 1; j < num_dots; j++) {
            for(int i1 = 0; i1 < num_dots; i1++) {
                for(int j1 = i1 + 1; j1 < num_dots; j1++) {
                    hypothesis_set.insert({i, j, i1, j1});
                }
            }
        }
    }
}

int main() {
    Point a, b, c, d;
    while (Point::rectangleDegenerate(a, b, c, d)) {
        a.generate_again();
        b.generate_again();
        c.generate_again();
        d.generate_again();
    }
    Rectangle rectangle(a, b, c, d);



    return 0;
}
