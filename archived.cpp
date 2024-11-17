#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <set>

double EPSILON = 0.00001;
int NUM_RECTANGLES = 2000;
double LOWER_BOUND = 0;
double UPPER_BOUND = 100;

// generisi tacke
// tacke generisem pomocu random funkcija i ako je degenerisan slucaj odmah odbacim
static void generate_rectangles(int n) {
    std::fstream rectanglesFile("rectangles.txt", std::ios::out | std::ios::trunc);

    double lower_bound = 0;
    double upper_bound = 500;
    std::random_device rd;
    std::mt19937 re(rd());
    std::uniform_real_distribution<double> unifDistrib(lower_bound, upper_bound);

    double ax;
    double ay;
    double bx;
    double by;
    double cx;
    double cy;
    double dx;
    double dy;

    for (int i = 0; i < n; i++) {
        // A nema nikakva ogranicenja
        ax = unifDistrib(re);
        ay = unifDistrib(re);

        // B ne sme biti identicno sa A
        bx = unifDistrib(re);
        by = unifDistrib(re);
        while (std::fabs(bx - ax) < EPSILON && std::fabs(by - ay) < EPSILON) {
            bx = unifDistrib(re);
            by = unifDistrib(re);
        }

        // C ne sme biti kolinearno sa AB
        // ekvivalentno, AC i BC nisu linearno zavisni
        // ako jesu, generisemo ponovo
        cx = unifDistrib(re);
        cy = unifDistrib(re);

        while(true) {
            if (std::fabs((ax - cx)*(by - cy) - (bx - cx)*(ay - cy)) > EPSILON) {
                break;
            }
            cx = unifDistrib(re);
            cy = unifDistrib(re);
        }

        // D ne sme biti kolinearno ni sa AB, ni sa BC ni sa CA
        // ekvivalentno, AD lin nez sa BD, itd.
        dx = unifDistrib(re);
        dy = unifDistrib(re);

        while(true) {
            if (std::fabs((ax - dx)*(by - dy) - (bx - dx)*(ay - dy)) > EPSILON &&
                std::fabs((cx - dx)*(by - dy) - (bx - dx)*(cy - dy)) > EPSILON &&
                std::fabs((ax - dx)*(cy - dy) - (cx - dx)*(ay - dy)) > EPSILON) {
                break;
                }
            dx = unifDistrib(re);
            dy = unifDistrib(re);
        }

        // za sada, cetvorougao ne mora da bude prost! upisimo temena u fajl.
        rectanglesFile << ax << " " << ay << " " << bx << " " << by << " " << cx
            << " " << cy << " " << dx << " " << dy << " ";

        // odmah ubacujem i koordinate sredista stranica
        rectanglesFile << (ax + bx)/2.0 << " " << (ay + by)/2.0 << " "
            << (bx + cx)/2.0 << " " << (by + cy)/2.0 << " "
            << (cx + dx)/2.0 << " " << (cy + dy)/2.0 << " "
            << (dx + ax)/2.0 << " " << (dy + ay)/2.0 << std::endl;

    }
    rectanglesFile.close();
}

static bool has_next_line(std::fstream& file) {
    char c = file.peek();
    return (c != EOF);

}
static void test_parallel() {
    std::fstream rectanglesFile("rectangles.txt", std::ios::in);
    std::fstream hypoFile("hypoParallel.txt", std::ios::out | std::ios::trunc);

    // skup parova duzi == SKUP HIPOTEZA
    std::set<std::vector<int>> skupHipoteza;
    int num_dots = 8;
    int num_of_pairs = 0;

    for(int i = 0; i < num_dots; i++) {
        for(int j = i + 1; j < num_dots; j++) {
            for(int i1 = 0; i1 < num_dots; i1++) {
                for(int j1 = i1 + 1; j1 < num_dots; j1++) {
                    // dodajem ((i, j), (i1, j1)) akko je prvi par leksikografski veci
                    // jer mi ne trebaju duplikati, poredak je nebitan
                    if (i > i1) {
                        skupHipoteza.insert({i, j, i1, j1});
                        num_of_pairs ++;
                    } else if (i == i1) {
                        if (j > j1) {
                            skupHipoteza.insert({i, j, i1, j1});
                            num_of_pairs ++;
                        }
                    }
                }
            }
        }
    }

    while (has_next_line(rectanglesFile)) {
        double vertices[16];
        rectanglesFile >> vertices[0] >> vertices[1] >> vertices[2] >> vertices[3]
            >> vertices[4] >> vertices[5] >> vertices[6] >> vertices[7]
            >> vertices[8] >> vertices[9] >> vertices[10] >> vertices[11]
            >> vertices[12] >> vertices[13] >> vertices[14] >> vertices[15];

        /*
        std::cout << "obradjujem sledeci cetvorougao: \n";
        std::cout << vertices[0] << " " << vertices[1] << " " << vertices[2] << " " << vertices[3] << " "
            << vertices[4] << " " << vertices[5] << " " << vertices[6] << " " << vertices[7] << " "
            << vertices[8] << " " << vertices[9] << " " << vertices[10] << " " << vertices[11] << " "
            << vertices[12] << " " << vertices[13] << " " << vertices[14] << " " << vertices[15] << std::endl;
        */

        for (auto it = skupHipoteza.begin(); it != skupHipoteza.end(); ) {
            // proveriti da li hipoteza vazi za ucitana temena
            // tj proveriti da li za ucitane vertices vazi
            // hipoteza[0] hipoteza[1] || hipoteza[2] hipoteza[3]
            // ako ne vazi, izbaciti hipotezu iz skupa
            int i = it->at(0);
            int j = it->at(1);
            int i1 = it->at(2);
            int j1 = it->at(3);

            if (std::fabs((vertices[2*i] - vertices[2*j])*(vertices[2*i1+1] - vertices[2*j1+1]) -
                        (vertices[2*i1] - vertices[2*j1])*(vertices[2*i+1] - vertices[2*j+1])) < 1) {
                it++;
                // std::cout << i << " " << j << " || " << i1 << " " << j1 << std::endl;
            } else {

                // ako je det veliko, nisu paralelne, pa izbacujemo hipotezu.
                // tu hipotezu onda ne proveravam ni za naredne cetvorouglove.
                it = skupHipoteza.erase(it);
            }
        }
    }
    // pisem one hipoteze koje su prezivele sve cetvorouglove
    for (auto hipoteza : skupHipoteza) {
        hypoFile << hipoteza[0] << " " << hipoteza[1] << " || " << hipoteza[2] << " " << hipoteza[3] << std::endl;
    }

    hypoFile.close();
    rectanglesFile.close();
}

static void test_2x() {
    std::fstream rectanglesFile("rectangles.txt", std::ios::in);
    std::fstream hypoFile("hypo2x.txt", std::ios::out | std::ios::trunc);

    // skup parova duzi == SKUP HIPOTEZA
    std::set<std::vector<int>> skupHipoteza;
    int num_dots = 8;

    for(int i = 0; i < num_dots; i++) {
        for(int j = i + 1; j < num_dots; j++) {
            for(int i1 = 0; i1 < num_dots; i1++) {
                for(int j1 = i1 + 1; j1 < num_dots; j1++) {
                    skupHipoteza.insert({i, j, i1, j1});
                }
            }
        }
    }

    while (has_next_line(rectanglesFile)) {
        double vertices[16];
        rectanglesFile >> vertices[0] >> vertices[1] >> vertices[2] >> vertices[3]
            >> vertices[4] >> vertices[5] >> vertices[6] >> vertices[7]
            >> vertices[8] >> vertices[9] >> vertices[10] >> vertices[11]
            >> vertices[12] >> vertices[13] >> vertices[14] >> vertices[15];

        for (auto it = skupHipoteza.begin(); it != skupHipoteza.end(); ) {
            // proveriti da li hipoteza vazi za ucitana temena
            // ako ne vazi, izbaciti hipotezu iz skupa
            int i = it->at(0);
            int j = it->at(1);
            int i1 = it->at(2);
            int j1 = it->at(3);

            if (std::fabs(
                    (vertices[2*i] - vertices[2*j])*(vertices[2*i] - vertices[2*j]) +
                    (vertices[2*i+1] - vertices[2*j+1])*(vertices[2*i+1] - vertices[2*j+1]) -
                    4*(vertices[2*i1] - vertices[2*j1])*(vertices[2*i1] - vertices[2*j1]) -
                    4*(vertices[2*i1+1] - vertices[2*j1+1])*(vertices[2*i1+1] - vertices[2*j1+1])
                ) < 10)
            {
                it++;
            } else {
                it = skupHipoteza.erase(it);
            }
        }
    }
    // pisem one hipoteze koje su prezivele sve cetvorouglove
    for (auto hipoteza : skupHipoteza) {
        hypoFile << hipoteza[0] << " " << hipoteza[1] << " = 2 x " << hipoteza[2] << " " << hipoteza[3] << std::endl;
    }

    hypoFile.close();
    rectanglesFile.close();
}

static void test_3x() {
    std::fstream rectanglesFile("rectangles.txt", std::ios::in);
    std::fstream hypoFile("hypo3x.txt", std::ios::out | std::ios::trunc);

    // skup parova duzi == SKUP HIPOTEZA
    std::set<std::vector<int>> skupHipoteza;
    int num_dots = 8;

    for(int i = 0; i < num_dots; i++) {
        for(int j = i + 1; j < num_dots; j++) {
            for(int i1 = 0; i1 < num_dots; i1++) {
                for(int j1 = i1 + 1; j1 < num_dots; j1++) {
                    skupHipoteza.insert({i, j, i1, j1});
                }
            }
        }
    }

    while (has_next_line(rectanglesFile)) {
        double vertices[16];
        rectanglesFile >> vertices[0] >> vertices[1] >> vertices[2] >> vertices[3]
            >> vertices[4] >> vertices[5] >> vertices[6] >> vertices[7]
            >> vertices[8] >> vertices[9] >> vertices[10] >> vertices[11]
            >> vertices[12] >> vertices[13] >> vertices[14] >> vertices[15];

        for (auto it = skupHipoteza.begin(); it != skupHipoteza.end(); ) {
            // proveriti da li hipoteza vazi za ucitana temena
            // ako ne vazi, izbaciti hipotezu iz skupa
            int i = it->at(0);
            int j = it->at(1);
            int i1 = it->at(2);
            int j1 = it->at(3);

            if (std::fabs(
                    (vertices[2*i] - vertices[2*j])*(vertices[2*i] - vertices[2*j]) +
                    (vertices[2*i+1] - vertices[2*j+1])*(vertices[2*i+1] - vertices[2*j+1]) -
                    9*(vertices[2*i1] - vertices[2*j1])*(vertices[2*i1] - vertices[2*j1]) -
                    9*(vertices[2*i1+1] - vertices[2*j1+1])*(vertices[2*i1+1] - vertices[2*j1+1])
                ) < 10)
            {
                it++;
            } else {
                it = skupHipoteza.erase(it);
            }
        }
    }
    // pisem one hipoteze koje su prezivele sve cetvorouglove
    for (auto hipoteza : skupHipoteza) {
        hypoFile << hipoteza[0] << " " << hipoteza[1] << " = 2 x " << hipoteza[2] << " " << hipoteza[3] << std::endl;
    }

    hypoFile.close();
    rectanglesFile.close();
}