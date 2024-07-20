#include "main.h"
#include <fenv.h>
#include <fstream>
#include <iostream>
#include <string>

int read_abcd(const char *filename, double *a, double *b, double *c, double *d)
{

    std::ifstream file(filename);
    std::vector<double> params;
    if (!file.good()) {
        printf("Failed to open file");
        return -2;
    }

    while (file.good() && (params.size() < 4)) {
        std::string line;
        std::getline(file, line);
        int pos_start = -1;
        for (size_t i = 0; i < line.size(); i++) {
            char ch = line[i];
            if (ch == '#') {
                break;
            }
            switch (ch) {
            case ' ':
            case '\n':
                if (pos_start >= 0) {
                    line[i] = '\0';
                    params.push_back(std::stod(&line[0] + pos_start));
                }
                pos_start = -1;
                break;
            default:
                if (pos_start < 0) {
                    pos_start = i;
                }
            }
        }
        if (pos_start >= 0) {
            params.push_back(std::stod(&line[0] + pos_start));
        }
    }

    if (params.size() < 4) {
        printf("Bad file format");
        return -2;
    }

    *a = params[0];
    *b = params[1];
    *c = params[2];
    *d = params[3];
    return 0;
}

int main(int argc, char **argv)
{
    feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW | FE_UNDERFLOW);
    setbuf(stdout, NULL);

    const char *help_message = "0 - change function number\n"
                               "1 - change plot type\n"
                               "2 - increase scale\n"
                               "3 - decrease scale\n"
                               "4 - increase points count\n"
                               "5 - decrease points count\n"
                               "6 - increase middle point error\n"
                               "7 - decrease middle point error\n"
                               "8 - rotate clockwise\n"
                               "9 - rotate counter-clockwise\n";

    double a = -2;
    double b = 2;
    double c = -2;
    double d = 2;
    int n = 4;
    int m = 4;
    int k = 0;
    double eps = 0.1;
    int p = 1;

    if ((argc != 1)) {
        if (argc == 7) {
            int ret = read_abcd(argv[1], &a, &b, &c, &d);
            if (ret != 0) {
                return ret;
            }

            n = std::stoi(argv[2]);
            m = std::stoi(argv[3]);
            k = std::stoi(argv[4]);
            eps = std::stod(argv[5]);
            p = std::stoi(argv[6]);

        } else {
            printf("Wrong number of input arguments");
            return -1;
        }
    }

    printf("INSTRUCTIONS\n");
    printf(help_message);
    printf("END of INSTRUCTIONS\n");
    QApplication app(argc,
                     argv); // создаём приложение, инициализация оконной системы
    QMessageBox msg;

    // создаём виджет класса Scene3D
    Graph sur(NULL, a, b, c, d, n, m, k, eps, p);
    msg.setText(help_message);
    sur.setWindowTitle("Interpol3D"); // название окна
    // размеры (nWidth, nHeight) окна
    sur.resize(500, 500);
    msg.exec();
    sur.show(); // изобразить виджет
    return app.exec();
}
