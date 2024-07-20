#include "surface.h"
#include "array.h"
#include "function.h"
#include <fstream>

#define DEVELOP

GLfloat VertexArray1[DEF_X * DEF_Y][3];
GLfloat VertexArray2[DEF_X * DEF_Y][3];
GLfloat VertexArray3[DEF_X * DEF_Y][3];
GLfloat ColorArray1[DEF_X * DEF_Y][3];
GLfloat ColorArray2[DEF_X * DEF_Y][3];
GLfloat ColorArray3[DEF_X * DEF_Y][3];
GLuint IndexArray[(DEF_X - 1) * (DEF_Y - 1)][3];

struct FuncInfo {
    double (*f)(double, double);
    double (*dfdx)(double, double);
    double (*dfdy)(double, double);
    double (*d2fdxdy)(double, double);
    int min_scale;
    const char *text;
};

const float MAX_SCALE_Z = 1e3;

const FuncInfo FUNCTIONS[] = {
    {f1, df1_dx, df1_dy, d2f1_dxdy, 0, "f(x, y)=1"},
    {f2, df2_dx, df2_dy, d2f2_dxdy, 0, "f(x, y)=x"},
    {f3, df3_dx, df3_dy, d2f3_dxdy, 0, "f(x, y)=y"},
    {f4, df4_dx, df4_dy, d2f4_dxdy, 0, "f(x, y)=x+y"},
    {f5, df5_dx, df5_dy, d2f5_dxdy, -24, "f(x, y)=sqrt(x^2+y^2)"},
    {f6, df6_dx, df6_dy, d2f6_dxdy, -19, "f(x, y)=x^2+y^2"},
    {f7, df7_dx, df7_dy, d2f7_dxdy, -3, "f(x, y)=exp(x^2-y^2)"},
    {f8, df8_dx, df8_dy, d2f8_dxdy, 0, "f(x, y)=1/((25(x^2+y^2)+1)"},
};

// double (* FUNCTIONS [])(double, double) = {f1, f2, f3, f4, f5, f6, f7, f8};
const int FUNC_COUNT = 8;

Graph::Graph(QWidget *parent, int a_, int b_, int c_, int d_, int n_, int m_,
             int k_, double eps_, int p_)
    : QGLWidget(parent), a(a_), b(b_), c(c_), d(d_), n(n_), m(m_), mf(nullptr),
      xRot(-90), yRot(0), zRot(0), zTra(0), scale_xy(1), scale_z(1), eps(eps_)
{
    maxpog = 0;
    chh = 0;
    graph_id = 0;
    dif_id = 0;
    graphf_id = 1;
    ch = 0;
    error_factor = 0;
    f_max_abs = 0.;
    f_interp_abs_max = 0.;
    f_error_abs_max = 0.;
    f_max_val = &f_max_abs;
    scale_xy = (GLfloat)getScaleXY();
    scale = 0;
    scale_z = 0.8;
    graph_mode = 0;

    func_id = k_;
    thread_count = p_;

    validateParameters();
    defineVertexShape();
}

/*virtual*/ void Graph::initializeGL() // инициализация
{
    qglClearColor(Qt::white); // цвет для очистки буфера изображения - здесь
    // просто фон окна
    glEnable(GL_DEPTH_TEST); // устанавливает режим проверки глубины пикселей
    glShadeModel(GL_FLAT); // отключает режим сглаживания цветов

    getVertexArray(); // определить массив вершин
    getColorArray(); // определить массив цветов вершин
    getIndexArray(); // определить массив индексов вершин

    glEnableClientState(GL_VERTEX_ARRAY); // активизация массива вершин
    glEnableClientState(GL_COLOR_ARRAY); // активизация массива цветов вершин
}

/*virtual*/ void Graph::resizeGL(int nWidth, int nHeight) // окно виджета
{
    glMatrixMode(GL_PROJECTION); // устанавливает текущей проекционную матрицу
    glLoadIdentity(); // присваивает проекционной матрице единичную матрицу

    // отношение высоты окна виджета к его ширине
    GLfloat ratio = (GLfloat)nHeight / (GLfloat)nWidth;

    // мировое окно
    if (nWidth >= nHeight)
        glOrtho(-1.0 / ratio, 1.0 / ratio, -1.0, 1.0, -10.0,
                1.0); // параметры видимости ортогональной проекции
    else
        glOrtho(-1.0, 1.0, -1.0 * ratio, 1.0 * ratio, -10.0,
                1.0); // параметры видимости ортогональной проекции
    // плоскости отсечения (левая, правая, верхняя, нижняя, передняя, задняя)

    // glFrustum(-1.0, 1.0, -1.0, 1.0, 1.0, 10.0); // параметры видимости
    // перспективной проекции плоскости отсечения (левая, правая, верхняя,
    // нижняя, ближняя, дальняя)

    // поле просмотра
    glViewport(0, 0, (GLint)nWidth, (GLint)nHeight);
}

/*virtual*/ void Graph::paintGL() // рисование
{
    // glClear(GL_COLOR_BUFFER_BIT); // окно виджета очищается текущим цветом
    // очистки
    glClear(GL_COLOR_BUFFER_BIT |
            GL_DEPTH_BUFFER_BIT); // очистка буфера изображения и глубины

    glMatrixMode(GL_MODELVIEW); // устанавливает положение и
    // ориентацию матрице моделирования
    glLoadIdentity(); // загружает единичную матрицу моделирования

    // последовательные преобразования
    glScalef(scale_xy, scale_z, scale_xy); // масштабирование
    // трансляция
    glTranslatef(0.0f, zTra, 0.0f);
    // поворот вокруг оси X
    glRotatef(xRot, 1.0f, 0.0f, 0.0f);
    // поворот вокруг оси Y
    glRotatef(yRot, 0.0f, 1.0f, 0.0f);
    // поворот вокруг оси Z
    glRotatef(zRot, 0.0f, 0.0f, 1.0f);

    // рисование осей координат
    drawAxis();
    drawFigure1(); // нарисовать фигуру
    drawFigure2();
    drawFigure3();
    drawText();
}

void Graph::drawText()
{
    int textX = 5;
    int textY = 20;
    char buffer[64];
    const FuncInfo *info = FUNCTIONS + func_id;
    sprintf(buffer, "k=%d: %s", func_id, info->text);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "a = %.3f", a);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "b = %.3f", b);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "c = %.3f", c);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "d = %.3f", d);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "N = %d", n);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "M = %d", m);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "Rotation angle = %.1f\u00B0", zRot);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "F_max = %.3f", *f_max_val);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "p = %d", error_factor);
    renderText(textX, textY, QString(buffer));
    textY += 20;
    sprintf(buffer, "scale = %d", scale);
    renderText(textX, textY, QString(buffer));
}

/*virtual*/ void Graph::mousePressEvent(QMouseEvent *pe) // нажатие клавиши мыши
{
    // при нажатии пользователем кнопки мыши переменной ptrMousePosition будет
    // присвоена координата указателя мыши
    ptrMousePosition = pe->pos();

    // ptrMousePosition = (*pe).pos(); // можно и так написать
}

/*virtual*/ void
Graph::mouseReleaseEvent(QMouseEvent *pe) // отжатие клавиши мыши
{
    // некоторые функции, которые должны выполняться при отжатии клавиши мыши
    ptrMousePosition = pe->pos();
}

/*virtual*/ void
Graph::mouseMoveEvent(QMouseEvent *pe) // изменение положения стрелки мыши
{
    int y = pe->y() - ptrMousePosition.y();
    // вычисление углов поворота
    xRot += 180 / scale_xy * (GLfloat)(y) / height();
    int x = pe->x() - ptrMousePosition.x();
    zRot += 180 / scale_xy * (GLfloat)(x) / width();

    ptrMousePosition = pe->pos();
    normalizeRotationAngle();

    updateGL(); // обновление изображения
}

/*virtual*/ void Graph::wheelEvent(QWheelEvent *pe) // вращение колёсика мыши
{
    if ((pe->delta()) > 0)
        scale_plus();
    else if ((pe->delta()) < 0)
        scale_minus();

    updateGL(); // обновление изображения
}

void Graph::changeViewType()
{
    graph_mode = (graph_mode + 1) % 3;
    graph_id = 0;
    graphf_id = 0;
    dif_id = 0;
    switch (graph_mode) {
    case 0:
        graphf_id = 1;
        f_max_val = &f_max_abs;
        break;
    case 1:
        graph_id = 1;
        f_max_val = &f_interp_abs_max;
        break;
    case 2:
        dif_id = 1;
        f_max_val = &f_error_abs_max;
        break;
    }
}

/*virtual*/ void
Graph::keyPressEvent(QKeyEvent *pe) // нажатие определенной клавиши
{
    switch (pe->key()) {
    case Qt::Key_1:
        changeViewType(); // сменить график
        // для отображения (функция, аппроксимация, погрешность)
        break;
    case Qt::Key_2:
        if (scale_axes(1)) {
            scale_plus(2.f);
            getVertexArray();
        }
        break;
    case Qt::Key_3:
        if (scale_axes(-1)) {
            scale_minus(2.f);
            getVertexArray();
        }
        break;
    case Qt::Key_Q:
        chang_cd_m();
        chang_matr();
        getVertexArray();
        break;
    case Qt::Key_W:
        chang_cd_p();
        chang_matr();
        getVertexArray();
        break;
    case Qt::Key_4:
        changePointsNum(true);
        getVertexArray();
        break;
    case Qt::Key_5:
        changePointsNum(false);
        getVertexArray();
        break;
    case Qt::Key_6:
        error_factor -= 1;
        // TODO: update error in the middle point only
        //        addError();
        getVertexArray();
        break;
    case Qt::Key_7:
        error_factor += 1;
        // TODO: update error in the middle point only
        //        addError();
        getVertexArray();
        break;
    case Qt::Key_0:
        chang_graphfunction(); // заменить отображаемую функцию
        getVertexArray();
        break;
    case Qt::Key_Plus:
        scale_plus(); // приблизить сцену
        break;
    case Qt::Key_C:
        chang_ch();
        getVertexArray(); // приблизить сцену
        break;
    case Qt::Key_Equal:
        scale_plus(); // приблизить сцену
        break;

    case Qt::Key_Minus:
        scale_minus(); // удалиться от сцены
        break;

    case Qt::Key_Up:
        rotate_up(); // повернуть сцену вверх
        break;

    case Qt::Key_Down:
        rotate_down(); // повернуть сцену вниз
        break;

    case Qt::Key_Left:
        rotate_left(); // повернуть сцену влево
        break;

    case Qt::Key_8:
        rotate_left(30.f); // повернуть сцену влево
        break;

    case Qt::Key_9:
        rotate_right(30.f); // повернуть сцену вправо
        break;

    case Qt::Key_Right:
        rotate_right(); // повернуть сцену вправо
        break;

    case Qt::Key_Z:
        translate_down(); // транслировать сцену вниз
        break;

    case Qt::Key_X:
        translate_up(); // транслировать сцену вверх
        break;

    case Qt::Key_Space: // клавиша пробела
        defaultScene(); // возвращение значений по умолчанию
        break;

    case Qt::Key_Escape: // клавиша "эскейп"
        // завершает приложение
        this->close();
        break;
    }

    normalizeRotationAngle();
    scale_z = getScaleZ();
    printf("|F_max| = %f\n", *f_max_val);

    updateGL(); // обновление изображения
}

void Graph::chang_matr()
{
    //    if (!Init(a, b, c, d))
    //    {
    //        printf("Not enough memory.\n");
    //        Finalize();
    //    }
    // Input(ch);
    // Calc(ch);
}

void Graph::chang_ch()
{
    chh = (chh + 1) % 2;
    ch = (ch + 1) % 2;
    printf("n=%d\n", n);
    printf("ch=%d\n", ch);
    printf("maxpogr=%e\n", maxpog);
}

void Graph::chang_ab_p()
{
    a -= 3;
    b += 3;

    printf("[a,b]=[%f,%f]\n", a, b);
    printf("maxpogr=%e\n", maxpog);
}

void Graph::chang_cd_p()
{
    c -= 3;
    d += 3;

    printf("[c,d]=[%f,%f]\n", c, d);
    printf("maxpogr=%e\n", maxpog);
}

void Graph::chang_ab_m()
{
    if ((a < -3) && (b > 3)) {
        a += 3;
        b -= 3;
        printf("[a,b]=[%f,%f]\n", a, b);
        printf("maxpogr=%e\n", maxpog);
    }
}
void Graph::chang_cd_m()
{
    if (c < -3)
        c += 3;
    if (d > 3)
        d -= 3;

    printf("[c,d]=[%f,%f]\n", c, d);
    printf("maxpogr=%e\n", maxpog);
}

void Graph::chang_graphfunction()
{
    func_id = (func_id + 1) % FUNC_COUNT;

    int min_scale = FUNCTIONS[func_id].min_scale;
    if ((min_scale != 0) && (scale <= min_scale)) {
        scale = min_scale + 1;
        scale_xy = getScaleXY() * pow(2, scale);
    }
}

void Graph::chang_graph()
{
    graph_id = (graph_id + 1) % 2;
}

void Graph::chang_dif()
{
    dif_id = (dif_id + 1) % 2;
}

void Graph::changePointsNum(bool increase)
{
    if (increase) {
        n *= 2;
        m *= 2;
    } else {
        n /= 2;
        m /= 2;
    }

    if ((n < 4)) {
        n = 4;
    }

    if ((m < 4)) {
        m = 4;
    }

    defineVertexShape();
}

void Graph::scale_plus(float k) // приблизить сцену
{
    scale_xy = scale_xy * k;
}

void Graph::scale_minus(float k) // удалиться от сцены
{
    scale_xy = scale_xy / k;
}

bool Graph::scale_axes(int k)
{
    if (k < 0) {
        int new_scale = scale + k;
        int min_scale = FUNCTIONS[func_id].min_scale;
        if ((min_scale < 0) && (new_scale <= min_scale)) {
            return false;
        }
    }

    scale += k;
    return true;
}

void Graph::rotate_up() // повернуть сцену вверх
{
    xRot += 1.0;
}

void Graph::rotate_down() // повернуть сцену вниз
{
    xRot -= 1.0;
}

void Graph::rotate_left(float angle) // повернуть сцену влево
{
    zRot += angle;
}

void Graph::rotate_right(float angle) // повернуть сцену вправо
{
    zRot -= angle;
}

void Graph::translate_down() // транслировать сцену вниз
{
    zTra -= 0.05;
}

void Graph::translate_up() // транслировать сцену вверх
{
    zTra += 0.05;
}

void Graph::defaultScene() // наблюдение сцены по умолчанию
{
    xRot = -90;
    yRot = 0;
    zRot = 0;
    zTra = 0;
    scale_xy = 1;
    scale_z = 1;
}

void Graph::drawAxis() // построить оси координат
{
    glLineWidth(3.0f); // устанавливаю ширину линии приближённо в пикселях
    // до вызова команды ширина равна 1 пикселю по умолчанию

    glColor4f(1.00f, 0.00f, 0.00f,
              1.0f); // устанавливается цвет последующих примитивов
    // ось x красного цвета
    // построение линии
    glBegin(GL_LINES);
    // первая точка
    glVertex3f(1.0f, 0.0f, 0.0f);
    glVertex3f(-1.0f, 0.0f, 0.0f); // вторая точка
    glEnd();

    QColor halfGreen(0, 128, 0, 255);
    qglColor(halfGreen);
    glBegin(GL_LINES);
    // ось y зеленого цвета
    glVertex3f(0.0f, 1.0f, 0.0f);
    glVertex3f(0.0f, -1.0f, 0.0f);

    glColor4f(0.00f, 0.00f, 1.00f, 1.0f);
    // ось z синего цвета
    glVertex3f(0.0f, 0.0f, 1.0f);
    glVertex3f(0.0f, 0.0f, -1.0f);
    glEnd();
}

void Graph::addError()
{
    double mid_point_error = error_factor * 0.1f * f_max_abs;
    int i = nx / 2;
    int j = ny / 2;

    double dx = getDX();
    double dy = getDY();
    const FuncInfo *finfo = &FUNCTIONS[func_id];
    double f_val = finfo->f(a + i * dx, c + j * dy);
    VertexArray1[i + nx * j][2] = f_val + mid_point_error;
}

double Graph::evaluateFunction(double x1, double y1, int i, int j, double dx,
                               double dy)
{
    const FuncInfo *finfo = &FUNCTIONS[func_id];
    double f_val = finfo->f(x1 + i * dx, y1 + j * dy);
    bool is_error_point = false;
    if (nx > n) {
        is_error_point =
            (i == (n / 2) * (nx / (n - 1))) && (j == (m / 2) * (ny / (m - 1)));
    } else {
        is_error_point = (i == nx / 2) && (j == ny / 2);
    }
    if (is_error_point) {
        double mid_point_error = error_factor * 0.1f * f_max_abs;
        f_val += mid_point_error;
    }
    return f_val;
}

double Graph::getDX()
{

    double dx = (b - a) / (nx - 1);
    return dx;
}

double Graph::getDY()
{

    double dy = (d - c) / (ny - 1);
    return dy;
}

double Graph::getScaleZ()
{
    double scale_z = 0.8 / *f_max_val;
    if (scale_z > MAX_SCALE_Z) {
        scale_z = 1;
    }
    return scale_z;
}

double Graph::getScaleXY()
{
    double eps = 1e-6;
    double scale_x = 1 / (b - a + eps);
    double scale_y = 1 / (d - c + eps);
    double scale = scale_x;
    if (scale_y < scale_x) {
        scale = scale_y;
    }
    return 1.5 * scale;
}

void Graph::defineVertexShape()
{
    int target_size = 64;
    //    nx = 64; ny = 64;
    if (n > target_size) {
        int scale_x = (n / target_size);
        int scale_y = (m / target_size);

        nx = n / scale_x;
        ny = m / scale_y;

    } else {
        nx = (target_size / (n - 1)) * (n - 1);
        ny = (target_size / (m - 1)) * (m - 1);
    }
    getIndexArray();
}

void Graph::validateParameters()
{
    if ((func_id < 0) || (func_id >= FUNC_COUNT)) {
        func_id = 0;
    }
    if (a >= b) {
        a = -2;
        b = 2;
    }
    if (d <= c) {
        c = -2;
        d = 2;
    }
    if (thread_count < 0) {
        thread_count = 1;
    }

    eps = abs(eps);
    if (n < 3) {
        n = 3;
    }

    if (m < 3) {
        m = 3;
    }
}

void Graph::normalizeRotationAngle()
{
    zRot -= (int)(zRot / 360.f) * 360.f;
}

double abs_max(double f_val, double f_max_abs)
{
    double f_val_abs = abs(f_val);
    if (f_val_abs > f_max_abs) {
        f_max_abs = f_val_abs;
    }
    return f_max_abs;
}

void Graph::getVertexArray()
{

    double x_center = (a + b) / 2;
    double y_center = (c + d) / 2;

    double s = pow(2., scale);
    double x_size = (b - a) / s;
    double y_size = (d - c) / s;

    double x1 = x_center - x_size / 2;
    double x2 = x_center + x_size / 2;

    double y1 = y_center - y_size / 2;
    double y2 = y_center + y_size / 2;

    try {

        Array buffer(n * m * P * P + n + 2 + m);
        //    fillArray(sampleMin, sampleMax, xdata, nx);
        //    fillArray(sampleMin, sampleMax, ydata, ny);
        //    std::ofstream out("/home/vlatse/Temp/G.txt");

        double mid_point_error = error_factor * 0.1f * f_max_abs;

        const FuncInfo *finfo = &FUNCTIONS[func_id];
        double (*d2fdxdy)(double, double) = finfo->d2fdxdy;
        initG(buffer.data, finfo->f, finfo->dfdx, finfo->dfdy, d2fdxdy, n, m,
              x1, x2, y1, y2, mid_point_error);

        GLfloat dx = (x2 - x1) / (nx - 1);
        GLfloat dy = (y2 - y1) / (ny - 1);

        int i;
        int j;
        f_max_abs = 1e-12;
        f_error_abs_max = 1e-12;
        f_interp_abs_max = 1e-12;

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                double f_val = evaluateFunction(x1, y1, i, j, dx, dy);
                f_max_abs = abs_max(f_val, f_max_abs);

                VertexArray1[i + nx * j][0] = x1 + i * dx;
                VertexArray1[i + nx * j][1] = y1 + j * dy;
                VertexArray1[i + nx * j][2] = f_val;
                VertexArray2[i + nx * j][0] = x1 + i * dx;
                VertexArray2[i + nx * j][1] = y1 + j * dy;
                double xi = x1 + i * dx;
                double yi = y1 + j * dy;
                double fij_i = compute(buffer.data, n, m, xi, yi);
                f_interp_abs_max = abs_max(fij_i, f_interp_abs_max);

                VertexArray2[i + nx * j][2] = fij_i;
                VertexArray3[i + nx * j][0] = x1 + i * dx;
                VertexArray3[i + nx * j][1] = y1 + j * dy;
                double error = f_val - fij_i;
                VertexArray3[i + nx * j][2] = error;
                f_error_abs_max = abs_max(error, f_error_abs_max);
            }
        }

        maxpog = 0;
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                if (fabs(VertexArray3[i + nx * j][2]) > fabs(maxpog)) {
                    maxpog = VertexArray3[i + nx * j][2];
                }
            }
        }

    } catch (std::bad_alloc &) {
        QMessageBox messageBox;
        messageBox.critical(0, "Error", "Not enough memory");
        n = 8;
        m = 8;
        defineVertexShape();
    }
}

void Graph::getColorArray() // определить массив цветов вершин
{
    for (int i = 0; i < DEF_X; i++) {
        for (int j = 0; j < DEF_Y; j++) {
            ColorArray1[i + DEF_X * j][0] = 0.3f;
            ColorArray1[i + DEF_X * j][1] = 0.3f;
            ColorArray1[i + DEF_X * j][2] = 1.0f;
            ColorArray2[i + DEF_X * j][0] = 1.0f;
            ColorArray2[i + DEF_X * j][1] = 1.0f;
            ColorArray2[i + DEF_X * j][2] = 0.0f;
            ColorArray3[i + DEF_X * j][0] = 0.5f;
            ColorArray3[i + DEF_X * j][1] = 0.5f;
            ColorArray3[i + DEF_X * j][2] = 0.5f;
        }
    }
}

void Graph::getIndexArray() // определить массив индексов
{
    for (int j = 0; j < (ny - 1); j++) {
        for (int i = 0; i < (nx - 1); i++) {
            IndexArray[i + (nx - 1) * j][0] = j * nx + i;
            IndexArray[i + (nx - 1) * j][1] = (j + 1) * nx + i;
            IndexArray[i + (nx - 1) * j][2] = (j + 1) * nx + i + 1;
        }
    }
}

void Graph::drawFigure1() // построить фигуру
{
    if (graphf_id == 1) {
        glVertexPointer(3, GL_FLOAT, 0, VertexArray1);
        glColorPointer(3, GL_FLOAT, 0, ColorArray1);
        glDrawElements(GL_TRIANGLES, (nx - 1) * (ny - 1) * 3, GL_UNSIGNED_INT,
                       IndexArray);
    }
}

void Graph::drawFigure2()
{
    if (graph_id == 1) {
        glVertexPointer(3, GL_FLOAT, 0, VertexArray2);
        glColorPointer(3, GL_FLOAT, 0, ColorArray2);
        glDrawElements(GL_TRIANGLES, (nx - 1) * (ny - 1) * 3, GL_UNSIGNED_INT,
                       IndexArray);
    }
}

void Graph::drawFigure3()
{
    if (dif_id == 1) {
        glVertexPointer(3, GL_FLOAT, 0, VertexArray3);
        glColorPointer(3, GL_FLOAT, 0, ColorArray3);
        glDrawElements(GL_TRIANGLES, (nx - 1) * (ny - 1) * 3, GL_UNSIGNED_INT,
                       IndexArray);
    }
}
