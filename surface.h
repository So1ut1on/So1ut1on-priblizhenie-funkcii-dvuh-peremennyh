#ifndef SURFACE_H
#define SURFACE_H
#include "interpolation.h"
#include <QGLWidget>
#include <QtCore>
#include <QtGui>
#include <QtOpenGL>

#define DEF_X 128
#define DEF_Y 128

class Graph : public QGLWidget
{
  private:
    double a;
    double b;
    double c;
    double d;

    int n;
    int m;

    int func_id;
    int graph_id;
    int graphf_id;
    int graph_mode;

    double maxpog;
    int dif_id;
    int ch;
    int chh;

    int nx;
    int ny;
    double f_max_abs;
    double f_interp_abs_max;
    double f_error_abs_max;
    double *f_max_val;
    int error_factor;

    double *mf;
    GLfloat xRot;
    GLfloat yRot;
    GLfloat zRot;
    GLfloat zTra;
    //    GLfloat nSca;
    GLfloat scale_xy;
    int scale;
    GLfloat scale_z;

    double eps;
    int thread_count;

    QPoint ptrMousePosition;

    void scale_plus(float k = 1.1);
    void scale_minus(float k = 1.1);
    bool scale_axes(int k);
    void chang_ch();
    void rotate_up();
    void rotate_down();
    void rotate_left(float angle = 1.0);
    void rotate_right(float angle = 1.0);
    void translate_down();
    void translate_up();
    void defaultScene();
    void chang_graph();
    void chang_dif();
    void changePointsNum(bool increase);
    void chang_matr();
    void chang_graphfunction();
    void chang_ab_m();
    void chang_ab_p();
    void chang_cd_p();
    void chang_cd_m();

    void changeViewType();

    void drawAxis();
    void addError();

    double evaluateFunction(double x1, double y1, int i, int j, double dx,
                            double dy);
    void getVertexArray();
    void getColorArray();
    void getIndexArray();
    void drawFigure1();
    void drawFigure2();
    void drawFigure3();
    void drawText();

    double getDX();
    double getDY();
    double getScaleZ();
    double getScaleXY();

    void defineVertexShape();
    void validateParameters();

    void normalizeRotationAngle();

  protected:
    /*virtual*/ void initializeGL();
    /*virtual*/ void resizeGL(int nWidth, int nHeight);
    /*virtual*/ void paintGL();

    /*virtual*/ void mousePressEvent(QMouseEvent *pe);
    /*virtual*/ void mouseMoveEvent(QMouseEvent *pe);
    /*virtual*/ void mouseReleaseEvent(QMouseEvent *pe);
    /*virtual*/ void wheelEvent(QWheelEvent *pe);
    /*virtual*/ void keyPressEvent(QKeyEvent *pe);

  public:
    Graph(QWidget *parent, int a_, int b_, int c_, int d_, int n_, int m_,
          int k_, double eps_, int p_);
};
#endif
