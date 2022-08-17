#include <stdio.h>
#include <bits/stdc++.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <windows.h>
#include <GL/glut.h>
#include "bitmap_image.hpp"

#define PI (2 * acos(0.0))
#define WIDTH 7

#define RAD(t) (t * PI / 180.0)
#define INF 99999999
#define Debug true
#define pi (2 * acos(0.0))

using namespace std;
double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;
double rad;
double fixed;

double speed;

double rad_2;
double height_2;
double rotate_angle;
double theta;

// extern vector<Object*> objects;
// extern vector<Light> lights;

int level_rec = 0, pixel_count = 0, obj_count = 0, light_count = 0;
int window_width = 500, window_height = 500;
double fovY = 30;

int recursion_level = 0;
int bitmapImageCount;

struct point
{
    double x, y, z;
};
class Color
{
public:
    double R;
    double G;
    double B;

    Color()
    {
        R = 0;
        G = 0;
        B = 0;
    }
    Color(int r, int g, int b)
    {
        R = r;
        G = g;
        B = b;
    }
    friend istream &operator>>(istream &is, Color &p)
    {
        is >> p.R >> p.G >> p.B;
        return is;
    }
    friend ostream &operator<<(ostream &os, const Color &p)
    {
        os << p.R << " " << p.G << " " << p.B << endl;
        return os;
    }
    ~Color()
    {
    }
};
class Point
{
public:
    double x, y, z;

    Point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    Point(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void normalise()
    {

        double r = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        this->x = this->x / r;
        this->y = this->y / r;
        this->z = this->z / r;
    }
    friend istream &operator>>(istream &is, Point &p)
    {
        is >> p.x >> p.y >> p.z;
        return is;
    }
    friend ostream &operator<<(ostream &os, const Point &p)
    {
        os << p.x << " " << p.y << " " << p.z << endl;
        return os;
    }
    Point operator+(const Point &p2)
        const
    {
        Point res = Point(x + p2.x, y + p2.y, z + p2.z);
        return res;
    }
    Point operator*(const double d)
        const
    {
        Point res = Point(x * d, y * d, z * d);
        return res;
    }
    Point operator-(const Point &p2)
        const
    {
        Point res = Point(x - p2.x, y - p2.y, z - p2.z);
        return res;
    }
    static Point Cross(const Point &u, const Point &v)
    {
        return Point(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
    }
    static double Dot(const Point &u, const Point &v)
    {
        return u.x * v.x + u.y * v.y + u.z * v.z;
    }
    static Point Rodrigauge(const Point &x, const Point &a, const double angle)
    {

        auto cost = cos(RAD(angle));
        auto sint = sin(RAD(angle));

        return x * cost + a * ((1 - cost) * (Point::Dot(a, x))) + (Point::Cross(a, x)) * sint;
    }
    ~Point()
    {
        x = y = z = 0;
    }
};

class Vector
{
public:
    double x, y, z;
    Vector()
    {
        this->x = 0, this->y = 0, this->z = 0;
    }
    Vector(struct point p)
    {
        this->x = p.x, this->y = p.y, this->z = p.z;
    }
    Vector(const Vector &v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
    }
    double abs_value()
    {
        return sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z));
    }
    double Calculatedistance(const Vector &v)
    {
        return sqrt((this->x - v.x) * (this->x - v.x) + (this->y - v.y) * (this->y - v.y) + (this->z - v.z) * (this->z - v.z));
    }


    void normalize()
    {

        double r = sqrt(this->x * this->x + this->y * this->y + this->z * this->z);
        this->x = this->x / r;
        this->y = this->y / r;
        this->z = this->z / r;
    }
    Vector(double x, double y, double z)
    {
        this->x = x, this->y = y, this->z = z;
    }

    Vector operator+(const Vector v)
    {
        return Vector(this->x + v.x, this->y + v.y, this->z + v.z);
    }
    Vector operator-(const Vector v)
    {
        return Vector(this->x - v.x, this->y - v.y, this->z - v.z);
    }
    double operator*(const Vector v) // dot multiplication
    {
        return (this->x * v.x + this->y * v.y + this->z * v.z);
    }
    Vector operator*(const double d) // scaler
    {
        return Vector(this->x * d, this->y * d, this->z * d);
    }
    Vector operator^(const Vector v)
    {
        return Vector((this->y * v.z - v.y * this->z), (this->z * v.x - v.z * this->x), (this->x * v.y - v.x * this->y));
    }
    friend istream &operator>>(istream &is, Vector &v)
    {
        is >> v.x >> v.y >> v.z;
        return is;
    }
    friend ostream &operator<<(ostream &os, const Vector &v)
    {
        os << v.x << " " << v.y << " " << v.z << endl;
        return os;
    }
    ~Vector()
    {
        x = y = z = 0.0;
    }
};

class Light
{
public:
    Vector position;
    Color color;
    double radius;
    int segments;
    int stacks;

    Light()
    {
        radius = 0;
        segments = 0;
        stacks = 0;
    }
    Light(Vector p, Color c, double r, int s, int st)
    {
        position = p;
        color = c;
        radius = r;
        segments = s;
        stacks = st;
    }
    ~Light()
    {
    }
    void draw()
    {
        struct point Vector[stacks + 1][segments + 1];
        int i, j;
        double h, r;
        // generate points

        for (i = 0; i <= stacks; i++)
        {
            h = radius * sin(((double)i / (double)stacks) * (pi / 2));
            r = radius * cos(((double)i / (double)stacks) * (pi / 2));
            for (j = 0; j <= segments; j++)
            {
                // cout<<"i "<<i<<" j "<<j<<endl;
                Vector[i][j].x = r * cos(((double)j / (double)segments) * 2 * pi);
                Vector[i][j].y = r * sin(((double)j / (double)segments) * 2 * pi);
                Vector[i][j].z = h;
            }
        }

        // draw quads using generated points
        for (i = 0; i < stacks; i++)
        {
            glColor3f(color.R, color.G, color.B);
            for (j = 0; j < segments; j++)
            {
                glBegin(GL_QUADS);
                {
                    // upper hemisphere
                    glVertex3f(position.x + Vector[i][j].x, position.y + Vector[i][j].y, position.z + Vector[i][j].z);
                    glVertex3f(position.x + Vector[i][j + 1].x, position.y + Vector[i][j + 1].y, position.z + Vector[i][j + 1].z);
                    glVertex3f(position.x + Vector[i + 1][j + 1].x, position.y + Vector[i + 1][j + 1].y, position.z + Vector[i + 1][j + 1].z);
                    glVertex3f(position.x + Vector[i + 1][j].x, position.y + Vector[i + 1][j].y, position.z + Vector[i + 1][j].z);
                    // lower hemisphere
                    glVertex3f(position.x + Vector[i][j].x, position.y + Vector[i][j].y, position.z - Vector[i][j].z);
                    glVertex3f(position.x + Vector[i][j + 1].x, position.y + Vector[i][j + 1].y, position.z - Vector[i][j + 1].z);
                    glVertex3f(position.x + Vector[i + 1][j + 1].x, position.y + Vector[i + 1][j + 1].y, position.z - Vector[i + 1][j + 1].z);
                    glVertex3f(position.x + Vector[i + 1][j].x, position.y + Vector[i + 1][j].y, position.z - Vector[i + 1][j].z);
                }
                glEnd();
            }
        }
    }
};
class Ray
{
public:
    Vector R_O;
    Vector R_D;
    Ray()
    {
        R_O = Vector(0, 0, 0);
        R_D = Vector(0, 0, 0);
    }
    Ray(Vector R_O, Vector R_D)
    {
        this->R_O = R_O;
        this->R_D = R_D;
    }
    ~Ray()
    {
    }
};

Vector pos;
Vector u;
Vector r;
Vector l;

// vector <PointLight> pointLights;
// vector <SpotLight> spotLights;

class Coefficient
{
    double ambient, diffuse, specular, rec_ref_cof;

public:
    Coefficient()
    {
        ambient = 0, diffuse = 0, specular = 0, rec_ref_cof = 0;
    }
    Coefficient(double ambient, double diffuse, double specular, double rec_ref_cof)
    {
        this->ambient = ambient, this->diffuse = diffuse, this->specular = specular, this->rec_ref_cof = rec_ref_cof;
    }
    ~Coefficient()
    {
    }
};

class Object
{

public:
    Vector reference_point; // should have x, y, z
    double height, width, length;
    Color color;
    Coefficient coEfficients; // ambient, diffuse, specular, reflection coefficients
    int shine;                // exponent term of specular component
    Object()
    {
    }
    virtual void draw() = 0;
    virtual double intersect(Ray ray, Color &t, int level) = 0;
    void setColor(Color c)
    {
        this->color = c;
    }
    void setShine(double shine)
    {
        this->shine = shine;
    }
    void setCoEfficients(Coefficient c)
    {
        this->coEfficients = c;
    }
    ~Object()
    {
    }
};
vector<Object *> objects;
vector<Light> lights;
class Sphere : public Object
{
public:
    Vector center;
    double radius;
    int stacks;
    int segments;

    Sphere()
    {
        stacks = 30;
        segments = 30;
        this->radius = 0;
    }
    Sphere(Vector center, double radius)
    {
        stacks = 30;
        segments = 30;
        this->center = center;
        this->radius = radius;
    }

    void draw()
    {
        // cout << "In sphere draw" << endl;
        // cout << "center " << center << " radius " << radius << endl;
        //  write codes for drawing sphere
        struct point Vector[stacks + 1][segments + 1];
        int i, j;
        double h, r;
        // generate points

        for (i = 0; i <= stacks; i++)
        {
            h = radius * sin(((double)i / (double)stacks) * (pi / 2));
            r = radius * cos(((double)i / (double)stacks) * (pi / 2));
            for (j = 0; j <= segments; j++)
            {
                // cout<<"i "<<i<<" j "<<j<<endl;
                Vector[i][j].x = r * cos(((double)j / (double)segments) * 2 * pi);
                Vector[i][j].y = r * sin(((double)j / (double)segments) * 2 * pi);
                Vector[i][j].z = h;
            }
        }

        // draw quads using generated points
        for (i = 0; i < stacks; i++)
        {
            glColor3f(color.R, color.G, color.B);
            for (j = 0; j < segments; j++)
            {
                glBegin(GL_QUADS);
                {
                    // upper hemisphere
                    glVertex3f(center.x + Vector[i][j].x, center.y + Vector[i][j].y, center.z + Vector[i][j].z);
                    glVertex3f(center.x + Vector[i][j + 1].x, center.y + Vector[i][j + 1].y, center.z + Vector[i][j + 1].z);
                    glVertex3f(center.x + Vector[i + 1][j + 1].x, center.y + Vector[i + 1][j + 1].y, center.z + Vector[i + 1][j + 1].z);
                    glVertex3f(center.x + Vector[i + 1][j].x, center.y + Vector[i + 1][j].y, center.z + Vector[i + 1][j].z);
                    // lower hemisphere
                    glVertex3f(center.x + Vector[i][j].x, center.y + Vector[i][j].y, center.z - Vector[i][j].z);
                    glVertex3f(center.x + Vector[i][j + 1].x, center.y + Vector[i][j + 1].y, center.z - Vector[i][j + 1].z);
                    glVertex3f(center.x + Vector[i + 1][j + 1].x, center.y + Vector[i + 1][j + 1].y, center.z - Vector[i + 1][j + 1].z);
                    glVertex3f(center.x + Vector[i + 1][j].x, center.y + Vector[i + 1][j].y, center.z - Vector[i + 1][j].z);
                }
                glEnd();
            }
        }
    }
    double getTmin(Ray ray)
    {
        double a = 1;
        double b = (ray.R_D * (ray.R_O - center)) * 2;
        double c = ((ray.R_O - center) * (ray.R_O - center)) - radius * radius;
        double disc = b * b - 4 * a * c;
        double t;
        if (disc < 0)
        {
            return -1;
        }
        else if (disc > 0)
        {
            double e = sqrt(disc);
            double t0 = (-b - e) / (2 * a);
            double t1 = (-b + e) / (2 * a);
            if (t0 > 0.001)
            {
                t = t0;
            }
            else if (t1 > 0.001)
            {
                t = t1;
            }
            else
            {
                t = t1;
            }
            return t;
        }
        else
        {
            t = -b / (2 * a);
            return t;
        }
    }
    double intersect(Ray ray, Color &t, int level)
    {
        double tMin = getTmin(ray);


        return tMin;
    }
    ~Sphere()
    {
    }
};

void drawAxes()
{
    if (drawaxes == 1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f(160, 0, 0);
            glVertex3f(-160, 0, 0);

            glVertex3f(0, -160, 0);
            glVertex3f(0, 160, 0);

            glVertex3f(0, 0, 100);
            glVertex3f(0, 0, -100);
        }
        glEnd();
    }
}

void drawGrid()
{
    int i;
    if (drawgrid == 1)
    {
        glColor3f(0.6, 0.6, 0.6); // grey
        glBegin(GL_LINES);
        {
            for (i = -15; i <= 15; i++)
            {

                if (i == 0)
                    continue; // SKIP the MAIN axes

                // lines parallel to Y-axis
                glVertex3f(i * 10, -160, 0);
                glVertex3f(i * 10, 160, 0);

                // lines parallel to X-axis
                glVertex3f(-160, i * 10, 0);
                glVertex3f(160, i * 10, 0);
            }
        }
        glEnd();
    }
}

void capture()
{
    cout << "In capture" << endl;
    bitmap_image bitmapImage(pixel_count, pixel_count);

    for (int i = 0; i < pixel_count; i++)
    {
        for (int j = 0; j < pixel_count; j++)
        {
            Color c = Color(0, 0, 0);
            bitmapImage.set_pixel(i, j, c.R, c.G, c.B); // black
        }
    }
    double planeDistance = (window_height / 2.0) / tan(fovY * (pi / 180) / 2.0);
    Vector topLeft = pos + l * planeDistance - r * (window_width / 2.0) + u * (window_height / 2.0);
    double du = window_width / pixel_count;
    double dv = window_height / pixel_count;

    topLeft = topLeft + r * (0.5 * du) - u * (0.5 * dv);

    int nearest;
    double t, tMin;
    cout << "Done loading init" << endl;
    /*
    for (int i = 0; i < pixel_count; i++)
    {
        for (int j = 0; j < pixel_count; j++)
        {
            Vector cur_pixel = topLeft + r * (i * du) - u * (j * dv);
            Vector rayDirection = cur_pixel - pos;
            rayDirection.normalize();

            Ray ray = Ray(pos, rayDirection);
            tMin = 1000000;
            nearest = -1;
            for (int k = 0; k < objects.size(); k++)
            {
                Color c;
                t = objects[k]->intersect(ray, c, 0);
                if (t > 0 && t < tMin)
                {
                    tMin = t;
                    nearest = k;
                }
            }
            if (nearest != -1)
            {
                Color c = objects[nearest]->color;
                tMin = t = objects[nearest]->intersect(ray, c, 1);
                bitmapImage.set_pixel(i, j, c.R, c.G, c.B);
            }
        }
    }
    */
    for (int column = 0; column < pixel_count; column++)
    {
        for (int row = 0; row < pixel_count; row++)
        {
            /* calculating current pixel and casting ray from camera to (curPixel-camera) direction */
            Vector curPixel = topLeft + r * (column * du) - u * (row * dv);
            Ray ray(pos, curPixel - pos);

            /* finding nearest intersecting object (if available) */
            int nearest = INT_MAX;
            double t, tMin = INF;

            for (int i = 0; i < objects.size(); i++)
            {
                Color color; // color = black
                t = objects[i]->intersect(ray, color, 0);

                if (t > 0.0 && t < tMin)
                {
                    tMin = t;
                    nearest = i;
                }
            }

            /* finding color for current pixel */
            if (nearest != INT_MAX)
            {
                Color color; // color = black
                tMin = objects[nearest]->intersect(ray, color, 1);
                bitmapImage.set_pixel(column, row, (int)round(color.R * 255.0), (int)round(color.G * 255.0), (int)round(color.B * 255.0));
            }
        }
    }

    stringstream currentBitmapImageCount;
    currentBitmapImageCount << (++bitmapImageCount);

    bitmapImage.save_image("D:\\Study\\L-4 T-1_mine\\CSE_410_graphics_sessional\\Ray_tracing\\Code\\output" + currentBitmapImageCount.str() + ".bmp");
    cout << pos << ": bitmap image captured" << endl;
}

void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '0':
        capture();
        break;
    case '1':
        r = r * cos(pi / 180) + (u ^ r) * sin(pi / 180);
        l = l * cos(pi / 180) + (u ^ l) * sin(pi / 180);

        break;
    case '2':

        r = r * cos(-pi / 180) + (u ^ r) * sin(-pi / 180);
        l = l * cos(-pi / 180) + (u ^ l) * sin(-pi / 180);

        break;
    case '3':

        l = l * cos(pi / 180) + (r ^ l) * sin(pi / 180);
        u = u * cos(pi / 180) + (r ^ u) * sin(pi / 180);

        break;
    case '4':
        l = l * cos(-pi / 180) + (r ^ l) * sin(-pi / 180);
        u = u * cos(-pi / 180) + (r ^ u) * sin(-pi / 180);
        break;
    case '5':
        u = u * cos(pi / 180) + (l ^ u) * sin(pi / 180);
        r = r * cos(pi / 180) + (l ^ r) * sin(pi / 180);
        break;
    case '6':
        u = u * cos(-pi / 180) + (l ^ u) * sin(-pi / 180);
        r = r * cos(-pi / 180) + (l ^ r) * sin(-pi / 180);
        break;

    default:
        break;
    }
}

void specialKeyListener(int key, int x, int y)
{
    switch (key)
    {
    case GLUT_KEY_DOWN: // down arrow key
        // cameraHeight -= 3.0;
        pos = pos - l;
        break;
    case GLUT_KEY_UP: // up arrow key
        // cameraHeight += 3.0;
        pos = pos + l;
        break;

    case GLUT_KEY_RIGHT:
        // cameraAngle += 0.03;
        pos = pos + r;
        break;
        // pos.SetValue(pos.x + r.x , pos.y+r.y ,pos.z+r.z);
        break;
    case GLUT_KEY_LEFT:
        // cameraAngle -= 0.03;
        pos = pos - r;

        break;

    case GLUT_KEY_PAGE_UP:
        pos = pos + u;
        break;
    case GLUT_KEY_PAGE_DOWN:
        pos = pos - u;
        break;

    case GLUT_KEY_INSERT:
        break;

    case GLUT_KEY_HOME:

        break;
    case GLUT_KEY_END:

        break;

    default:
        break;
    }
}

void mouseListener(int button, int state, int x, int y) // x, y is the x-y of the screen (2D)
{
    switch (button)
    {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) // 2 times?? in ONE click? -- solution is checking DOWN or UP
        {
            drawaxes = 1 - drawaxes;
        }
        break;

    case GLUT_RIGHT_BUTTON:
        //........
        break;

    case GLUT_MIDDLE_BUTTON:
        //........
        break;

    default:
        break;
    }
}

void display()
{

    // clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0); // color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    // load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    // initialize the matrix
    glLoadIdentity();

    // now give three info
    // 1. where is the camera (viewer)?
    // 2. where is the camera looking?
    // 3. Which direction is the camera's UP direction?

    // gluLookAt(100,100,100,	0,0,0,	0,0,1);
    // gluLookAt(200 * cos(cameraAngle), 200 * sin(cameraAngle), cameraHeight, 0, 0, 0, 0, 0, 1);
    // gluLookAt(0,0,200,	0,0,0,	0,1,0);
    //  printf("Previous Position: %lf, %lf, %lf\n",  pos.x , pos.y ,pos.z);
    gluLookAt(pos.x, pos.y, pos.z,                   // 0,200,200
              pos.x + l.x, pos.y + l.y, pos.z + l.z, // 0 , 0 ,200
              u.x, u.y, u.z);

    // again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);

    /****************************
    / Add your objects from here
    ****************************/
    // add objects

    drawAxes();
    drawGrid();

    // glColor3f(1,0,0);
    // glTranslatef(25,25,25);
    // drawCube(50);

    int slice = 20;
    int stacks = 30;

    // draw_Sphere_To_Square( stacks,  slice);

    // drawWheel(slice,stacks);
    /*drawAxes();
    Point c(0,0,0);
    Sphere sphere(c, 10);
    sphere.draw();
*/

    for (int i = 0; i < objects.size(); i++)
    {
        // cout << "drawing " << i << "\n";
        objects[i]->draw();
    }

    for (int i = 0; i < lights.size(); i++)
    {
        // cout << "drawing " << i << "\n";
        lights[i].draw();
    }
    // glTranslatef(0,0,15);
    // drawOneEight_Sphere(r,slice, stacks);

    // drawCone(7.0,10, 29);
    // drawSS();

    // drawCircle(30,24);

    // drawCone(20,50,24);

    // drawSphere(30,24,20);

    // ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}

void animate()
{
    angle += 0.01;
    // codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    // codes for initialization

    fovY = 80.0;

    drawgrid = 1;
    drawaxes = 1;
    cameraHeight = 100.0;
    cameraAngle = 1.0;
    angle = 0;

    pos = Vector(100, 100, 0);

    u = Vector(0, 0, 1);
    r = Vector(-0.707, 0.707, 0);
    l = Vector(-0.707, -0.707, 0);

    rad_2 = 30;
    height_2 = 10;
    rotate_angle = 0;
    theta = 0;
    speed = 3;

    /*
    u = Vector(0, 1, 0);
    //r = {1/sqrt(2), 1/sqrt(2), 0};
    l = Vector(-200, -200, -200);
    */

    // clear the screen
    glClearColor(0, 0, 0, 0);

    /************************
    / set-up projection here
    ************************/
    // load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    // initialize the matrix
    glLoadIdentity();

    // give PERSPECTIVE parameters
    gluPerspective(100, 1, 1, 1000.0);
    // field of view in the Y (vertically)
    // aspect ratio that determines the field of view in the X direction (horizontally)
    // near distance
    // far distance
}
void loadData()
{

    ifstream scene;
    scene.open("D:\\Study\\L-4 T-1_mine\\CSE_410_graphics_sessional\\Ray_tracing\\Code\\scene_test.txt");
    // cout << "Starting" << endl;
    if (!scene.is_open())
    {
        // cout << "Error opening file" << endl;
        exit(1);
    }

    scene >> level_rec;
    scene >> pixel_count;
    scene >> obj_count;

    Object *object = NULL;
    Coefficient coeff;
    double shininess;
    Color color;
    for (int i = 0; i < obj_count; i++)
    {
        string obj_name;
        scene >> obj_name;
        if (obj_name == "sphere")
        {

            Vector center;
            double radius;

            double ambient, diffuse, specular, rec_ref_cof;

            scene >> center;
            scene >> radius;
            scene >> color;
            scene >> ambient >> diffuse >> specular >> rec_ref_cof >> shininess;
            // cout << "Sphere: " << center << " " << radius << " " << color << " )" << ambient << " " << diffuse << " " << specular << " " << rec_ref_cof << " " << shininess << endl;
            coeff = Coefficient(ambient, diffuse, specular, rec_ref_cof);

            object = new Sphere(center, radius);
        }
        else if (obj_name == "triangle")
        {
            Vector p1, p2, p3;

            Color color;
            double ambient, diffuse, specular, rec_ref_cof;
            double shininess;

            scene >> p1 >> p2 >> p3 >> color >> ambient >> diffuse >> specular >> rec_ref_cof >> shininess;
        }
        else if (obj_name == "general")
        {
        }
        object->setCoEfficients(coeff);
        object->setShine(shininess);
        object->setColor(color);

        objects.push_back(object);
    }

    object = NULL;

    scene >> light_count;
    for (int i = 0; i < light_count; i++)
    {
        Vector pos;
        Color color;
        scene >> pos >> color;
        Light light = Light(pos, color, 1, 10, 5);
        lights.push_back(light);
    }

    // cout << "Finished reading file\n";
}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST); // enable Depth Testing

    glutDisplayFunc(display); // display callback function
    glutIdleFunc(animate);    // what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    loadData();

    glutMainLoop(); // The main loop of OpenGL

    return 0;
}
