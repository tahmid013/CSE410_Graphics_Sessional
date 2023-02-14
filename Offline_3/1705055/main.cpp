#include <string>
#include <sstream>
#include <windows.h>
#include <GL/glut.h>
#include <iostream>
#include <bits/stdc++.h>
#include <fstream>
#include "1705055_bitmap_image.hpp"
#include "1705055_Header.hpp"
using namespace std;

#define PI (2 * acos(0.0))
#define WIDTH 7
#define RAD(t) (t * PI / 180.0)
#define INF 99999
#define Debug true

int window_width = 500;
int window_height = 500;
double fovY;

extern Vector pos;
Vector u;
Vector r;
Vector l;

int recursionLevel = 0;
extern int level_of_recursion;
int Pixel_Both_Dimension = 0;
int object_count = 0;
int light_count = 0;
int spotlight_count = 0;

extern vector<Object *> objects;
extern vector<PointLight> lights;
extern vector<SpotLight> spotlights;

int img_counter;

void drawAxes()
{

    glColor3f(1.0, 1.0, 1.0);
    glBegin(GL_LINES);
    {
        glVertex3f(100, 0, 0);
        glVertex3f(-100, 0, 0);

        glVertex3f(0, -100, 0);
        glVertex3f(0, 100, 0);

        glVertex3f(0, 0, 100);
        glVertex3f(0, 0, -100);
    }
    glEnd();
}

void capture()
{

    bitmap_image bitmapImage(Pixel_Both_Dimension, Pixel_Both_Dimension);

    //bitmap image initialization to black
    for (int i = 0; i < Pixel_Both_Dimension; i++)
    {
        for (int j = 0; j < Pixel_Both_Dimension; j++)
        {
            bitmapImage.set_pixel(i, j, 0, 0, 0); // color = black
        }
    }
    //calculate plane distance
    double planeDistance = window_height / (2.0 * tan(fovY / 2.0 * PI / 180.0));
    //calculate topleft
    Vector topLeft = pos + l * planeDistance - r * (window_width / 2.0) + u * (window_height / 2.0);

    //calculate du and dv
    double du = (window_width * 1.0) / Pixel_Both_Dimension;
    double dv = (window_height * 1.0) / Pixel_Both_Dimension;

    //Choose middle of the grid cell
    topLeft = topLeft + r * (du * 0.5) - u * (dv * 0.5);

    for (int i = 0; i < Pixel_Both_Dimension; i++)
    {
        for (int j = 0; j < Pixel_Both_Dimension; j++)
        {
            int nearest = INT_MAX;
            double t, tMin = INF;
            //calculate current pixel
            Vector cur_pixel = topLeft + r * (i * du) - u * (j * dv);
            // cast ray from eye to (current pixel - eye) direction
            Ray ray(pos, cur_pixel - pos);
            // for each object in the scene so that min positive
            //value saves the nearest object
            for (int k = 0; k < objects.size(); k++)
            {
                Color dummy__color; // color = black
                t = objects[k]->intersect(ray, dummy__color, 0);

                if (t > 0.0 && t < tMin)
                {
                    tMin = t;
                    nearest = k;
                }
            }
            //update pixel color if nearest object is found
            if (nearest != INT_MAX)
            {
                Color t_color;
                objects[nearest]->intersect(ray, t_color, 1);

                bitmapImage.set_pixel(i, j, (int)round(t_color.R * 255.0), (int)round(t_color.G * 255.0), (int)round(t_color.B * 255.0));
            }
        }
    }

    stringstream cur_img_count;
    cur_img_count << (++img_counter);

    bitmapImage.save_image("D:\\Study\\L-4_T-1_mine\\CSE_410_graphics_sessional\\Ray_tracing\\submission\\1705055\\output_" + cur_img_count.str() + ".bmp");
}

void keyboardListener(unsigned char key, int x, int y)
{
    switch (key)
    {
    case '0':
        capture();
        break;
    case '1':
        r = r * cos(PI / 180) + Vector::CrossMultiplication(u, r) * sin(PI / 180);
        l = l * cos(PI / 180) + Vector::CrossMultiplication(u, l) * sin(PI / 180);
        break;
    case '2':
        r = r * cos(-PI / 180) + Vector::CrossMultiplication(u, r) * sin(-PI / 180);
        l = l * cos(-PI / 180) + Vector::CrossMultiplication(u, l) * sin(-PI / 180);
        break;
    case '3':
        l = l * cos(PI / 180) + Vector::CrossMultiplication(r, l) * sin(PI / 180);
        u = u * cos(PI / 180) + Vector::CrossMultiplication(r, u) * sin(PI / 180);
        break;
    case '4':
        l = l * cos(-PI / 180) + Vector::CrossMultiplication(r, l) * sin(-PI / 180);
        u = u * cos(-PI / 180) + Vector::CrossMultiplication(r, u) * sin(-PI / 180);
        break;
    case '5':
        u = u * cos(PI / 180) + Vector::CrossMultiplication(l, u) * sin(PI / 180);
        r = r * cos(PI / 180) + Vector::CrossMultiplication(l, r) * sin(PI / 180);
        break;
    case '6':
        u = u * cos(-PI / 180) + Vector::CrossMultiplication(l, u) * sin(-PI / 180);
        r = r * cos(-PI / 180) + Vector::CrossMultiplication(l, r) * sin(-PI / 180);
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
    default:
        break;
    }
}

void mouseListener(int button, int state, int x, int y)
{
    switch (button)
    {
    case GLUT_LEFT_BUTTON:
        if (state == GLUT_DOWN) // 2 times?? in ONE click? -- solution is checking DOWN or UP
        {
            // drawaxes=1-drawaxes;
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
    /* clearing the display */
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0, 0, 0, 0); // color = black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);


    // load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    // initialize the matrix
    glLoadIdentity();

    // now give three info
    // 1. where is the camera (viewer)?
    // 2. where is the camera looking?
    // 3. Which direction is the camera's UP direction?

    // gluLookAt(100,100,100,	0,0,0,	0,0,1);
    // gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    // gluLookAt(0,0,200,	0,0,0,	0,1,0);
    //  printf("Previous Position: %lf, %lf, %lf\n",  pos.x , pos.y ,pos.z);
    gluLookAt(pos.x, pos.y, pos.z,                   // 0,200,200
              pos.x + l.x, pos.y + l.y, pos.z + l.z, // 0 , 0 ,200
              u.x, u.y, u.z);

    /* again, selecting MODEL-VIEW matrix */
    glMatrixMode(GL_MODELVIEW);

    /****************************
    / Add your objects from here
    ****************************/
    // add objects

    drawAxes();

    // adding objects
    for (int i = 0; i < objects.size(); i++)
    {
        objects[i]->draw();
    }

    // adding point lights
    for (int i = 0; i < lights.size(); i++)
    {
        lights[i].draw();
    }
    // adding spot lights
    for (int i = 0; i < spotlights.size(); i++)
    {

        spotlights[i].draw();
    }
    /* ADD this line in the end: if you use double buffer (i.e. GL_DOUBLE) */
    glutSwapBuffers();
}

void animate()
{
    // angle+=0.01;
    /* codes for any changes in Models, Camera */
    glutPostRedisplay();
}

void init()
{
    fovY = 80.0;

    u = Vector(0, 0, 1);
    r = Vector(-0.707, 0.707, 0);
    l = Vector(-0.707, -0.707, 0);

    pos = Vector(100, 100, 100);

    img_counter = 0;

    /* clearing the screen */
    glClearColor(0, 0, 0, 0); // color = black

    /* loading the PROJECTION matrix */
    glMatrixMode(GL_PROJECTION);

    /* initializing the matrix */
    glLoadIdentity();

    /* setting the camera perspective by providing necessary parameters */
    gluPerspective(fovY, 1.0, 1.0, 1000.0);
}


void loadData()
{
    ifstream scene;

    scene.open("D:\\Study\\L-4_T-1_mine\\CSE_410_graphics_sessional\\Ray_tracing\\submission\\1705055\\scene_test.txt");

    scene >> recursionLevel >> Pixel_Both_Dimension;
    cout << recursionLevel << Pixel_Both_Dimension << endl;

    scene >> object_count;
    cout <<  object_count << endl;
    Object *object = NULL;

    double shininess;
    for (int i = 0; i < object_count; i++)
    {
        string obj_name;
        scene >> obj_name;
        cout<<obj_name<<endl;
        if (obj_name == "sphere")
        {
            Vector center;
            double radius;

            scene >> center;
            scene >> radius;
            cout<<center<<radius<<endl;

            object = new Sphere(center, radius);
        }
        else if (obj_name == "triangle")
        {
            Vector p1, p2, p3;
            scene >> p1 ;
            scene >> p2 ;
            scene >> p3;
            cout<<p1<<endl;
            cout<<p2<<endl;
            cout<<p3<<endl;
            object = new Triangle(p1, p2, p3);
        }
        else if (obj_name == "general")
        {
            double A, B, C, D, E, F, G, H, I,J;
            Vector center;
            double height, width, depth;
            scene>>A>>B>>C>>D>>E>>F>>G>>H>>I>>J;
            cout<<A<<B<<C<<D<<E<<F<<G<<H<<I<<J<<endl;

            scene>>center;
            cout<<center<<endl;
            scene>>depth;
            cout<<depth<<endl;
            scene>>width;
            cout<<width<<endl;
            scene>>height;
            cout<<height<<endl;

            object = new General_Quadric_Surfaces(A, B, C, D, E, F, G, H, I,J, center, height, width, depth);
        }

        double shininess;
        Color color;
        double ambient, diffuse, specular, rec_ref_cof;
        scene >> color;
        cout<<color<<endl;
        scene >> ambient >> diffuse >> specular >> rec_ref_cof >> shininess;
        cout<<ambient<<diffuse<<specular<<rec_ref_cof<<shininess<<endl;

        object->ambient_coeff = ambient;
        object->diffuse_coeff = diffuse;
        object->specular_coeff = specular;
        object->recursive_coeff = rec_ref_cof;
        object->shine = shininess;

        object->setColor(color);

        objects.push_back(object);

    }

    object = NULL;

    scene >> light_count;
    cout << light_count << endl;
    cout<<"Lights: "<<light_count;
    for (int i = 0; i < light_count; i++)
    {
        Vector position;
        Color color;

        scene >> position;
        scene >> color;
        cout<<position<<color<<endl;
        PointLight light(position, color, 1.5, 10, 10);
        lights.push_back(light);
    }

    scene >> spotlight_count;

    for (int i = 0; i < spotlight_count; i++)
    {
        Vector position_;
        Color color;
        Vector direction;
        double angle;

        scene >> position_;
        scene >> color;
        scene >> direction;
        scene >> angle;

        SpotLight light(position_, color, 1.5, 10, 10, direction, angle);
        spotlights.push_back(light);
    }
    scene.close();

    //with floorwidth 500 and tile width 20
    object = new Floor(500.0, 20.0);
    // default altering color is white
    object->setColor(Color(1.0, 1.0, 1.0));

    //
    object->ambient_coeff= 0.30;
    object->diffuse_coeff= 0.30;
    object->specular_coeff= 0.30;
    object->recursive_coeff= 0.30;
    object->setShininess(25);

    objects.push_back(object);
    object = NULL;


}

int main(int argc, char **argv)
{
    glutInit(&argc, argv);
    glutInitWindowSize(window_width, window_height);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB); // depth, double buffer, rgb color

    glutCreateWindow("Assignment3: Ray Casting & Ray Tracing");

    init(); // initialization

    glEnable(GL_DEPTH_TEST); // enable depth testing

    glutDisplayFunc(display); // display(): callback function
    glutIdleFunc(animate);    // animate(): what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    /* carrying out atexit() function registration */
    /* reference: https://www.geeksforgeeks.org/atexit-function-in-c-c/ */

    /* loading and extracting inputs from input file */
    cout << "Starting...\n";
    loadData();

    glutMainLoop(); // the main loop of OpenGL

    //freeing memory
    for(int i=0;i<objects.size();i++)
    {
        delete objects[i];
    }
    lights.clear();
    spotlights.clear();

    return 0;

}
