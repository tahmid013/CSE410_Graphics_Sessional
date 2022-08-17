#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))

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
struct point
{
    double x,y,z;
};

class Vector
{
public:
    double x,y,z;
    Vector()
    {
        this->x = 0, this->y = 0,this->z = 0;
    }
    Vector(struct point p)
    {
        this->x = p.x, this->y = p.y,this->z = p.z;
    }
    Vector(const Vector& v)
    {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;

    }
    double abs_value(){
        return sqrt((this->x *this->x) +(this->y *this->y)+(this->z *this->z) );
    }

    Vector(double x, double y,double z)
    {
        this->x = x, this->y = y,this->z = z;
    }

    Vector operator+(const Vector v)
    {
        return Vector(this->x + v.x,this->y + v.y,this->z + v.z  );
    }
    Vector operator-(const Vector v)
    {
        return Vector(this->x - v.x,this->y - v.y,this->z - v.z  );
    }
    double operator*(const Vector v)  // dot multiplication
    {
        return (this->x * v.x + this->y * v.y  + this->z * v.z  );
    }
    Vector operator*(const double d)  // scaler
    {
        return Vector(this->x * d,this->y * d,this->z * d  );
    }
    Vector operator^(const Vector v)
    {
        return Vector((this->y*v.z - v.y*this->z ),(this->z*v.x - v.z*this->x ),(this->x*v.y - v.x*this->y ));
    }


    ~Vector()
    {
        x = y = z = 0.0;
    }

};

Vector dist;
Vector f;

Vector g;
Vector h;



void drawAxes()
{
    if(drawaxes==1)
    {
        glColor3f(1.0, 1.0, 1.0);
        glBegin(GL_LINES);
        {
            glVertex3f( 160,0,0);
            glVertex3f(-160,0,0);

            glVertex3f(0,-160,0);
            glVertex3f(0, 160,0);

            glVertex3f(0,0, 100);
            glVertex3f(0,0,-100);
        }
        glEnd();
    }
}


void drawGrid()
{
    int i;
    if(drawgrid==1)
    {
        glColor3f(0.6, 0.6, 0.6);	//grey
        glBegin(GL_LINES);
        {
            for(i=-15; i<=15; i++)
            {

                if(i==0)
                    continue;	//SKIP the MAIN axes

                //lines parallel to Y-axis
                glVertex3f(i*10, -160, 0);
                glVertex3f(i*10,  160, 0);

                //lines parallel to X-axis
                glVertex3f(-160, i*10, 0);
                glVertex3f( 160, i*10, 0);
            }
        }
        glEnd();
    }
}

void keyboardListener(unsigned char key, int x,int y)
{
    switch(key)
    {

    case 'w':
        theta = theta-speed;
        dist=dist+f*(rad_2*speed*pi/180);

        break;

    case 's':
        theta = theta+speed;
        dist=dist-f*(rad_2*speed*pi/180);

        break;

    case 'a':

        f = f*cos(pi/180) + (g^f)*sin(pi/180);
        h = h*cos(pi/180) + (g^h)*sin(pi/180);
        break;

    case 'd':
        f = f*cos(-pi/180) + (g^f)*sin(-pi/180);
        h = h*cos(-pi/180) + (g^h)*sin(-pi/180);
        break;


    default:
        break;
    }
}


void specialKeyListener(int key, int x,int y)
{
    switch(key)
    {
    case GLUT_KEY_DOWN:		//down arrow key
        cameraHeight -= 3.0;
        ;
        break;
    case GLUT_KEY_UP:		// up arrow key
        cameraHeight += 3.0;

        break;

    case GLUT_KEY_RIGHT:
        cameraAngle += 0.03;

        break;


    case GLUT_KEY_LEFT:
        cameraAngle -= 0.03;


        break;

    case GLUT_KEY_PAGE_UP:

        break;
    case GLUT_KEY_PAGE_DOWN:

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


void mouseListener(int button, int state, int x, int y) 	//x, y is the x-y of the screen (2D)
{
    switch(button)
    {
    case GLUT_LEFT_BUTTON:
        if(state == GLUT_DOWN) 		// 2 times?? in ONE click? -- solution is checking DOWN or UP
        {
            drawaxes=1-drawaxes;
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


void drawCylinder(double radius, double height,int stacks, int slices)
{

    struct point points[100][100];
    int i,j;
    double h,r;


    //generate points
    for(i=0; i<=stacks; i++)
    {
        h=height * ((double)i/(double)stacks);

        for(j=0; j<=slices; j++)
        {
            points[i][j].x=radius*cos(((double)j/(double)slices)*2*pi);
            points[i][j].y=radius*sin(((double)j/(double)slices)*2*pi);
            points[i][j].z=h;
        }
    }
    //draw quads using generated points


    glPushMatrix();
    for(i=0; i<stacks; i++)
    {

        for(j=0; j<slices; j++)
        {
            glBegin(GL_QUADS);
            {
                 glColor3f((double)j/(double)stacks,(double)j/(double)stacks,(double)j/(double)stacks);
                //upper hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
                glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
                glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
                glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere

            }
            glEnd();
        }
    }


    glBegin(GL_QUADS);
    {
        glVertex3f(radius, 0,0);
        glVertex3f(radius, 0,height);
        glVertex3f(-radius, 0,height);
        glVertex3f(-radius, 0,0);

    }
    glEnd();
    glBegin(GL_QUADS);
    {
        glVertex3f(0,radius, 0);
        glVertex3f(0,radius, height);
        glVertex3f(0,-radius, height);
        glVertex3f(0,-radius, 0);

    }
    glEnd();


    glPopMatrix();
}
void drawWheel(double stacks, double slices)
{


    Vector y = Vector(0,1,0);
    double calculated_angle_bwt_f_y=0;
    double dot_product =f.y;

    if(dot_product >1 )
        calculated_angle_bwt_f_y = 0;
    else if(dot_product <-1 )
       calculated_angle_bwt_f_y = 180;

    else{
    if(f.x >= 0)
        calculated_angle_bwt_f_y = -acos((f*y)) *180/pi;
    else if (f.x <= 0)
        calculated_angle_bwt_f_y = acos((f*y)) *180/pi;

    }



    glTranslatef(dist.x,dist.y,0);

    glRotatef(calculated_angle_bwt_f_y,g.x,g.y,g.z);


    glTranslatef(-height_2/2,0,0);
    glTranslatef(0,0,rad_2);
    glRotatef(90, 0,1,0);
    glRotatef(theta, 0,0,1);

    drawCylinder(rad_2, height_2,stacks,slices);

}

void display()
{

    //clear the display
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glClearColor(0,0,0,0);	//color black
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    /********************
    / set-up camera here
    ********************/
    //load the correct matrix -- MODEL-VIEW matrix
    glMatrixMode(GL_MODELVIEW);

    //initialize the matrix
    glLoadIdentity();

    //now give three info
    //1. where is the camera (viewer)?
    //2. where is the camera looking?
    //3. Which direction is the camera's UP direction?

    //gluLookAt(100,100,100,	0,0,0,	0,0,1);
    gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
    //gluLookAt(0,0,200,	0,0,0,	0,1,0);
    // printf("Previous Position: %lf, %lf, %lf\n",  pos.x , pos.y ,pos.z);





    //again select MODEL-VIEW
    glMatrixMode(GL_MODELVIEW);


    /****************************
    / Add your objects from here
    ****************************/
    //add objects

    drawAxes();
    drawGrid();

    //glColor3f(1,0,0);
    //glTranslatef(25,25,25);
    //drawCube(50);

    int slice = 20;
    int stacks = 30;


    //draw_Sphere_To_Square( stacks,  slice);



    drawWheel(slice,stacks);






    //glTranslatef(0,0,15);
    //drawOneEight_Sphere(r,slice, stacks);


    //drawCone(7.0,10, 29);
    //drawSS();

    //drawCircle(30,24);

    //drawCone(20,50,24);

    //drawSphere(30,24,20);




    //ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
    glutSwapBuffers();
}


void animate()
{
    angle+=0.01;
    //codes for any changes in Models, Camera
    glutPostRedisplay();
}

void init()
{
    //codes for initialization
    drawgrid=1;
    drawaxes=1;
    cameraHeight=100.0;
    cameraAngle=1.0;
    angle=0;


    dist = Vector(0,0,0);
    f = Vector(0, 1, 0);
    g = Vector(0, 0, 1);
    h = Vector(1, 0, 0);

    rad_2 = 30;
    height_2 = 10;
    rotate_angle = 0;
    theta = 0;
    speed =3;


    /*
    u = Vector(0, 1, 0);
    //r = {1/sqrt(2), 1/sqrt(2), 0};
    l = Vector(-200, -200, -200);
    */

    //clear the screen
    glClearColor(0,0,0,0);

    /************************
    / set-up projection here
    ************************/
    //load the PROJECTION matrix
    glMatrixMode(GL_PROJECTION);

    //initialize the matrix
    glLoadIdentity();

    //give PERSPECTIVE parameters
    gluPerspective(100,	1,	1,	1000.0);
    //field of view in the Y (vertically)
    //aspect ratio that determines the field of view in the X direction (horizontally)
    //near distance
    //far distance
}

int main(int argc, char **argv)
{
    glutInit(&argc,argv);
    glutInitWindowSize(500, 500);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

    glutCreateWindow("My OpenGL Program");

    init();

    glEnable(GL_DEPTH_TEST);	//enable Depth Testing

    glutDisplayFunc(display);	//display callback function
    glutIdleFunc(animate);		//what you want to do in the idle time (when no drawing is occuring)

    glutKeyboardFunc(keyboardListener);
    glutSpecialFunc(specialKeyListener);
    glutMouseFunc(mouseListener);

    glutMainLoop();		//The main loop of OpenGL

    return 0;
}
