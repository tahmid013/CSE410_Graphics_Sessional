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


struct point
{
	double x,y,z;
};

class Vector
{
    public:
    double x,y,z;
    Vector(){this->x = 0, this->y = 0,this->z = 0;}
    Vector(struct point p)
    {
        this->x = p.x, this->x = p.y,this->z = p.z;
    }
     Vector(const Vector& v)
     {
         this->x = v.x;
         this->y = v.y;
         this->z = v.z;

     }

        Vector(double x, double y,double z)
        {
            this->x = x, this->y = y,this->z = z;
        }
        Vector operator+(const Vector v);
        Vector operator-(const Vector v);
        Vector operator*(const Vector v);
        Vector operator*(const double d);
        Vector operator^(const Vector v);

    ~Vector() {
        x = y = z = 0.0;
    }

};

    Vector Vector::operator+(const Vector v){
        return Vector(this->x + v.x ,this->y + v.y ,this->z + v.z  );
    }
   Vector Vector::operator-(const Vector v){
        return Vector(this->x - v.x ,this->y - v.y ,this->z - v.z  );
    }
    Vector Vector::operator*(const Vector v){ // dot multiplication
        return Vector(this->x * v.x ,this->y * v.y ,this->z * v.z  );
    }
    Vector Vector::operator*(const double d){ // scaler multiplication
        return Vector(this->x * d ,this->y * d ,this->z * d  );
    }
    Vector Vector::operator^(const Vector v){
        return Vector((this->y*v.z - v.y*this->z ),(this->z*v.x - v.z*this->x ),(this->x*v.y - v.x*this->y ));
    }

Vector pos;
Vector u;
Vector r;
Vector l;




void drawAxes()
{
	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( 100,0,0);
			glVertex3f(-100,0,0);

			glVertex3f(0,-100,0);
			glVertex3f(0, 100,0);

			glVertex3f(0,0, 100);
			glVertex3f(0,0,-100);
		}glEnd();
	}
}


void drawGrid()
{
	int i;
	if(drawgrid==1)
	{
		glColor3f(0.6, 0.6, 0.6);	//grey
		glBegin(GL_LINES);{
			for(i=-8;i<=8;i++){

				if(i==0)
					continue;	//SKIP the MAIN axes

				//lines parallel to Y-axis
				glVertex3f(i*10, -90, 0);
				glVertex3f(i*10,  90, 0);

				//lines parallel to X-axis
				glVertex3f(-90, i*10, 0);
				glVertex3f( 90, i*10, 0);
			}
		}glEnd();
	}
}
void drawOneForth_Cylinder(double radius,double height,int slices,int stacks)
{
   struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=height*sin(((double)i/(double)stacks)*(pi/2));
		r=radius;
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*0.5*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*0.5*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f(0.0f,1.0f,0.0f );
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);

			}glEnd();
		}
	}

}
void drawOneEight_Sphere(double radius,int slices,int stacks)
{
 	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*0.5*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*0.5*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f(1.0f,0.0f,0.0f );
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere

			}glEnd();
		}
	}
}
void drawSquare(double a)
{


    //glColor3f(1.0,0.0,0.0);
	glBegin(GL_QUADS);{
		glVertex3f( a, a,2);
		glVertex3f( a,-a,2);
		glVertex3f(-a,-a,2);
		glVertex3f(-a, a,2);
	}glEnd();
}
void SerVertices(point p){
    glVertex3f( p.x, p.y,p.z);
}
void drawCube(double a)
{
    point A ={-a/2,a/2,a/2} ;
    point B ={a/2,a/2,a/2} ;
    point C ={a/2,a/2,-a/2} ;
    point D ={-a/2,a/2,-a/2} ;
    point E ={a/2,-a/2,-a/2} ;
    point F ={a/2,-a/2,a/2} ;
    point G ={-a/2,-a/2,a/2} ;
    point H ={-a/2,-a/2,-a/2} ;


	glBegin(GL_QUADS);{
	    glColor3f(0.5f, 0.0f, 1.0f); // Color:     Plane:
		SerVertices(A);
		SerVertices(G);
		SerVertices(F);
		SerVertices(B);

	}glEnd();

	glBegin(GL_QUADS);{
	    glColor3f(0.5f, 0.0f, 1.0f); // Color:     Plane:
		SerVertices(B);
		SerVertices(F);
		SerVertices(E);
		SerVertices(C);

	}glEnd();

	glBegin(GL_QUADS);{
	    glColor3f(0.5f, 0.0f, 1.0f); // Color:     Plane:
		SerVertices(A);
		SerVertices(B);
		SerVertices(C);
		SerVertices(D);

	}glEnd();

	glBegin(GL_QUADS);{
	    glColor3f(0.5f, 0.0f, 1.0f); // Color:     Plane:
		SerVertices(G);
		SerVertices(F);
		SerVertices(E);
		SerVertices(H);

	}glEnd();

	glBegin(GL_QUADS);{
	    glColor3f(0.5f, 0.0f, 1.0f); // Color:     Plane:
		SerVertices(H);
		SerVertices(E);
		SerVertices(C);
		SerVertices(D);

	}glEnd();

	glBegin(GL_QUADS);{
	    glColor3f(0.5f, 0.0f, 1.0f); // Color:     Plane:
		SerVertices(G);
		SerVertices(H);
		SerVertices(D);
		SerVertices(A);

	}glEnd();


}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[100];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw segments using generated points
    for(i=0;i<segments;i++)
    {
        glBegin(GL_LINES);
        {
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}

void drawCone(double radius,double height,int segments)
{
    int i;
    double shade;
    struct point points[100];
    //generate points
    for(i=0;i<=segments;i++)
    {
        points[i].x=radius*cos(((double)i/(double)segments)*2*pi);
        points[i].y=radius*sin(((double)i/(double)segments)*2*pi);
    }
    //draw triangles using generated points
    for(i=0;i<segments;i++)
    {
        //create shading effect
        if(i<segments/2)shade=2*(double)i/(double)segments;
        else shade=2*(1.0-(double)i/(double)segments);
        glColor3f(shade,shade,shade);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(0,0,height);
			glVertex3f(points[i].x,points[i].y,0);
			glVertex3f(points[i+1].x,points[i+1].y,0);
        }
        glEnd();
    }
}


void drawSphere(double radius,int slices,int stacks)
{
	struct point points[100][100];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=r*cos(((double)j/(double)slices)*2*pi);
			points[i][j].y=r*sin(((double)j/(double)slices)*2*pi);
			points[i][j].z=h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
                //lower hemisphere
                glVertex3f(points[i][j].x,points[i][j].y,-points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,-points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,-points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,-points[i+1][j].z);
			}glEnd();
		}
	}
}
void drawSS()
{
    glColor3f(1,0,0);
    drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    drawSquare(15);

    glPushMatrix();
    {
        glRotatef(angle,0,0,1);
        glTranslatef(60,0,0);
        glRotatef(2*angle,0,0,1);
        glColor3f(0,0,1);
        drawSquare(10);
    }
    glPopMatrix();

    glRotatef(3*angle,0,0,1);
    glTranslatef(40,0,0);
    glRotatef(4*angle,0,0,1);
    glColor3f(1,1,0);
    drawSquare(5);
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

		case '1':
		    r = r*cos(pi/180) + (u^r)*sin(pi/180);
		    l = l*cos(pi/180) + (u^l)*sin(pi/180);

			break;
        case '2':

		    r = r*cos(-pi/180) + (u^r)*sin(-pi/180);
		    l = l*cos(-pi/180) + (u^l)*sin(-pi/180);

			break;
        case '3':

		    l = l*cos(pi/180) + (r^l)*sin(pi/180);
		    u = u*cos(pi/180) + (r^u)*sin(pi/180);

			break;
        case '4':
            l = l*cos(-pi/180) + (r^l)*sin(-pi/180);
		    u = u*cos(-pi/180) + (r^u)*sin(-pi/180);
			break;
        case '5':
            u = u*cos(pi/180) + (l^u)*sin(pi/180);
		    r = r*cos(pi/180) + (l^r)*sin(pi/180);
			break;
        case '6':
            u = u*cos(-pi/180) + (l^u)*sin(-pi/180);
		    r = r*cos(-pi/180) + (l^r)*sin(-pi/180);
			break;


		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			//cameraHeight -= 3.0;
            pos=pos-l ;


			break;
		case GLUT_KEY_UP:		// up arrow key
			//cameraHeight += 3.0;
            pos=pos+l ;
            break;

		case GLUT_KEY_RIGHT:
			//cameraAngle += 0.03;
			pos=pos+r ;
            break;
			//pos.SetValue(pos.x + r.x , pos.y+r.y ,pos.z+r.z);
			break;
		case GLUT_KEY_LEFT:
			//cameraAngle -= 0.03;
			pos=pos-r ;

			break;

		case GLUT_KEY_PAGE_UP:
		    pos=pos+u ;
			break;
		case GLUT_KEY_PAGE_DOWN:
		    pos=pos-u ;
			break;

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
		    if(rad>0)
                rad--;
			break;
		case GLUT_KEY_END:
		    if(rad<fixed)
                rad++;
			break;

		default:
			break;
	}
}


void mouseListener(int button, int state, int x, int y){	//x, y is the x-y of the screen (2D)
	switch(button){
		case GLUT_LEFT_BUTTON:
			if(state == GLUT_DOWN){		// 2 times?? in ONE click? -- solution is checking DOWN or UP
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



void display(){

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
	//gluLookAt(200*cos(cameraAngle), 200*sin(cameraAngle), cameraHeight,		0,0,0,		0,0,1);
	//gluLookAt(0,0,200,	0,0,0,	0,1,0);
		   // printf("Previous Position: %lf, %lf, %lf\n",  pos.x , pos.y ,pos.z);

	gluLookAt(pos.x,pos.y,pos.z,   // 0,200,200
                pos.x + l.x, pos.y + l.y , pos.z + l.z,  //0 , 0 ,200
                u.x,u.y,u.z);




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

    int slice = 24;
    int stacks = 20;


    double r = rad;
    double h = fixed - r;

    int square_half_a = h;




    glPushMatrix();
    glTranslatef(square_half_a,square_half_a,0);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    //
    glTranslatef(-square_half_a,-square_half_a,0);
    glRotatef(180, 0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();

    glTranslatef(square_half_a,-square_half_a,0);
    glRotatef(270, 0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();

    glTranslatef(-square_half_a,square_half_a,0);
    glRotatef(90, 0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

//                        /*      up down       */
    glPushMatrix();
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90, 0,0,1);
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180, 0,0,1);
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glRotatef(270, 0,0,1);
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();



    /*..........Down.......*/

    glPushMatrix();
    glRotatef(180,0,1,0);

    glPushMatrix();
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glRotatef(90, 0,0,1);
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glRotatef(180, 0,0,1);
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glRotatef(270, 0,0,1);
    glTranslatef(square_half_a,0,square_half_a);
    glRotatef(-45,0,1,0);
    glRotatef(-90,1,0,0);
    glRotatef(-45,0,0,1);
    drawOneForth_Cylinder(r,h,slice, stacks);
    glPopMatrix();
    glPopMatrix();




/*..................sphere..........*/



    glPushMatrix();
    glTranslatef(square_half_a,square_half_a,square_half_a);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-square_half_a,square_half_a,square_half_a);
    glRotatef(90,0,0,1);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-square_half_a,-square_half_a,square_half_a);
    glRotatef(180,0,0,1);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(square_half_a,-square_half_a,square_half_a);
    glRotatef(270,0,0,1);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

     /*..............Down sphere..........*/


    glPushMatrix();
    glRotatef(180,0,1,0);

    glPushMatrix();
    glTranslatef(square_half_a,square_half_a,square_half_a);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-square_half_a,square_half_a,square_half_a);
    glRotatef(90,0,0,1);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(-square_half_a,-square_half_a,square_half_a);
    glRotatef(180,0,0,1);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

    glPushMatrix();
    glTranslatef(square_half_a,-square_half_a,square_half_a);
    glRotatef(270,0,0,1);
    drawOneEight_Sphere(r,slice, stacks);
    glPopMatrix();

    glPopMatrix();




    /*.....white square....*/
    double a = square_half_a+r;


    glPushMatrix();

    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90, 0,1,0);
    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();
    glPushMatrix();

    glColor3f(1.0f,1.0f,1.0f);
    glRotatef(180, 0,1,0);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();
    glPushMatrix();

    glColor3f(1.0f,1.0f,1.0f);
    glRotatef(270, 0,1,0);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();



        /*Rotated White region*/
    glPushMatrix();
    glRotatef(90,1,0,0);

    glPushMatrix();

    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();

    glPushMatrix();
    glRotatef(90, 0,1,0);
    glColor3f(1.0f,1.0f,1.0f);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();
    glPushMatrix();

    glColor3f(1.0f,1.0f,1.0f);
    glRotatef(180, 0,1,0);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();
    glPushMatrix();

    glColor3f(1.0f,1.0f,1.0f);
    glRotatef(270, 0,1,0);
    glBegin(GL_QUADS);{
		glVertex3f( square_half_a, square_half_a,a);
		glVertex3f( square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a,-square_half_a,a);
		glVertex3f(-square_half_a, square_half_a,a);
	}glEnd();
    glPopMatrix();

    glPopMatrix();


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


void animate(){
	angle+=0.01;
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=1;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

	rad = 10;
	fixed =50;


	pos = Vector(100, 100, 0);

    u = Vector(0, 0, 1);
    r = Vector(-0.707, 0.707, 0);
    l = Vector(-0.707, -0.707, 0);





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

int main(int argc, char **argv){
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
