#include<iostream>
#include<bits/stdc++.h>
#include<fstream>

#define PI (2 * acos(0.0))
#define WIDTH 7

#define RAD(t) (t * PI / 180.0)
#define INF 99999999
#define Debug true

#include "bitmap_image.hpp"
#include <time.h>
using namespace std;

int TriagleCount  = 0;

class Point
{
public:
    double x,y,z;

    Point()
    {
        this->x = 0;
        this->y = 0;
        this->z = 0;
    }
    Point(double x,double y,double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void normalise()
    {

        double r = sqrt(this->x*this->x + this->y*this->y + this->z*this->z);
        this->x = this->x / r;
        this->y = this->y / r;
        this->z = this->z / r;
    }
    friend istream  &operator>>(istream &is, Point &p)
    {
        is >> p.x >> p.y >> p.z ;
        return is;
    }
    friend ostream& operator<<(ostream& os, const Point& p)
    {
        os << p.x <<" "<< p.y <<" "<< p.z <<endl;
        return os;
    }
    Point operator+(const Point &p2)
    const
    {
        Point res =  Point(x + p2.x, y + p2.y, z + p2.z );
        return res;
    }
    Point operator*(const double d)
    const
    {
        Point res =  Point(x * d, y * d, z * d );
        return res;
    }
    Point operator-(const Point &p2)
    const
    {
        Point res = Point(x - p2.x, y - p2.y, z - p2.z );
        return res;
    }
    static Point Cross(const Point &u, const Point &v)
    {
        return Point(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x);
    }
    static double Dot(const Point &u, const Point &v)
    {
        return u.x*v.x + u.y*v.y + u.z*v.z;
    }
    static Point Rodrigauge(const Point &x, const Point &a, const double angle)
    {

        auto cost = cos(RAD(angle));
        auto sint = sin(RAD(angle));

        return x*cost  + a*((1 - cost) * (Point::Dot(a, x) )) +  (Point::Cross(a, x))*sint;

    }

};

class Color
{
public:
    int R;
    int G;
    int B;

};

class Triangle
{
public:
    Point points[3];
    Color color;

};




class Matrix
{
    vector<vector<double>> a;
public:
    Matrix(int r,int c)
    {
        vector<double> row(c,0);
        for(int i =0 ; i<r; i++)
        {
            a.push_back(row);
        }
    }

    void MakeIdentityMatrix()
    {
        int r = a.size();
        for(int i = 0; i<r; i++)
        {
            a[i][i] = 1;
        }
    }
    static Matrix ColumMatrix(Point &p)
    {

        Matrix c = Matrix(4,1);

        c.a[0][0] = p.x;
        c.a[1][0] = p.y;
        c.a[2][0] = p.z;
        c.a[3][0] = 1;

        return c;

    }
    static Matrix Identity(int r = 4, int c = 4)
    {
        Matrix I(r,c);

        for(int i=0; i<r; i++)
        {
            I.a[i][i] = 1;
        }

        return I;

    }
    static Matrix CreateMatrix(Point& p1, Point& p2, Point& p3 )
    {

        Matrix c = Matrix(4,4);

        c.a[0][0] = p1.x, c.a[0][1] = p1.y, c.a[0][2] = p1.z;
        c.a[1][0] = p2.x, c.a[1][1] = p2.y, c.a[1][2] = p2.z;
        c.a[2][0] = p3.x, c.a[2][1] = p3.y, c.a[2][2] = p3.z;

        c.a[3][3] = 1;

        return c;

    }
    static Matrix TranslateMatrix(double tx,double ty,double tz)
    {
        Matrix z(4,4);

        for(int i = 0; i<4; i++)
        {
            z.a[i][i] = 1;
        }
        z.a[0][3] = tx;
        z.a[1][3] = ty;
        z.a[2][3] = tz;
        return z;
    }
    static Matrix ScaleMatrics(double sx,double sy,double sz)
    {
        Matrix z(4,4);

        z.a[3][3] = 1;

        z.a[0][0] = sx;
        z.a[1][1] = sy;
        z.a[2][2] = sz;

        return z;
    }
    static Matrix ProjectionMatrix(double near, double far, double r_, double t)
    {
        Matrix res = Matrix(4, 4);

        res.a[0][0] = near / r_;
        res.a[1][1] = near / t;

        res.a[2][2] = -(far + near) / (far - near);
        res.a[2][3] = -(2 * far * near) / (far - near);

        res.a[3][2] = -1;

        return res;
    }
    static Matrix ModifyPoint(Matrix &m, Point &p)
    {

        return m*(Matrix::ColumMatrix(p));

    }
    static Matrix RotateMatrix(double angle, double ax,double ay,double az)
    {
        Point a = Point(ax, ay,az);
        a.normalise();

        Point i = Point(1,0,0);
        Point j = Point(0,1,0);
        Point k = Point(0,0,1);
        Point c1 = Point::Rodrigauge(i, a, angle);

        Point c2 = Point::Rodrigauge(j, a, angle);
        Point c3 = Point::Rodrigauge(k, a, angle);



        Matrix res = Matrix::Identity();

        res.a[0][0] = c1.x, res.a[0][1] = c2.x, res.a[0][2] = c3.x;
        res.a[1][0] = c1.y, res.a[1][1] = c2.y, res.a[1][2] = c3.y;
        res.a[2][0] = c1.z, res.a[2][1] = c2.z, res.a[2][2] = c3.z;

        return res;
    }

    static Point ToPoint(const Matrix& m )
    {

        if(m.a.size() == 4 && m.a[0].size() == 1)
        {
            return Point(m.a[0][0] /m.a[3][0], m.a[1][0]/m.a[3][0], m.a[2][0]/m.a[3][0]);
        }
        else
        {
            cout<<"Invalid conversion";
            return Point(0,0,0);
        }
    }
    Matrix operator*(const Matrix &m2)
    const
    {
        Matrix temp =Matrix(a.size(),m2.a[0].size());

        for(int i=0; i<a.size(); i++)
        {
            for(int j=0; j<m2.a[i].size(); j++)
            {
                temp.a[i][j] = 0;
                for(int k =0; k<m2.a.size(); k++)
                {
                    temp.a[i][j] += this->a[i][k]* m2.a[k][j];
                }
            }
        }
        return temp;
    }
    Matrix operator+(const Matrix &m2)
    const
    {
        auto r = 4;
        auto c = 4;

        Matrix temp(4,4);

        for(int i=0; i<4; i++)
        {
            for(int j=0; j<4; j++)
            {
                temp.a[i][j] = this->a[i][j] + m2.a[i][j];
            }
        }
        return temp;
    }
    Matrix operator*(const double d)
    const
    {
        int r = this->a.size();
        int c = this->a[0].size();
        Matrix temp =Matrix(r,c);

        for(int i=0; i<r; i++)
        {
            for(int j=0; j<c; j++)
            {
                temp.a[i][j] = this->a[i][j]* d;
            }
        }
        return temp;
    }
    friend ostream& operator<<(ostream& os, const Matrix& m)
    {

        for(int i=0; i<m.a.size(); i++)
        {
            for(int j = 0; j<m.a[i].size(); j++)
            {
                os << m.a[i][j] << " ";
            }
            os << endl;
        }
        return os;

    }
    void print()
    {
        cout<<"\n--------------------------------------\n";
        for(int i =0 ; i<a.size(); i++)
        {
            for(int j =0 ; j<a[i].size(); j++)
            {
                cout<<a[i][j]<<" ";

            }
            cout<<endl;
        }
        cout<<"\n--------------------------------------\n";
    }
};

bool Check_Intersect(double yp, double y1, double y2)
{
    if(((y1 > yp )&&(yp >y2) )|| ((y1 < yp )&&(yp < y2) ))
        return true;
    else
        return false;

}

double GetIntersectingPoint_X(double yp, Point p1, Point p2)
{

    return p1.x + ((yp- p1.y)*(p1.x - p2.x) / (p1.y - p2.y));

}

double GetIntersectingPoint_Z(double yp, Point p1, Point p2)
{

    return p1.z + ((yp- p1.y)*(p1.z - p2.z) / (p1.y - p2.y));

}
int main()
{
    srand(time(0));

    FILE *fp = freopen("scene.txt", "r", stdin);
    ofstream stage1, stage2, stage3;
    stage1.open("stage1.txt");
    stage2.open("stage2.txt");
    stage3.open("stage3.txt");

    stage1 << setprecision(7) << fixed;
    stage2 << setprecision(7) << fixed;
    stage3 << setprecision(7) << fixed;


    Point eye, look, up;

    double fovY, aspectRatio, near, far;
    cin>>eye>>look>>up ;
    cin>>fovY>>aspectRatio>>near>>far;

    Point l,r,u;
    l =look - eye;


    l.normalise();



    r = Point::Cross(l,up);
    r.normalise();

    u = Point::Cross(r,l);

    Matrix trM = Matrix::TranslateMatrix(-eye.x, -eye.y, -eye.z);

    Point m_l = Point(-l.x, -l.y, -l.z);

    Matrix rot_Matrix = Matrix::CreateMatrix(r,u,m_l);

    Matrix ViewTransformMatrix = rot_Matrix * trM;

    /*------------stage 3-----------------*/

    double fovX = fovY * aspectRatio ;

    double t = near * tan(RAD(fovY /2));
    double r_ = near * tan(RAD(fovX /2));

    Matrix projectionMatrix = Matrix::ProjectionMatrix(near, far, r_, t);

    Matrix m = Matrix(4,4);
    m.MakeIdentityMatrix();


    stack<Matrix> S;
    stack<stack<Matrix>> lastState;

    S.push(m);

    string s;

    while(true)
    {
        cin>>s;

        if(s=="triangle")
        {
            Point p[3];
            cin>>p[0]>>p[1]>>p[2];

            for(int i =0; i<3; i++)
            {
                Point model = Matrix::ToPoint(Matrix::ModifyPoint(S.top(), p[i]));
                //ViewTransformMatrix.print();
                Point view = Matrix::ToPoint(Matrix::ModifyPoint(ViewTransformMatrix, model));


                Point projection = Matrix::ToPoint(Matrix::ModifyPoint(projectionMatrix, view));

                stage1 << model;
                stage2 << view;
                stage3 << projection;
            }
            stage1<<endl;
            stage2 << endl;
            stage3 << endl;
            TriagleCount++;
        }
        else if(s=="translate")
        {
            double tx,ty,tz;
            cin>>tx>>ty>>tz;
            S.push(S.top()*(Matrix::TranslateMatrix(tx,ty,tz)));
        }

        else if(s=="scale")
        {
            double sx,sy,sz;

            cin>>sx>>sy>>sz;

            S.push(S.top()*(Matrix::ScaleMatrics(sx,sy,sz)));
            //(Matrix::ScaleMatrics(sx,sy,sz)).print();
        }

        else if(s=="rotate")
        {
            double angle,ax,ay,az;
            cin>>angle>>ax>>ay>>az;

            S.push(S.top()*(Matrix::RotateMatrix(angle,ax,ay,az)));

        }

        else if(s=="push")
        {
            lastState.push(S);
        }

        else if(s=="pop")
        {
            S = lastState.top();
            lastState.pop();
        }
        else if(s=="end")
        {
            break;
        }
    }

    stage1.close();
    stage2.close();
    stage3.close();
    fclose(fp);

    //************** Stage 4 *****************

    //cout<< "Starting stage 4\n";


    ifstream conf;
    conf.open("config.txt");

    //Reading from stage3.txt

    int screenWidth, screenHeight;
    double left_lim_X, right_lim_X, bot_lim_Y, top_limit_Y, front_lim_Z, rear_lim_Z;
    conf>>screenWidth>>screenHeight;
    conf>>left_lim_X;
    right_lim_X = -left_lim_X;
    conf>>bot_lim_Y;
    top_limit_Y = -bot_lim_Y;
    conf >> front_lim_Z >> rear_lim_Z;

    conf.close();

    // Necessary Initialization

    double dx = (right_lim_X - left_lim_X) /screenWidth ;
    double dy = (top_limit_Y- bot_lim_Y) /screenHeight ;

    double topY = top_limit_Y - dy/2 ;
    double leftX = left_lim_X + dx/2 ;

    double botY = bot_lim_Y + dy/2 ;
    double rightX = right_lim_X - dx/2 ;


    cout<< screenHeight <<" "<<screenWidth <<endl;

    //Allocating ZBuffer
    double ** Zbuffer = new double*[screenHeight];

    for(int i= 0; i <screenHeight ; i++)
    {
        Zbuffer[i] = new double[screenWidth];
    }

    //initializong  ZBuffer
    for(int i = 0; i < screenHeight ; i++ )
    {
        for(int j = 0; j < screenWidth ; j++)
        {
            Zbuffer[i][j] = rear_lim_Z;
        }
    }

    //alocating FrameBuffer

    Color ** FrameBuffer = new Color*[screenHeight];

    for(int i= 0; i <screenHeight ; i++)
    {
        FrameBuffer[i] = new Color[screenWidth];
    }

    //initializong  FrameBuffer
    for(int i = 0; i< screenHeight ; i ++)
    {
        for(int j = 0; j< screenWidth ; j++)
        {
            FrameBuffer[i][j].R = 0;
            FrameBuffer[i][j].G = 0;
            FrameBuffer[i][j].B = 0;

        }
    }



    ifstream stage3_;
    stage3_.open("stage3.txt");

//  Data Structure;
    Triangle triangle[TriagleCount];
    cout<<"tr count : "<<TriagleCount<<endl;

    //getting triangle points from stage 3
    for(int i = 0; i < TriagleCount; i++)
    {

        stage3_>>triangle[i].points[0];
        stage3_>>triangle[i].points[1];
        stage3_>>triangle[i].points[2];

        triangle[i].color.R = rand()%255;
        triangle[i].color.G = rand()%255;
        triangle[i].color.B = rand()%255;

    }
    stage3_.close();


    for(int i = 0; i < TriagleCount; i++)
    {

        //cout<<"\n**********   Triangle "<<i<<" *********  \n";
        //cout<<triangle[i].points[0];
        //cout<<triangle[i].points[1];
        //cout<<triangle[i].points[2];

        //Cliping Vertically
        double max_y = max(triangle[i].points[0].y, max(triangle[i].points[1].y, triangle[i].points[2].y));
        double min_y = min(triangle[i].points[0].y, min(triangle[i].points[1].y, triangle[i].points[2].y));

        double top_scan_line, bottom_scan_line;


        if( max_y  > topY )
            top_scan_line = topY;
        else
            top_scan_line = max_y;

        if( min_y > botY )
            bottom_scan_line = min_y;
        else
            bottom_scan_line = botY;

        //Scane Vertically

        for(double yp = top_scan_line ; yp >= bottom_scan_line ; yp -= dy )
        {

            double x_a=5, x_b=5, z_a,z_b;
            int flag = 0;

            for(int j = 0; j<3; j++)
            {

                //ensures that ys will not intersect with 3 linea ...
                if(Check_Intersect(yp,triangle[i].points[j].y, triangle[i].points[(j+1)%3].y))
                {

                    double tobe_checked_x = GetIntersectingPoint_X(yp, triangle[i].points[j], triangle[i].points[(j+1)%3]);
                    double tobe_checked_z = GetIntersectingPoint_Z(yp, triangle[i].points[j], triangle[i].points[(j+1)%3]);

                    if(flag == 0) //
                    {
                        x_a =tobe_checked_x;
                        z_a = tobe_checked_z;
                        flag = 1;
                    }
                    else
                    {
                        x_b  = tobe_checked_x;
                        z_b = tobe_checked_z;
                    };
                }
            }


            //Cliping Horizontally

            int left_intersect_column = min(round(((x_b - leftX) /dx)),round(((x_a - leftX) /dx) ));
            int right_intersect_column = max(round(((x_b - leftX) /dx)),round(((x_a - leftX) /dx) ));

            double constant = (z_b - z_a ) / (x_b - x_a) ;

            double Z_p ;

            if(x_a > x_b)
                Z_p = z_b;
            else
                Z_p = z_a;

            int row = round( (topY-yp) /dy );

            //Scan column
            for(int c = left_intersect_column; c <= right_intersect_column; c++ )
            {

                Z_p += constant;

                //only update within border and set color

                if(c >0 && c< screenWidth)
                {
                    if(Zbuffer[c][row]>Z_p)
                    {
                        Zbuffer[c][row] = Z_p;
                        FrameBuffer[c][row].R = triangle[i].color.R;
                        FrameBuffer[c][row].G = triangle[i].color.G;
                        FrameBuffer[c][row].B = triangle[i].color.B;
                    }
                }
            }
        }
    }


    // set pixel color

    bitmap_image image(screenWidth, screenHeight);
    for (int i = 0; i < screenHeight; i++)
    {
        for(int j=0; j<screenWidth; j++)
        {
            image.set_pixel(i, j, FrameBuffer[i][j].R,FrameBuffer[i][j].G,FrameBuffer[i][j].B);
        }
    }

    image.save_image("output.bmp");

    //printing Zbuffer
    ofstream z_text;
    z_text.open("z_buffer.txt");
    for(int i=0; i<screenHeight; i++)
    {
        int fl = 0;
        for(int j=0; j<screenWidth; j++)
        {
            if(Zbuffer[i][j]<rear_lim_Z)
            {
                fl = 1;
                z_text << Zbuffer[i][j]<<"\t";
            }
        }
        if(fl == 1)
        {
            z_text <<endl;
        }

    }
    z_text.close();

    //Free Mem
    for(int  i = 0; i <screenHeight ; i++)
    {
        delete[] Zbuffer[i];
        delete[] FrameBuffer[i];
    }
    delete Zbuffer;
    delete FrameBuffer;



}
