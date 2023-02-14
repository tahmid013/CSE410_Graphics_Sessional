
#include <fstream>
#include <vector>
#include <iostream>
#include <bits/stdc++.h>
using namespace std;

#define PI 2 * acos(0.0)
#define INF 99999

// object and light count
extern int light_count;
extern int object_count;

// Vector class
class Vector
{

public:
    double x;
    double y;
    double z;
    Vector()
    {
        x = y = z = 0.0;
    }
    Vector(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }
    void normalize()
    {
        double r = sqrt((this->x * this->x) + (this->y * this->y) + (this->z * this->z));
        this->x = this->x / r;
        this->y = this->y / r;
        this->z = this->z / r;
    }
    double Calculatedistance(const Vector &v)
    {
        return sqrt((this->x - v.x) * (this->x - v.x) + (this->y - v.y) * (this->y - v.y) + (this->z - v.z) * (this->z - v.z));
    }

    Vector operator+(const Vector v)
    {
        return Vector(this->x + v.x, this->y + v.y, this->z + v.z);
    }
    Vector operator-(const Vector v)
    {
        return Vector(this->x - v.x, this->y - v.y, this->z - v.z);
    }

    Vector operator*(const double d) // scaler
    {
        return Vector(this->x * d, this->y * d, this->z * d);
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
    // dot mutiplication
    static double DotMultiplication(Vector v1, Vector v2)
    {
        return ((v1.x * v2.x) + (v1.y * v2.y) + (v1.z * v2.z));
    }
    // cross mutiplication
    static Vector CrossMultiplication(Vector v1, Vector v2)
    {
        return Vector((v1.y * v2.z - v2.y * v1.z), (v1.z * v2.x - v2.z * v1.x), (v1.x * v2.y - v2.x * v1.y));
    }
    ~Vector()
    {
    }
};

// Ray class having origin and direction
class Ray
{

public:
    Vector r_O;
    Vector r_D;
    Ray()
    {
        r_O = Vector(0, 0, 0);
        r_D = Vector(0, 0, 0);
    }
    Ray(Vector r_O, Vector r_D)
    {
        this->r_O = r_O;
        this->r_D = r_D;
        this->r_D.normalize(); // as it is only direction vector so normalize it
    }

    ~Ray()
    {
    }
};

// Color class with operator overloading for easy calculation
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

    Color(double R, double G, double B)
    {
        this->R = R;
        this->G = G;
        this->B = B;
    }

    Color operator+(const Color c)
    {
        return Color(this->R + c.R, this->G + c.G, this->B + c.B);
    }
    Color operator-(const Color c)
    {
        return Color(this->R - c.R, this->G - c.G, this->B - c.B);
    }
    Color operator*(const Color c) // copmonent wise multiplication
    {
        return Color(this->R * c.R, this->G * c.G, this->B * c.B);
    }
    Color operator*(const double d) // scaler multiplication
    {
        return Color(this->R * d, this->G * d, this->B * d);
    }
    void Clip() // cliping values if goes out of boundary
    {
        if (this->R < 0)
        {
            this->R = 0;
        }
        if (this->G < 0)
        {
            this->G = 0;
        }
        if (this->B < 0)
        {
            this->B = 0;
        }
        if (this->R > 1)
        {
            this->R = 1;
        }
        if (this->G > 1)
        {
            this->G = 1;
        }
        if (this->B > 1)
        {
            this->B = 1;
        }
    }

    friend ifstream &operator>>(ifstream &input, Color &color)
    {
        input >> color.R >> color.G >> color.B;
        return input;
    }

    friend ostream &operator<<(ostream &output, Color &color)
    {
        output << color.R << " " << color.G << " " << color.B << endl;
        return output;
    }

    ~Color() {}
};

// PointLight class having position and color
// small sphere like
class PointLight
{

public:
    Vector position;
    Color color;

    double rad;
    int segments;
    int stacks;
    PointLight()
    {
        rad = 0.0;
        segments = 0;
        stacks = 0;
    }

    PointLight(Vector position, Color color, double rad, int segments, int stacks)
    {
        this->position = position;
        this->color = color;

        this->rad = rad;
        this->segments = segments;
        this->stacks = stacks;
    }

    void draw()
    {
        // sphere drawing code

        Vector points[stacks + 1][segments + 1];
        double height, _radius;

        for (int i = 0; i <= stacks; i++)
        {
            height = rad * sin(((double)i / (double)stacks) * (PI / 2));
            _radius = rad * cos(((double)i / (double)stacks) * (PI / 2));

            for (int j = 0; j <= segments; j++)
            {
                points[i][j] = Vector(_radius * cos(((double)j / (double)segments) * 2 * PI), _radius * sin(((double)j / (double)segments) * 2 * PI), height);
            }
        }

        glColor3f(color.R, color.G, color.B); // setting respective color

        for (int i = 0; i < stacks; i++)
        {
            for (int j = 0; j < segments; j++)
            {
                glBegin(GL_QUADS);
                {
                    // upper hemisphere
                    glVertex3f((position + points[i][j]).x, (position + points[i][j]).y, (position + points[i][j]).z);
                    glVertex3f((position + points[i][j + 1]).x, (position + points[i][j + 1]).y, (position + points[i][j + 1]).z);
                    glVertex3f((position + points[i + 1][j + 1]).x, (position + points[i + 1][j + 1]).y, (position + points[i + 1][j + 1]).z);
                    glVertex3f((position + points[i + 1][j]).x, (position + points[i + 1][j]).y, (position + points[i + 1][j]).z);

                    // lower hemisphere
                    glVertex3f((position + points[i][j]).x, (position + points[i][j]).y, (position - points[i][j]).z);
                    glVertex3f((position + points[i][j + 1]).x, (position + points[i][j + 1]).y, (position - points[i][j + 1]).z);
                    glVertex3f((position + points[i + 1][j + 1]).x, (position + points[i + 1][j + 1]).y, (position - points[i + 1][j + 1]).z);
                    glVertex3f((position + points[i + 1][j]).x, (position + points[i + 1][j]).y, (position - points[i + 1][j]).z);
                }
                glEnd();
            }
        }
    }

    ~PointLight() {}
};

// same as PointLight class
// with extra direction vector and cutoff angle

class SpotLight
{

public:
    Vector position;
    Color color;
    Vector direction;
    double cutoffAngle;

    double rad;
    int segments;
    int stacks;
    SpotLight()
    {
        rad = 0.0;
        segments = 0;
        stacks = 0;
        direction = Vector(0, 0, 0);
        cutoffAngle = 0.0;
    }

    SpotLight(Vector position, Color color, double rad, int segments, int stacks, Vector direction, double cutoffAngle)
    {
        this->position = position;
        this->color = color;
        this->direction = direction;
        this->cutoffAngle = cutoffAngle;
        this->rad = rad;
        this->segments = segments;
        this->stacks = stacks;
    }
    // same as spotlight class
    void draw()
    {
        Vector points[stacks + 1][segments + 1];
        double height, _radius;

        for (int i = 0; i <= stacks; i++)
        {
            height = rad * sin(((double)i / (double)stacks) * (PI / 2));
            _radius = rad * cos(((double)i / (double)stacks) * (PI / 2));

            for (int j = 0; j <= segments; j++)
            {
                points[i][j] = Vector(_radius * cos(((double)j / (double)segments) * 2 * PI), _radius * sin(((double)j / (double)segments) * 2 * PI), height);
            }
        }

        glColor3f(color.R, color.G, color.B);

        for (int i = 0; i < stacks; i++)
        {
            for (int j = 0; j < segments; j++)
            {
                glBegin(GL_QUADS);
                {
                    glVertex3f((position + points[i][j]).x, (position + points[i][j]).y, (position + points[i][j]).z);
                    glVertex3f((position + points[i][j + 1]).x, (position + points[i][j + 1]).y, (position + points[i][j + 1]).z);
                    glVertex3f((position + points[i + 1][j + 1]).x, (position + points[i + 1][j + 1]).y, (position + points[i + 1][j + 1]).z);
                    glVertex3f((position + points[i + 1][j]).x, (position + points[i + 1][j]).y, (position + points[i + 1][j]).z);

                    /* lower hemisphere */
                    glVertex3f((position + points[i][j]).x, (position + points[i][j]).y, (position - points[i][j]).z);
                    glVertex3f((position + points[i][j + 1]).x, (position + points[i][j + 1]).y, (position - points[i][j + 1]).z);
                    glVertex3f((position + points[i + 1][j + 1]).x, (position + points[i + 1][j + 1]).y, (position - points[i + 1][j + 1]).z);
                    glVertex3f((position + points[i + 1][j]).x, (position + points[i + 1][j]).y, (position - points[i + 1][j]).z);
                }
                glEnd();
            }
        }
    }

    ~SpotLight() {}
};

// Parent class for all objects in the scene
//  all objects coefficents , color and shinning
class Object
{

public:
    Color color;
    double ambient_coeff, diffuse_coeff, specular_coeff, recursive_coeff;
    int shine;
    Object()
    {
        color = Color(0, 0, 0);
        ambient_coeff = 0;
        diffuse_coeff = 0;
        specular_coeff = 0;
        recursive_coeff = 0;
        shine = 0;
    }

    Color getColor() const
    {
        return color;
    }

    void setColor(Color color)
    {
        this->color = color;
    }

    int getShininy() const
    {
        return shine;
    }

    void setShininess(int shine)
    {
        this->shine = shine;
    }

    void AmbientLightCalc(Color &, Color);
    void RecursiveReflCalc(Color &, Color);

    virtual void draw() = 0;
    virtual double intersect(Ray, Color &, int) = 0;

    ~Object()
    {
        color = Color(0.0, 0.0, 0.0);
        /* destructor */
    }
};

void Object::AmbientLightCalc(Color &color, Color intersectionPointColor)
{
    color = intersectionPointColor * ambient_coeff;
}

void Object::RecursiveReflCalc(Color &color, Color reflected_color)
{
    color = color + reflected_color * recursive_coeff;
}

// global position vector
Vector pos;

extern int recursionLevel;
vector<Object *> objects;
vector<PointLight> lights;
vector<SpotLight> spotlights;

// Sphere class inherits from Object class
class Sphere : public Object
{
public:
    Vector center;
    double radius;
    int stacks;
    int segments;

    Sphere()
    {
        stacks = 25;
        segments = 25;
        this->radius = 0;
    }
    Sphere(Vector center, double radius)
    {
        stacks = 25;
        segments = 25;
        this->center = center;
        this->radius = radius;
    }

    void draw()
    {

        Vector points[stacks + 1][segments + 1];
        int i, j;
        double h, r;

        for (i = 0; i <= stacks; i++)
        {
            h = radius * sin(((double)i / (double)stacks) * (PI / 2));
            r = radius * cos(((double)i / (double)stacks) * (PI / 2));
            for (j = 0; j <= segments; j++)
            {
                points[i][j].x = r * cos(((double)j / (double)segments) * 2 * PI);
                points[i][j].y = r * sin(((double)j / (double)segments) * 2 * PI);
                points[i][j].z = h;
            }
        }

        for (i = 0; i < stacks; i++)
        {
            glColor3f(color.R, color.G, color.B);
            for (j = 0; j < segments; j++)
            {
                glBegin(GL_QUADS);
                {
                    glVertex3f(center.x + points[i][j].x, center.y + points[i][j].y, center.z + points[i][j].z);
                    glVertex3f(center.x + points[i][j + 1].x, center.y + points[i][j + 1].y, center.z + points[i][j + 1].z);
                    glVertex3f(center.x + points[i + 1][j + 1].x, center.y + points[i + 1][j + 1].y, center.z + points[i + 1][j + 1].z);
                    glVertex3f(center.x + points[i + 1][j].x, center.y + points[i + 1][j].y, center.z + points[i + 1][j].z);

                    glVertex3f(center.x + points[i][j].x, center.y + points[i][j].y, center.z - points[i][j].z);
                    glVertex3f(center.x + points[i][j + 1].x, center.y + points[i][j + 1].y, center.z - points[i][j + 1].z);
                    glVertex3f(center.x + points[i + 1][j + 1].x, center.y + points[i + 1][j + 1].y, center.z - points[i + 1][j + 1].z);
                    glVertex3f(center.x + points[i + 1][j].x, center.y + points[i + 1][j].y, center.z - points[i + 1][j].z);
                }
                glEnd();
            }
        }
    }
    // get tMin for intersection
    double getTmin(Ray ray)
    {
        double a = 1;
        double b = Vector::DotMultiplication(ray.r_D, (ray.r_O - center)) * 2;
        double c = Vector::DotMultiplication((ray.r_O - center), (ray.r_O - center)) - radius * radius;
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

    double intersect(Ray ray, Color &color, int level)
    {
        double tMin = getTmin(ray);

        if (level == 0)
        {
            return tMin;
        }

        // store objects color to intersection point color
        Color intersectionPointColor = getColor();

        // get intersection point by ray equation
        Vector intersectingPoint = ray.r_O + ray.r_D * tMin;

        // as it is sphere , normal is outgoing from center of sphere to intersection point
        Vector normal = (intersectingPoint - center);
        normal.normalize();
        if (pos.Calculatedistance(center) <= radius) // if camera is inside sphere , normal is opposite
            normal = normal * -1;

        // Ambient coefficient calculation
        color = (intersectionPointColor * ambient_coeff);

        // for all the point lights do,
        for (int i = 0; i < lights.size(); i++)
        {
            // Casting ray from lights position to intersecting point
            Ray ray1 = Ray(lights[i].position, intersectingPoint - lights[i].position);

            // checking the intersection point is in shadow
            double t, t_min_sh = INF;
            for (int j = 0; j < objects.size(); j++)
            {
                Color temp(0, 0, 0);
                double t = objects[j]->intersect(ray1, temp, 0);
                if (t > 0 && t < t_min_sh)
                    t_min_sh = (t);
            }
            // calculate shadow intersection point by ray eqn
            Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
            double epsilon = 0.0000001;

            // if  in shadow , diffuse and specular coefficient need not to be calculated
            if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
            {
                continue;
            }

            // calculating lambert value (L . N)
            double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
            lambert = max(lambert, 0.0);
            // finding reflected ray by formula ,  2 (L.N)N – L
            Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
            // calculate (R.V)
            double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
            PhongVal = max(PhongVal, 0.0);

            Color AmbientComponent = color;
            //  Ip kd * (L.N)
            Color DiffuseComponent = (lights[i].color * intersectionPointColor) * lambert;
            //  Ip ks * (R.v)^n
            Color SpecularComponent = (lights[i].color * intersectionPointColor) * pow(PhongVal, shine);
            // set color
            color = AmbientComponent + DiffuseComponent + SpecularComponent;
        }

        // almost same as above , except for cutoff angle checking
        for (int i = 0; i < spotlights.size(); i++)
        {

            Ray ray1 = Ray(spotlights[i].position, intersectingPoint - spotlights[i].position);
            // calculate angle between light ray and light direction of spotlight(taken as input)
            double angle = acos(Vector::DotMultiplication(ray1.r_D, spotlights[i].direction)) * 180 / PI;
            // if angle is less than cutoff angle , only then add light
            if (angle > spotlights[i].cutoffAngle)
            {
                continue;
            }

            // shadow
            double t, t_min_sh = INF;
            for (int j = 0; j < objects.size(); j++)
            {
                Color temp(0, 0, 0);
                double t = objects[j]->intersect(ray1, temp, 0);
                if (t > 0 && t < t_min_sh)
                    t_min_sh = (t);
            }
            Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
            double epsilon = 0.0000001;

            if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
            {
                continue;
            }

            double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
            lambert = max(lambert, 0.0);
            Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
            double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
            PhongVal = max(PhongVal, 0.0);

            Color AmbientComponent = color;
            Color DiffuseComponent = (spotlights[i].color * intersectionPointColor) * lambert;
            Color SpecularComponent = (spotlights[i].color * intersectionPointColor) * pow(PhongVal, shine);

            color = AmbientComponent + DiffuseComponent + SpecularComponent;
        }

        if (level >= recursionLevel)
        {
            return tMin;
        }
        // finding reflected ray by formula ,  L- 2 (L.N)N
        Vector reflected_ray_direction = ray.r_D - (normal * (2 * Vector::DotMultiplication(ray.r_D, normal)));
        reflected_ray_direction.normalize();
        // get ray from intersection point to reflected ray direction
        Ray reflected_ray(intersectingPoint + reflected_ray_direction, reflected_ray_direction);

        int near_indx = INT_MAX;
        double t, t_min = INF;
        // for all the objects  check
        // find nearest object index that intersects
        for (int j = 0; j < objects.size(); j++)
        {
            Color temp;
            double t = objects[j]->intersect(reflected_ray, temp, 0);
            if (t > 0 && t < t_min)
            {
                t_min = (t);
                near_indx = j;
            }
        }

        // pass this color to the tracked objects intersect func
        // to get changed color
        Color reflected_color;

        if (near_indx != INT_MAX)
        {
            objects[near_indx]->intersect(reflected_ray, reflected_color, level + 1);
        }

        // calculate recursive reflection coefficient of the reflected ray properties
        RecursiveReflCalc(color, reflected_color);

        // Do clipping  if goes out of bound
        color.Clip();

        return tMin;
    }
    ~Sphere() {}
};

// Triangle class
class Triangle : public Object
{
public:
    Vector p1, p2, p3;
    Triangle()
    {
        p1 = Vector(0, 0, 0);
        p2 = Vector(0, 0, 0);
        p3 = Vector(0, 0, 0);
    }
    Triangle(Vector p1, Vector p2, Vector p3)
    {
        this->p1 = p1;
        this->p2 = p2;
        this->p3 = p3;
    }
    void draw()
    {
        glColor3f(getColor().R, getColor().G, getColor().B);
        glBegin(GL_TRIANGLES);
        glVertex3f(p1.x, p1.y, p1.z);
        glVertex3f(p2.x, p2.y, p2.z);
        glVertex3f(p3.x, p3.y, p3.z);
        glEnd();
    }

    double getTmin(Ray ray)
    {
        // finding barycentric equations variable by cramers rule
        double A_det, beta, alpha, gamma, t_det, tMin;
        double m[3][3];
        m[0][0] = p1.x - p2.x;
        m[0][1] = p1.x - p3.x;
        m[0][2] = ray.r_D.x;
        m[1][0] = p1.y - p2.y;
        m[1][1] = p1.y - p3.y;
        m[1][2] = ray.r_D.y;
        m[2][0] = p1.z - p2.z;
        m[2][1] = p1.z - p3.z;
        m[2][2] = ray.r_D.z;
        A_det = m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);

        double m2[3][3];
        m2[0][0] = p1.x - ray.r_O.x;
        m2[0][1] = p1.x - p3.x;
        m2[0][2] = ray.r_D.x;
        m2[1][0] = p1.y - ray.r_O.y;
        m2[1][1] = p1.y - p3.y;
        m2[1][2] = ray.r_D.y;
        m2[2][0] = p1.z - ray.r_O.z;
        m2[2][1] = p1.z - p3.z;
        m2[2][2] = ray.r_D.z;
        beta = m2[0][0] * (m2[1][1] * m2[2][2] - m2[1][2] * m2[2][1]) - m2[0][1] * (m2[1][0] * m2[2][2] - m2[1][2] * m2[2][0]) + m2[0][2] * (m2[1][0] * m2[2][1] - m2[1][1] * m2[2][0]);
        if (A_det != 0)
        {
            beta = beta / A_det;
        }

        double m3[3][3];
        m3[0][0] = p1.x - p2.x;
        m3[0][1] = p1.x - ray.r_O.x;
        m3[0][2] = ray.r_D.x;
        m3[1][0] = p1.y - p2.y;
        m3[1][1] = p1.y - ray.r_O.y;
        m3[1][2] = ray.r_D.y;
        m3[2][0] = p1.z - p2.z;
        m3[2][1] = p1.z - ray.r_O.z;
        m3[2][2] = ray.r_D.z;
        gamma = m3[0][0] * (m3[1][1] * m3[2][2] - m3[1][2] * m3[2][1]) - m3[0][1] * (m3[1][0] * m3[2][2] - m3[1][2] * m3[2][0]) + m3[0][2] * (m3[1][0] * m3[2][1] - m3[1][1] * m3[2][0]);
        if (A_det != 0)
            gamma = gamma / A_det;

        double m4[3][3];
        m4[0][0] = p1.x - p2.x;
        m4[0][1] = p1.x - p3.x;
        m4[0][2] = p1.x - ray.r_O.x;
        m4[1][0] = p1.y - p2.y;
        m4[1][1] = p1.y - p3.y;
        m4[1][2] = p1.y - ray.r_O.y;
        m4[2][0] = p1.z - p2.z;
        m4[2][1] = p1.z - p3.z;
        m4[2][2] = p1.z - ray.r_O.z;
        t_det = m4[0][0] * (m4[1][1] * m4[2][2] - m4[1][2] * m4[2][1]) - m4[0][1] * (m4[1][0] * m4[2][2] - m4[1][2] * m4[2][0]) + m4[0][2] * (m4[1][0] * m4[2][1] - m4[1][1] * m4[2][0]);
        if (A_det != 0)
            t_det = t_det / A_det;

        // does not intersect
        if (A_det == 0)
        {
            tMin = -1;
        }
        else
        {
            if ((beta + gamma < 1) && (beta > 0) && (gamma > 0))
            {
                tMin = t_det;
            }
            else
            { // does not in boundary
                tMin = -1;
            }
        }
        return tMin;
    }

    double intersect(Ray ray, Color &color, int level)
    {

        double tMin = getTmin(ray);

        if (level == 0)
        {
            return tMin;
        }

        // store objects color to intersection point color
        Color intersectionPointColor = getColor();
        // get intersection point by ray equation
        Vector intersectingPoint = ray.r_O + ray.r_D * tMin;

        // get normal
        Vector normal = Vector::CrossMultiplication((p2 - p1), (p3 - p1));
        normal.normalize();
        // finding right normal direction
        if ((Vector::DotMultiplication(ray.r_D * (-1.0), normal)) <= 0.0)
        {
            normal = normal * (-1.0);
        }

        // Ambient coefficient calculation
        color = (intersectionPointColor * ambient_coeff);

        // for all the point lights do,
        for (int i = 0; i < light_count; i++)
        {
            // Casting ray from lights position to intersecting point
            Ray ray1 = Ray(lights[i].position, intersectingPoint - lights[i].position);
            // checking the intersection point is in shadow
            double t, t_min_sh = INF;
            for (int j = 0; j < object_count; j++)
            {
                Color temp(0, 0, 0);
                double t = objects[j]->intersect(ray1, temp, 0);
                if (t > 0 && t < t_min_sh)
                    t_min_sh = (t);
            }
            // calculate shadow intersection point by ray eqn
            Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
            double epsilon = 0.0000001;
            // if  in shadow , diffuse and specular coefficient need not to be calculated
            if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
            {
                continue;
            }

            // calculating lambert value (L . N)
            double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
            lambert = max(lambert, 0.0);
            // finding reflected ray by formula ,  2 (L.N)N – L
            Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
            // calculate (R.V)
            double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
            PhongVal = max(PhongVal, 0.0);

            Color AmbientComponent = color;
            //  Ip kd * (L.N)
            Color DiffuseComponent = (lights[i].color * intersectionPointColor) * lambert;
            //  Ip ks * (R.v)^n
            Color SpecularComponent = (lights[i].color * intersectionPointColor) * pow(PhongVal, shine);
            // set color
            color = AmbientComponent + DiffuseComponent + SpecularComponent;
        }
        // almost same as above , except for cutoff angle checking
        for (int i = 0; i < spotlights.size(); i++)
        {

            Ray ray1 = Ray(spotlights[i].position, intersectingPoint - spotlights[i].position);
            // calculate angle between light ray and light direction of spotlight(taken as input)
            double angle = acos(Vector::DotMultiplication(ray1.r_D, spotlights[i].direction)) * 180 / PI;
            // if angle is less than cutoff angle , only then add light
            if (angle > spotlights[i].cutoffAngle)
            {
                continue;
            }

            // shadow
            double t, t_min_sh = INF;
            for (int j = 0; j < objects.size(); j++)
            {
                Color temp(0, 0, 0);
                double t = objects[j]->intersect(ray1, temp, 0);
                if (t > 0 && t < t_min_sh)
                    t_min_sh = (t);
            }
            Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
            double epsilon = 0.0000001;

            if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
            {
                continue;
            }

            double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
            lambert = max(lambert, 0.0);
            Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
            double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
            PhongVal = max(PhongVal, 0.0);

            Color AmbientComponent = color;
            Color DiffuseComponent = (spotlights[i].color * intersectionPointColor) * lambert;
            Color SpecularComponent = (spotlights[i].color * intersectionPointColor) * pow(PhongVal, shine);

            color = AmbientComponent + DiffuseComponent + SpecularComponent;
        }

        if (level >= recursionLevel)
        {
            return tMin;
        }
        // finding reflected ray by formula ,  L- 2 (L.N)N
        Vector reflected_ray_direction = ray.r_D - (normal * (2 * Vector::DotMultiplication(ray.r_D, normal)));
        reflected_ray_direction.normalize();
        // get ray from intersection point to reflected ray direction
        Ray reflected_ray(intersectingPoint + reflected_ray_direction, reflected_ray_direction);

        int near_indx = INT_MAX;
        double t, t_min = INF;
        // for all the objects  check
        // find nearest object index that intersects
        for (int j = 0; j < objects.size(); j++)
        {
            Color temp;
            double t = objects[j]->intersect(reflected_ray, temp, 0);
            if (t > 0 && t < t_min)
            {
                t_min = (t);
                near_indx = j;
            }
        }
        // pass this color to the tracked objects intersect func
        // to get changed color
        Color reflected_color;

        if (near_indx != INT_MAX)
        {
            objects[near_indx]->intersect(reflected_ray, reflected_color, level + 1);
        }
        // calculate recursive reflection coefficient of the reflected ray properties
        RecursiveReflCalc(color, reflected_color);
        // Do clipping  if goes out of bound
        color.Clip();

        return tMin;
    }
};
// General quadratic surface class
class General_Quadric_Surfaces : public Object
{
public:
    double A, B, C, D, E, F, G, H, I, J;
    Vector center;
    double height, width, depth;

    General_Quadric_Surfaces(double a, double b, double c, double d, double e, double f, double g, double h, double i, double j, Vector center, double height, double width, double depth)
    {
        A = a;
        B = b;
        C = c;
        D = d;
        E = e;
        F = f;
        G = g;
        H = h;
        I = i;
        J = j;
        this->center = center;
        this->height = height;
        this->width = width;
        this->depth = depth;
    }

    General_Quadric_Surfaces()
    {
        A = 0;
        B = 0;
        C = 0;
        D = 0;
        E = 0;
        F = 0;
        G = 0;
        H = 0;
        I = 0;
        J = 0;
        center = Vector(0, 0, 0);
        height = 0;
        width = 0;
        depth = 0;
    }
    void draw()
    {
    }
    double intersect(Ray ray, Color &color, int level)
    {

        double aQ, bQ, cQ, tMin, tMax;
        double xd = ray.r_D.x;
        double yd = ray.r_D.y;
        double zd = ray.r_D.z;

        double xo = ray.r_O.x;
        double yo = ray.r_O.y;
        double zo = ray.r_O.z;

        // calculate aQ, bQ, cQ

        aQ = A * (xd * xd) + B * (yd * yd) + C * (zd * zd) + D * (xd * yd) + E * (xd * zd) + F * (yd * zd);

        bQ = 2 * (A * xo * xd + B * yo * yd + C * zo * zd) + D * (xo * yd + yo * xd) + E * (xo * zd + zo * xd) + F * (yo * zd + yd * zo) + G * xd + H * yd + I * zd;

        cQ = A * (xo * xo) + B * (yo * yo) + C * (zo * zo) + D * (xo * yo) + E * (xo * zo) + F * (yo * zo) + G * xo + H * yo + I * zo + J;

        double desc = bQ * bQ - 4 * aQ * cQ;

        // if descriminant is less than 0 , no intersection
        if (desc < 0)
        {
            return -1;
        }

        double t0 = (-bQ - sqrt(desc)) / (2 * aQ);
        double t1 = (-bQ + sqrt(desc)) / (2 * aQ);

        Vector intersecting_point_1 = ray.r_O + ray.r_D * t0;
        Vector intersecting_point_2 = ray.r_O + ray.r_D * t1;

        double min_x = center.x;
        double max_x = center.x + depth;
        double min_y = center.y;
        double max_y = center.y + width;
        double min_z = center.z;
        double max_z = center.z + height;

        // tracking by flag
        bool t0_valid = true;
        bool t1_valid = true;

        // set flags to false if the intersecting point is not in the bounding box in respective side
        if (depth > 0)
        {
            if (intersecting_point_1.x < min_x || intersecting_point_1.x > max_x)
            {
                t0_valid = false;
            }
        }
        if (width > 0)
        {
            if (intersecting_point_1.y < min_y || intersecting_point_1.y > max_y)
            {
                t0_valid = false;
            }
        }
        if (height > 0)
        {
            if (intersecting_point_1.z < min_z || intersecting_point_1.z > max_z)
            {
                t0_valid = false;
            }
        }

        if (depth > 0)
        {
            if (intersecting_point_2.x < min_x || intersecting_point_2.x > max_x)
            {
                t1_valid = false;
            }
        }
        if (width > 0)
        {
            if (intersecting_point_2.y < min_y || intersecting_point_2.y > max_y)
            {
                t1_valid = false;
            }
        }
        if (height > 0)
        {
            if (intersecting_point_2.z < min_z || intersecting_point_2.z > max_z)
            {
                t1_valid = false;
            }
        }

        // if both points are valid, return the closer one
        if (t0_valid && t1_valid)
        {
            if (t0 < t1)
            {
                tMin = t0;
                tMax = t1;
            }
            else
            {
                tMin = t1;
                tMax = t0;
            }
        }
        else if (t0_valid)
        {
            tMin = t0;
            tMax = t0;
        }
        else if (t1_valid)
        {
            tMin = t1;
            tMax = t1;
        }
        else
        {
            return -1;
        }

        if (level == 0)
        {
            return tMin;
        }
        Color intersectionPointColor = getColor();

        Vector intersectingPoint = ray.r_O + ray.r_D * tMin;

        // calculate normal component

        double xnormal, ynormal, znormal;
        xnormal = 2 * A * intersectingPoint.x + D * intersectingPoint.y + E * intersectingPoint.z + G;
        ynormal = 2 * B * intersectingPoint.y + D * intersectingPoint.x + F * intersectingPoint.z + H;
        znormal = 2 * C * intersectingPoint.z + E * intersectingPoint.x + F * intersectingPoint.y + I;
        Vector normal(xnormal, ynormal, znormal);

        // adjust normal
        normal = ((Vector::DotMultiplication(ray.r_D * (-1.0), normal)) > 0.0) ? normal : normal * (-1.0);

        // Ambient coefficient calculation
        color = (intersectionPointColor * ambient_coeff);

        // for all the point lights do,
        for (int i = 0; i < light_count; i++)
        {
            // Casting ray from lights position to intersecting point
            Ray ray1 = Ray(lights[i].position, intersectingPoint - lights[i].position);
            // checking the intersection point is in shadow
            double t, t_min_sh = INF;
            for (int j = 0; j < object_count; j++)
            {
                Color temp(0, 0, 0);
                double t = objects[j]->intersect(ray1, temp, 0);
                if (t > 0 && t < t_min_sh)
                    t_min_sh = (t);
            }
            // calculate shadow intersection point by ray eqn
            Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
            double epsilon = 0.0000001;

            // if  in shadow , diffuse and specular coefficient need not to be calculated
            if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
            {
                continue;
            }

            // calculating lambert value (L . N)
            double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
            lambert = max(lambert, 0.0);
            // finding reflected ray by formula ,  2 (L.N)N – L
            Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
            // calculate (R.V)
            double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
            PhongVal = max(PhongVal, 0.0);

            Color AmbientComponent = color;
            //  Ip kd * (L.N)
            Color DiffuseComponent = (lights[i].color * intersectionPointColor) * lambert;
            //  Ip ks * (R.v)^n
            Color SpecularComponent = (lights[i].color * intersectionPointColor) * pow(PhongVal, shine);
            // set color
            color = AmbientComponent + DiffuseComponent + SpecularComponent;
        }
        // almost same as above , except for cutoff angle checking
        for (int i = 0; i < spotlights.size(); i++)
        {

            Ray ray1 = Ray(spotlights[i].position, intersectingPoint - spotlights[i].position);
            // calculate angle between light ray and light direction of spotlight(taken as input)
            double angle = acos(Vector::DotMultiplication(ray1.r_D, spotlights[i].direction)) * 180 / PI;
            // if angle is less than cutoff angle , only then add light
            if (angle > spotlights[i].cutoffAngle)
            {
                continue;
            }

            double t, t_min_sh = INF;
            for (int j = 0; j < objects.size(); j++)
            {
                Color temp(0, 0, 0);
                double t = objects[j]->intersect(ray1, temp, 0);
                if (t > 0 && t < t_min_sh)
                    t_min_sh = (t);
            }
            Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
            double epsilon = 0.0000001;

            if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
            {
                continue;
            }

            // calculating lambert value (L . N)
            double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
            lambert = max(lambert, 0.0);
            // finding reflected ray by formula ,  2 (L.N)N – L
            Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
            // calculate (R.V)
            double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
            PhongVal = max(PhongVal, 0.0);

            Color AmbientComponent = color;
            //  Ip kd * (L.N)
            Color DiffuseComponent = (spotlights[i].color * intersectionPointColor) * lambert;
            //  Ip ks * (R.v)^n
            Color SpecularComponent = (spotlights[i].color * intersectionPointColor) * pow(PhongVal, shine);
            // set color
            color = AmbientComponent + DiffuseComponent + SpecularComponent;
        }

        if (level >= recursionLevel)
        {
            return tMin;
        }
        Vector reflected_ray_direction = ray.r_D - (normal * (2 * Vector::DotMultiplication(ray.r_D, normal)));
        reflected_ray_direction.normalize();
        Ray reflected_ray(intersectingPoint + reflected_ray_direction, reflected_ray_direction);

        int near_indx = INT_MAX;
        double t, t_min = INF;
        for (int j = 0; j < objects.size(); j++)
        {
            Color temp;
            double t = objects[j]->intersect(reflected_ray, temp, 0);
            if (t > 0 && t < t_min)
            {
                t_min = (t);
                near_indx = j;
            }
        }

        Color reflected_color; // color = black

        if (near_indx != INT_MAX)
        {
            objects[near_indx]->intersect(reflected_ray, reflected_color, level + 1);
        }

        RecursiveReflCalc(color, reflected_color);

        color.Clip();

        return tMin;
    }
    ~General_Quadric_Surfaces(){};
};

class Floor : public Object
{

public:
    double floorWidth;
    double tileWidth;

    Floor()
    {
        floorWidth = tileWidth = 0.0;
    }

    Floor(double floorWidth, double tileWidth)
    {
        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;
    }

    void draw();
    double intersect(Ray, Color &, int);
    ~Floor()
    {
    }
};
void drawSquare(Vector leftBottomCorner, double tileWidth)
{
    // draw a tile

    glBegin(GL_QUADS);
    {
        glVertex3f(leftBottomCorner.x, leftBottomCorner.y, leftBottomCorner.z);
        glVertex3f(leftBottomCorner.x + tileWidth, leftBottomCorner.y, leftBottomCorner.z);
        glVertex3f(leftBottomCorner.x + tileWidth, leftBottomCorner.y + tileWidth, leftBottomCorner.z);
        glVertex3f(leftBottomCorner.x, leftBottomCorner.y + tileWidth, leftBottomCorner.z);
    }
    glEnd();
}
void Floor::draw()
{
    int row_number = (int)floorWidth / tileWidth;
    int col_number = (int)floorWidth / row_number;

    for (int i = 0; i < row_number; i++)
    {
        for (int j = 0; j < col_number; j++)
        {
            // draw square of color white for even sum
            if ((i + j) % 2 == 0)
            {
                glColor3f(getColor().R, getColor().G, getColor().B);
            }

            else
            {
                // draw square of color black for odd sum
                glColor3f(0, 0, 0);
            }

            Vector leftBottomCorner(-floorWidth / 2.0 + tileWidth * j, -floorWidth / 2.0 + tileWidth * i, 0.0);

            drawSquare(leftBottomCorner, tileWidth);
        }
    }
}

double Floor::intersect(Ray ray, Color &color, int level)
{
    // normal is in positive Z direction
    Vector normal(0.0, 0.0, 1.0);
    // with angle more than 90 degree, normal is in negative Z direction
    if (Vector::DotMultiplication(pos, normal) <= 0.0)
    {
        normal = normal * (-1.0);
    }

    double tMin = INF;

    if (Vector::DotMultiplication(normal, (ray.r_D)) != 0.0)
    {
        // t = -n.RO / n.RD
        //  D = 0
        tMin = Vector::DotMultiplication(normal, ray.r_O) * (-1.0) / Vector::DotMultiplication(normal, ray.r_D);
    }

    if (tMin > 0.0 && tMin < INF) // if inside range
    {

        Vector intersectionPoint = ray.r_O + (ray.r_D) * tMin;

        // discard out of ranged
        if (!(intersectionPoint.x > -floorWidth / 2.0 && intersectionPoint.x < floorWidth / 2.0))
        {
            tMin = INF;
        }
        if (!(intersectionPoint.y > -floorWidth / 2.0 && intersectionPoint.y < floorWidth / 2.0))
        {
            tMin = INF;
        }
    }

    if (level == 0)
    {
        return tMin;
    }

    Vector intersectingPoint = ray.r_O + ray.r_D * tMin;
    Vector init = Vector(-floorWidth / 2.0, -floorWidth / 2.0, 0.0);
    // vector from center to intersection point
    Vector ref_pos = intersectingPoint - init;
    // set respective color

    // default black
    Color intersectionPointColor = Color(0, 0, 0);
    // if sum positive toggle to white
    if ((int)(floor(ref_pos.x / tileWidth) + floor(ref_pos.y / tileWidth)) % 2 == 0)
    {
        intersectionPointColor = getColor();
    }

    // Ambient coefficient calculation
    color = (intersectionPointColor * ambient_coeff);
    // for all the point lights do,
    for (int i = 0; i < light_count; i++)
    {
        // Casting ray from lights position to intersecting point
        Ray ray1 = Ray(lights[i].position, intersectingPoint - lights[i].position);
        // checking the intersection point is in shadowray
        double t, t_min_sh = INF;
        for (int j = 0; j < object_count; j++)
        {
            Color temp(0, 0, 0);
            double t = objects[j]->intersect(ray1, temp, 0);
            if (t > 0 && t < t_min_sh)
                t_min_sh = (t);
        }

        // calculate shadow intersection point by ray eqn
        Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
        double epsilon = 0.0000001;
        // if  in shadow , diffuse and specular coefficient need not to be calculated
        if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
        {
            continue;
        }

        // calculating lambert value (L . N)
        double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
        lambert = max(lambert, 0.0);
        // finding reflected ray by formula ,  2 (L.N)N – L
        Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
        // calculate (R.V)
        double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
        PhongVal = max(PhongVal, 0.0);

        Color AmbientComponent = color;
        //  Ip kd * (L.N)
        Color DiffuseComponent = (lights[i].color * intersectionPointColor) * lambert;
        //  Ip ks * (R.v)^n
        Color SpecularComponent = (lights[i].color * intersectionPointColor) * pow(PhongVal, shine);
        // set color
        color = AmbientComponent + DiffuseComponent + SpecularComponent;
    }
    for (int i = 0; i < spotlights.size(); i++)
    {

        Ray ray1 = Ray(spotlights[i].position, intersectingPoint - spotlights[i].position);
        // calculate angle between light ray and light direction of spotlight(taken as input)
        double angle = acos(Vector::DotMultiplication(ray1.r_D, spotlights[i].direction)) * 180 / PI;

        // if angle is less than cutoff angle , only then add light
        if (angle > spotlights[i].cutoffAngle)
        {
            continue;
        }

        // shadow
        double t, t_min_sh = INF;
        for (int j = 0; j < objects.size(); j++)
        {
            Color temp(0, 0, 0);
            double t = objects[j]->intersect(ray1, temp, 0);
            if (t > 0 && t < t_min_sh)
                t_min_sh = (t);
        }
        Vector intersectShadowpoint = ray1.r_O + ray1.r_D * t_min_sh;
        double epsilon = 0.0000001;

        if (ray1.r_O.Calculatedistance(intersectingPoint) > ray1.r_O.Calculatedistance(intersectShadowpoint) + epsilon)
        {
            continue;
        }

        double lambert = Vector::DotMultiplication((ray1.r_D * (-1.0)), normal);
        lambert = max(lambert, 0.0);
        Ray reflected_ray = Ray(intersectingPoint, ray1.r_D - (normal * (2 * Vector::DotMultiplication(ray1.r_D, normal))));
        double PhongVal = Vector::DotMultiplication(ray.r_D * (-1.0), reflected_ray.r_D);
        PhongVal = max(PhongVal, 0.0);

        Color AmbientComponent = color;
        Color DiffuseComponent = (spotlights[i].color * intersectionPointColor) * lambert;
        Color SpecularComponent = (spotlights[i].color * intersectionPointColor) * pow(PhongVal, shine);

        color = AmbientComponent + DiffuseComponent + SpecularComponent;
    }

    if (level >= recursionLevel)
    {
        return tMin;
    }
    // finding reflected ray by formula ,  L- 2 (L.N)N
    Vector reflected_ray_direction = ray.r_D - (normal * (2 * Vector::DotMultiplication(ray.r_D, normal)));
    reflected_ray_direction.normalize();
    // get ray from intersection point to reflected ray direction
    Ray reflected_ray(intersectingPoint + reflected_ray_direction, reflected_ray_direction);

    int near_indx = INT_MAX;
    double t, t_min = INF;

    // for all the objects  check
    // find nearest object index that intersects
    for (int j = 0; j < objects.size(); j++)
    {
        Color temp;
        double t = objects[j]->intersect(reflected_ray, temp, 0);
        if (t > 0 && t < t_min)
        {
            t_min = (t);
            near_indx = j;
        }
    }

    // pass this color to the tracked objects intersect func
    // to get changed color
    Color reflected_color; // color = black

    if (near_indx != INT_MAX)
    {
        objects[near_indx]->intersect(reflected_ray, reflected_color, level + 1);
    }

    // calculate recursive reflection coefficient of the reflected ray properties
    RecursiveReflCalc(color, reflected_color);

    // Do clipping  if goes out of bound
    color.Clip();

    return tMin;
}
