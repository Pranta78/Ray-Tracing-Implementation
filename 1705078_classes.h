#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#define eps 1e-20
#define pi (2*acos(0.0))
#define root_2 sqrt(2.0)

#define AMBIENT 0
#define DIFFUSE 1
#define SPECULAR 2
#define REFLECTION 0

class Vector3D
{
public:
    double x, y, z;

    Vector3D()  {}

    Vector3D(double x, double y, double z)
    {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    //Vector addition
	Vector3D operator+(const Vector3D& a) const
    {
        return Vector3D(a.x+x, a.y+y, a.z+z);
    }

    //Vector subtraction
    Vector3D operator-(const Vector3D& a) const
    {
        return Vector3D(x-a.x, y-a.y, z-a.z);
    }

    //Scalar multiplication
    Vector3D operator*(double a) const
    {
        return Vector3D(x*a, y*a, z*a);
    }

    //Scalar division
    Vector3D operator/(double a) const
    {
        return Vector3D(x/a, y/a, z/a);
    }

    //dot product
    double operator%(const Vector3D& a) const
    {
        return (x*a.x + y*a.y + z*a.z);
    }

    //Cross product
    Vector3D operator^(const Vector3D& a) const
    {
        return Vector3D(a.z*y - a.y*z, a.x*z - a.z*x, a.y*x - a.x*y);
    }

    double distance()
    {
        return sqrt(x*x + y*y + z*z);
    }

    void normalize()
    {
        double dist = sqrt(x*x + y*y + z*z);
        x /= dist;
        y /= dist;
        z /= dist;
    }

    friend Vector3D operator*(double a, const Vector3D& v)
    {
        return Vector3D(a*v.x, a*v.y, a*v.z);
    }
};

class Color
{
public:
    double r, g, b;

    Color()  {}

    Color(double x, double y, double z)
    {
        this->r = x;
        this->g = y;
        this->b = z;
    }

	Color operator+(const Color& a) const
    {
        return Color(a.r+r, a.g+g, a.b+b);
    }

    Color operator*(double a) const
    {
        return Color(r*a, g*a, b*a);
    }

    friend Color operator*(double a, const Color& v)
    {
        return Color(a*v.r, a*v.g, a*v.b);
    }

    Color operator^(const Color& a) const
    {
        return Color(r*a.r, g*a.g, b*a.b);
    }
};

class PointLight
{
public:
    Vector3D light_pos;
    //double color[3];
    Color color;

    PointLight()    {}

    PointLight(Vector3D light_pos, double r, double g, double b)
    {
        this->light_pos = light_pos;
//        this->color[0] = r;
//        this->color[1] = g;
//        this->color[2] = b;
        color.r = r;
        color.g = g;
        color.b = b;
    }

    void draw()
    {
        drawSphere(5, 100, 100);
    }

    void drawSphere(double radius, int stacks, int slices)
    {
        //Vector3D points[stacks+1][slices+1];
        Vector3D points[101][101];
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

        glPushMatrix();
        {
            glTranslatef(light_pos.x, light_pos.y, light_pos.z);

            glColor3f(color.r, color.g, color.b);

            //draw quads using generated points
            for(i=0;i<stacks;i++)
            {
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
        }glPopMatrix();
    }
};

class SpotLight
{
public:
    PointLight point_light;
    Vector3D light_direction;
    double cutoff_angle;

    SpotLight(PointLight point_light, Vector3D direction, double cutoff_angle)
    {
        this->point_light = point_light;
        this->light_direction = direction;
        this->cutoff_angle = cutoff_angle;
    }
};

//Vector3D cameraPosition = {100, 100, 0};
//Vector3D cameraPosition = {128.482, 96.9838, 52.5};
//Vector3D cameraPosition = {12.1908, -113.561, 55.1955};
//Vector3D cameraPosition = {63.4918, -92.1987, 39};
//Vector3D cameraPosition = {-37.3664, 107.725, 52.5};
//Vector3D cameraPosition = {151.476, 85.2762, 52.5};
//Vector3D cameraPosition = {-9.43898, -89.0654, 52.5};
//Vector3D cameraPosition = {8.73062, -72.1185, 45.0012};
//Vector3D u = {0, 0, 1};
//Vector3D r = {-1.0/root_2, 1.0/root_2, 0};
//Vector3D l = {-1.0/root_2, -1.0/root_2, 0};

//Vector3D cameraPosition = {12.1908, -113.561, 44.6955};
//Vector3D u = {0, 0, 1};
//Vector3D r = {0.999041, -0.0437915, 0};
//Vector3D l = {0.0437915, 0.999041, 0};

//Vector3D cameraPosition = {0, 0, -100};
//Vector3D u = {0, 1, 0};
//Vector3D r = {1, 0, 0};
//Vector3D l = {0, 0, -1};

Vector3D cameraPosition = {21.7096, 67.5, -39.2561};
Vector3D u = {0.0243135, 0.994905, 0.0978404};
Vector3D r = {-0.942009, -0.00996671, 0.335439};
Vector3D l = {0.334705, -0.100322, 0.936967};

vector <PointLight> pointLights;
vector <SpotLight> spotLights;

class Ray
{
public:
    Vector3D start;
    Vector3D dir;

    Ray(Vector3D start, Vector3D dir)
    {
        this->start = start;

        //normalize
        double dist = sqrt(dir.x*dir.x + dir.y*dir.y + dir.z*dir.z);

        this->dir = dir / dist;
    }
};

class Matrix3x3
{
public:
    double matrix[3][3];

    Matrix3x3()     {}

    void copyMatrix(double matrix[3][3])
    {
        for(int i=0; i<3; i++)
            for(int j=0; j<3; j++)
                this->matrix[i][j] = matrix[i][j];
    }

    double determinant()
    {
        return matrix[0][0]*(matrix[1][1]*matrix[2][2] - matrix[2][1]*matrix[1][2])
             - matrix[0][1]*(matrix[1][0]*matrix[2][2] - matrix[1][2]*matrix[2][0])
             + matrix[0][2]*(matrix[1][0]*matrix[2][1] - matrix[2][0]*matrix[1][1]);
    }
};

class Object
{
public:
    Vector3D reference_point;
    double height, width, length;
    //double color[3];
    Color color;
    Color pixelColor;
    double coefficients[4]; //ambient, diffuse, specular, reflection coefficient
    int shine;  //exponent term of specular component

    Object()    {}

    Object(Vector3D reference_point)
    {
        this->reference_point = reference_point;
    }

    Object(Vector3D reference_point, double length, double width, double height)
    {
        this->reference_point = reference_point;
        this->length = length;
        this->width = width;
        this->height = height;
    }

    virtual void draw() {}
    virtual double intersect(Ray *r, double *color, int level)
    {
        return -1.0;
    }

    void setColor(double red, double green, double blue)
    {
//        color[0] = red;
//        color[1] = green;
//        color[2] = blue;

        color.r = red;
        color.g = green;
        color.b = blue;
    }

    void setShine(int shine)
    {
        this->shine = shine;
    }

    void setCoEfficients(double ambient, double diffuse, double specular, double reflection_coefficient)
    {
        coefficients[0] = ambient;
        coefficients[1] = diffuse;
        coefficients[2] = specular;
        coefficients[3] = reflection_coefficient;
    }

    void print()
    {
        cout << "Object Color: r = " << color.r << ", g = " << color.g << ", b = " << color.b << "\n";
    }
};

vector <Object*> objects;

class Sphere : public Object
{
    Vector3D center;
    double radius;

public:
    Sphere(Vector3D center, double radius)
    {
        reference_point = center;
        this->center = reference_point;
        this->radius = radius;
        length = this->radius;
    }

    void draw()
    {
        drawSphere(radius, 100, 100);
    }

    void drawSphere(double radius, int stacks, int slices)
    {
        //Vector3D points[stacks+1][slices+1];
        Vector3D points[101][101];
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

        glPushMatrix();
        {
            //translate to sphere's center
            glTranslatef(center.x, center.y, center.z);

            //glColor3f(color[0], color[1], color[2]);
            glColor3f(color.r, color.g, color.b);

            //draw quads using generated points
            for(i=0;i<stacks;i++)
            {
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
        }glPopMatrix();
    }

    double intersect(Ray *r, double *color, int level)
    {
        //quadratic equation: a*t^2 + b*t + c = 0
        //Ray equation: start + t * dir
        //to accommodate for spheres not centered at origin, R0 is changed
        //R0 = R0 - Center
        Vector3D R0 = r->start - center;
        Vector3D Rd = r->dir;
        //a = 1, b = 2*Rd*R0, c = R0*R0 - r*r
        double a = 1;
        double b = 2 * Rd % R0;
        double c = R0 % R0 - radius * radius;
        double d_square = b*b - 4*a*c;

        double t_min = -1.0;

        //if d^2 < 0, ray does not intersect the sphere
        if(d_square >= 0.0)
        {
            double d = sqrt(d_square);

            double t1 = (-b + d) / (2.0 * a);
            double t2 = (-b - d) / (2.0 * a);

            //take closest positive
            if(t1 >= 0 && t2 >= 0)
                //return min(t1, t2);
                t_min = min(t1, t2);

            else if(t2 < 0.0)
                //return t1;
                t_min = t1;

            else if(t1 < 0.0)
                //return t2;
                t_min = t2;

//            double d = sqrt(d_square);
//            double t = (-b - d) / (2.0 * a);
//
//            if(t > 0.0)
//                t_min = t;
//            else
//                t_min = (-b + d) / (2.0 * a);
        }

        color[0] = min(max(this->color.r, 0.0), 1.0);
        color[1] = min(max(this->color.g, 0.0), 1.0);
        color[2] = min(max(this->color.b, 0.0), 1.0);

        //return -1.0;
        if(level == 0)
            return t_min;

        Vector3D intersectPoint = r->start + r->dir * t_min;
        Color intersectPointColor = this->color;

        //ambient component
        this->pixelColor = intersectPointColor * coefficients[AMBIENT];
        color[0] = min(max(this->pixelColor.r, 0.0), 1.0);
        color[1] = min(max(this->pixelColor.g, 0.0), 1.0);
        color[2] = min(max(this->pixelColor.b, 0.0), 1.0);

        Vector3D normal  = intersectPoint - center;
        normal.normalize();

        for(auto pl : pointLights)
        {
            //direction of the incident ray on the surface
            Vector3D rayl = intersectPoint - pl.light_pos;
            rayl.normalize();

            //find if rayl is obscured by any object
            Ray *ray = new Ray(pl.light_pos, rayl);
            double *dummyColor = new double[3];
            double tMin = 1e50;

            double curT = this->intersect(ray, dummyColor, 0);
            //double curT = 1.0;

            Vector3D rayIntersectPoint = ray->start + ray->dir * curT;

            //intersection point is not reachable by the ray
            if(abs(rayIntersectPoint.x - intersectPoint.x) > 1e-12 || abs(rayIntersectPoint.y - intersectPoint.y) > 1e-12 || abs(rayIntersectPoint.z - intersectPoint.z) > 1e-12)
            {
                //return t_min;
                continue;
            }

            //cout << "Cur t = " << curT << "\n";

            for(auto object : objects)
            {
                double iterT = object->intersect(ray, dummyColor, 0);

                if(iterT >= 0)
                    if(tMin > iterT)
                    {
                        tMin = iterT;
                    }
            }

            //another object is closer to the ray i.e. current object obscured
            if(tMin < curT) //|| curT < 0
            {
                //return t_min;
                continue;
            }

            //diffuse component
            //double lambertValue = (normal % rayl) / (normal.distance() * rayl.distance());
            Vector3D rayInv = pl.light_pos - intersectPoint;
            rayInv.normalize();

            double lambertValue = (normal % rayInv);
            lambertValue = max(lambertValue, 0.0);

            this->pixelColor = this->pixelColor + (((coefficients[DIFFUSE] * lambertValue) * pl.color) ^ intersectPointColor);

            //specular component
            //Vector3D rayv = intersectPoint - cameraPosition;
//            Vector3D rayv = cameraPosition - intersectPoint;
//            rayv.normalize();
//
//            //Vector3D rayr = 2.0 * (rayl % normal) * normal - rayl;
//            //Vector3D rayr = rayInv - ((2.0 * (rayInv % normal)) * normal);
//            //Vector3D rayr = rayl - ((2.0 * (rayl % normal)) * normal);
//            Vector3D rayr = ((2.0 * (rayl % normal)) * normal) - rayl;
//            //Vector3D rayr = ((2.0 * (rayInv % normal)) * normal) - rayInv;
//            rayr.normalize();
//
//            //double phongValue = (rayv % rayr) / (rayv.distance() * rayr.distance());
//            //double phongValue = (rayv % rayr);
//            double phongValue = (rayInv % rayr);

            //cout << "Camera Position: x = " << cameraPosition.x << ", y = " << cameraPosition.y << ", z = " << cameraPosition.z << "\n";
            Vector3D rayv = cameraPosition - intersectPoint;
            //Vector3D rayv = intersectPoint - cameraPosition;
            rayv.normalize();

            //Vector3D rayr = ((2.0 * (rayInv % normal)) * normal) - rayInv;
            //Vector3D rayr = rayl - ((2.0 * (rayl % normal)) * normal);
            Vector3D rayr = rayl - (normal * ((rayl % normal) * 2.0));
            rayr.normalize();

            double phongValue = (rayr % rayv);

            phongValue = max(pow(phongValue, this->shine), 0.0);

            //this->pixelColor = this->pixelColor + pl.color * coefficients[SPECULAR] * phongValue ^ intersectPointColor;
            //this->pixelColor = this->pixelColor + (((coefficients[SPECULAR] * phongValue) * pl.color) ^ intersectPointColor);
            this->pixelColor = this->pixelColor + (((coefficients[SPECULAR] * phongValue) * pl.color));

            color[0] = min(max(this->pixelColor.r, 0.0), 1.0);
            color[1] = min(max(this->pixelColor.g, 0.0), 1.0);
            color[2] = min(max(this->pixelColor.b, 0.0), 1.0);
        }

        for(auto spl : spotLights)
        {
            Vector3D rayl = intersectPoint - spl.point_light.light_pos;

            //check it the angle between rayl and light_direction is smaller than cutoff angle
            double theta = (180.0 / pi) * cos((rayl % spl.light_direction) / (rayl.distance() * spl.light_direction.distance()));

            if(theta > spl.cutoff_angle)
            {
                return t_min;
            }

            //find if rayl is obscured by any object
            Ray *ray = new Ray(spl.point_light.light_pos, rayl);
            double *dummyColor = new double[3];
            double tMin = 1e50;

            double curT = this->intersect(ray, dummyColor, 0);

            for(auto object : objects)
            {
                double iterT = object->intersect(ray, dummyColor, 0);

                if(iterT >= 0)
                    if(tMin > iterT)
                    {
                        tMin = min(tMin, iterT);
                    }
            }

            //another object is closer to the ray i.e. current object obscured
            if(tMin < curT || curT < 0)
            {
                //only ambient components need to be calculated
                return t_min;
            }

            //diffuse component
            double lambertValue = (normal % rayl) / (normal.distance() * rayl.distance());
            lambertValue = max(lambertValue, 0.0);

            this->pixelColor = this->pixelColor + spl.point_light.color * coefficients[DIFFUSE] * lambertValue ^ intersectPointColor;

            //specular component
            Vector3D rayv = intersectPoint - cameraPosition;
            Vector3D rayr = 2.0 * (rayl % normal) * normal - rayl;

            double phongValue = (rayv % rayr) / (rayv.distance() * rayr.distance());
            phongValue = max(pow(phongValue,this->shine*1.0), 0.0);

            this->pixelColor = this->pixelColor + spl.point_light.color * coefficients[SPECULAR] * phongValue ^ intersectPointColor;

            color[0] = this->pixelColor.r;
            color[1] = this->pixelColor.g;
            color[2] = this->pixelColor.b;
        }

        return t_min;
    }
};

class Triangle : public Object
{
    Vector3D point1, point2, point3;

public:
    Triangle(Vector3D point1, Vector3D point2, Vector3D point3)
    {
        this->point1 = point1;
        this->point2 = point2;
        this->point3 = point3;
    }

    void draw()
    {
        //glColor3f(color[0], color[1], color[2]);
        glColor3f(color.r, color.g, color.b);

        glBegin(GL_TRIANGLES);
        {
            glVertex3f(point1.x, point1.y, point1.z);
			glVertex3f(point2.x, point2.y, point2.z);
			glVertex3f(point3.x, point3.y, point3.z);
        }
        glEnd();
    }

    double intersect(Ray *r, double *color, int level)
    {
        color[0] = this->color.r;
        color[1] = this->color.g;
        color[2] = this->color.b;

        double ax = point1.x;
        double ay = point1.y;
        double az = point1.z;

        double bx = point2.x;
        double by = point2.y;
        double bz = point2.z;

        double cx = point3.x;
        double cy = point3.y;
        double cz = point3.z;

        Vector3D R0 = r->start;
        Vector3D Rd = r->dir;

        double R0x = R0.x;
        double R0y = R0.y;
        double R0z = R0.z;

        double Rdx = Rd.x;
        double Rdy = Rd.y;
        double Rdz = Rd.z;

        Matrix3x3 A;
        double matrix1[3][3] = {{ax-bx, ax-cx, Rdx},
                                {ay-by, ay-cy, Rdy},
                                {az-bz, az-cz, Rdz}};
        A.copyMatrix(matrix1);
        double detA = A.determinant();

        Matrix3x3 t_coefficient;
        double matrix2[3][3] = {{ax-bx, ax-cx, ax-R0x},
                                {ay-by, ay-cy, ay-R0y},
                                {az-bz, az-cz, az-R0z}};
        t_coefficient.copyMatrix(matrix2);
        double detT = t_coefficient.determinant() / detA;

        Matrix3x3 beta_coefficient;
        double matrix3[3][3] = {{ax-R0x, ax-cx, Rdx},
                                {ay-R0y, ay-cy, Rdy},
                                {az-R0z, az-cz, Rdz}};
        beta_coefficient.copyMatrix(matrix3);
        double beta = beta_coefficient.determinant() / detA;

        Matrix3x3 gamma_coefficient;
        double matrix4[3][3] = {{ax-bx, ax-R0x, Rdx},
                                {ay-by, ay-R0y, Rdy},
                                {az-bz, az-R0z, Rdz}};
        gamma_coefficient.copyMatrix(matrix4);
        double gamma = gamma_coefficient.determinant() / detA;

        double t_min = -1.0;

        //beta+gamma < 1 and beta>0 and gamma>0 and t>0
        if( ! (beta+gamma < 1.0 && beta > 0 && gamma > 0 && detT > 0))
        {
            //return -1.0;
        }

        //return detT;
        else
            t_min = detT;

        if(level == 0)
            return t_min;

        Vector3D intersectPoint = r->start + r->dir * t_min;
        Color intersectPointColor = this->color;

        //ambient component
        this->pixelColor = intersectPointColor * coefficients[AMBIENT];
        color[0] = this->pixelColor.r;
        color[1] = this->pixelColor.g;
        color[2] = this->pixelColor.b;

        //cross product of (b-a) and (c-a)
        Vector3D normal = (point2 - point1) ^ (point3 - point1);

        for(auto pl : pointLights)
        {
            Vector3D rayl = intersectPoint - pl.light_pos;

            //find if rayl is obscured by any object
            Ray *ray = new Ray(pl.light_pos, rayl);
            double *dummyColor = new double[3];
            double tMin = 1e50;

            double curT = this->intersect(ray, dummyColor, 0);

            for(auto object : objects)
            {
                double iterT = object->intersect(ray, dummyColor, 0);

                if(iterT >= 0)
                    if(tMin > iterT)
                    {
                        tMin = min(tMin, iterT);
                    }
            }

            //another object is closer to the ray i.e. current object obscured
            if(tMin < curT || curT < 0)
            {
                //only ambient components need to be calculated
                return t_min;
            }

            //diffuse component
            double lambertValue = (normal % rayl) / (normal.distance() * rayl.distance());
            lambertValue = max(lambertValue, 0.0);

            this->pixelColor = this->pixelColor + pl.color * coefficients[DIFFUSE] * lambertValue ^ intersectPointColor;

            //specular component
            Vector3D rayv = intersectPoint - cameraPosition;
            Vector3D rayr = 2.0 * (rayl % normal) * normal - rayl;

            double phongValue = (rayv % rayr) / (rayv.distance() * rayr.distance());
            phongValue = max(pow(phongValue,this->shine*1.0), 0.0);

            this->pixelColor = this->pixelColor + pl.color * coefficients[SPECULAR] * phongValue ^ intersectPointColor;

            color[0] = this->pixelColor.r;
            color[1] = this->pixelColor.g;
            color[2] = this->pixelColor.b;
        }

        for(auto spl : spotLights)
        {
            Vector3D rayl = intersectPoint - spl.point_light.light_pos;

            //check it the angle between rayl and light_direction is smaller than cutoff angle
            double theta = (180.0 / pi) * cos((rayl % spl.light_direction) / (rayl.distance() * spl.light_direction.distance()));

            if(theta > spl.cutoff_angle)
            {
                return t_min;
            }

            //find if rayl is obscured by any object
            Ray *ray = new Ray(spl.point_light.light_pos, rayl);
            double *dummyColor = new double[3];
            double tMin = 1e50;

            double curT = this->intersect(ray, dummyColor, 0);

            for(auto object : objects)
            {
                double iterT = object->intersect(ray, dummyColor, 0);

                if(iterT >= 0)
                    if(tMin > iterT)
                    {
                        tMin = min(tMin, iterT);
                    }
            }

            //another object is closer to the ray i.e. current object obscured
            if(tMin < curT || curT < 0)
            {
                //only ambient components need to be calculated
                return t_min;
            }

            //diffuse component
            double lambertValue = (normal % rayl) / (normal.distance() * rayl.distance());
            lambertValue = max(lambertValue, 0.0);

            this->pixelColor = this->pixelColor + spl.point_light.color * coefficients[DIFFUSE] * lambertValue ^ intersectPointColor;

            //specular component
            Vector3D rayv = intersectPoint - cameraPosition;
            Vector3D rayr = 2.0 * (rayl % normal) * normal - rayl;

            double phongValue = (rayv % rayr) / (rayv.distance() * rayr.distance());
            phongValue = max(pow(phongValue,this->shine*1.0), 0.0);

            this->pixelColor = this->pixelColor + spl.point_light.color * coefficients[SPECULAR] * phongValue ^ intersectPointColor;

            color[0] = this->pixelColor.r;
            color[1] = this->pixelColor.g;
            color[2] = this->pixelColor.b;
        }

        return t_min;
    }
};

class GeneralQuadric : public Object
{
    double a, b, c, d, e, f, g, h, i, j;

public:
    GeneralQuadric(Vector3D reference_point, double length, double width, double height, double a, double b, double c, double d, double e, double f, double g, double h, double i, double j) : Object(reference_point, length, width, height)
    {
        this->a = a;
        this->b = b;
        this->c = c;
        this->d = d;
        this->e = e;
        this->f = f;
        this->g = g;
        this->h = h;
        this->i = i;
        this->j = j;
    }

    void draw()     {}

    bool shouldBeDrawn(Ray *r, double t)
    {
        Vector3D intersectingPoint = r->start + r->dir * t;

        if(length != 0.0)   //clip
        {
            if((intersectingPoint.x > reference_point.x + length) || (intersectingPoint.x < reference_point.x))
                return false;
        }

        if(width != 0.0)   //clip
        {
            if((intersectingPoint.y > reference_point.y + width) || (intersectingPoint.y < reference_point.y))
                return false;
        }

        if(height != 0.0)   //clip
        {
            if((intersectingPoint.z > reference_point.z + height) || (intersectingPoint.z < reference_point.z))
                return false;
        }

        return true;
    }

    double intersect(Ray *r, double *color, int level)
    {
        //equation: a*x^2 + b*y^2 + c*z^2 + d*x*y + e*y*z + f*x*z + g*x + h*y + i*z + j = 0
        Vector3D R0 = r->start;
        Vector3D Rd = r->dir;

        double rx = R0.x;
        double ry = R0.y;
        double rz = R0.z;
        double dx = Rd.x;
        double dy = Rd.y;
        double dz = Rd.z;

        color[0] = this->color.r;
        color[1] = this->color.g;
        color[2] = this->color.b;

        double t_min = -1.0;

        //coefficient of t^2
        double aa = a*dx*dx + b*dy*dy + c*dz*dz + d*dx*dy + e*dx*dz + f*dy*dz;

        //coefficient of t
        double bb = 2.0*a*dx*rx + 2.0*b*dy*ry + 2.0*c*dz*rz + d*rx*dy + d*dx*ry + e*dx*rz + e*dz*rx + f*dy*rz + f*ry*dz + g*dx + h*dy + i*dz;

        //constant
        double cc = a*rx*rx + b*ry*ry + c*rz*rz + d*rx*ry + e*rx*rz + f*ry*rz + g*rx + h*ry + i*rz + j;

        double dd_square = bb * bb - 4.0 * aa * cc;

        //if dd^2 < 0, ray does not intersect the surface
        if(dd_square >= 0.0)
        {
            double dd = sqrt(dd_square);

            double t1 = (-bb + dd) / (2.0 * aa);
            double t2 = (-bb - dd) / (2.0 * aa);

            //take closest positive
            if(t1 >= 0 && t2 >= 0)
            //check which one falls within the reference cube
            {
                double tMin = min(t1, t2);
                double tMax = max(t1, t2);

                if(shouldBeDrawn(r, tMin))
                    //return tMin;
                    t_min = tMin;

                else if(shouldBeDrawn(r, tMax))
                    //return tMax;
                    t_min = tMax;

                //return -1.0;
            }

            else if(t2 < 0.0)
            //check which one falls within the reference cube
            {
                if(shouldBeDrawn(r, t1))
                    //return t1;
                    t_min = t1;

                //return -1.0;
            }

            else if(t1 < 0.0)
            //check which one falls within the reference cube
            {
                if(shouldBeDrawn(r, t2))
                    //return t2;
                    t_min = t2;

                //return -1.0;
            }
        }

        //return -1.0;
        if(level == 0)
            return t_min;

        Vector3D intersectPoint = r->start + r->dir * t_min;
        Color intersectPointColor = this->color;

        //ambient component
        this->pixelColor = intersectPointColor * coefficients[AMBIENT];
        color[0] = this->pixelColor.r;
        color[1] = this->pixelColor.g;
        color[2] = this->pixelColor.b;

        double ix = intersectPoint.x;
        double iy = intersectPoint.y;
        double iz = intersectPoint.z;

        double nx = 2.0*a*ix + d*iy + e*iz + g;
        double ny = 2.0*b*iy + d*ix + f*iz + h;
        double nz = 2.0*c*iz + e*ix + f*iy + i;

        Vector3D normal(nx, ny, nz);

        for(auto pl : pointLights)
        {
            Vector3D rayl = intersectPoint - pl.light_pos;

            //find if rayl is obscured by any object
            Ray *ray = new Ray(pl.light_pos, rayl);
            double *dummyColor = new double[3];
            double tMin = 1e50;

            double curT = this->intersect(ray, dummyColor, 0);

            for(auto object : objects)
            {
                double iterT = object->intersect(ray, dummyColor, 0);

                if(iterT >= 0)
                    if(tMin > iterT)
                    {
                        tMin = min(tMin, iterT);
                    }
            }

            //another object is closer to the ray i.e. current object obscured
            if(tMin < curT || curT < 0)
            {
                //only ambient components need to be calculated
                return t_min;
            }

            //diffuse component
            double lambertValue = (normal % rayl) / (normal.distance() * rayl.distance());
            lambertValue = max(lambertValue, 0.0);

            this->pixelColor = this->pixelColor + pl.color * coefficients[DIFFUSE] * lambertValue ^ intersectPointColor;

            //specular component
            Vector3D rayv = intersectPoint - cameraPosition;
            Vector3D rayr = 2.0 * (rayl % normal) * normal - rayl;

            double phongValue = (rayv % rayr) / (rayv.distance() * rayr.distance());
            phongValue = max(pow(phongValue,this->shine*1.0), 0.0);

            this->pixelColor = this->pixelColor + pl.color * coefficients[SPECULAR] * phongValue ^ intersectPointColor;

            color[0] = this->pixelColor.r;
            color[1] = this->pixelColor.g;
            color[2] = this->pixelColor.b;
        }

        for(auto spl : spotLights)
        {
            Vector3D rayl = intersectPoint - spl.point_light.light_pos;

            //check it the angle between rayl and light_direction is smaller than cutoff angle
            double theta = (180.0 / pi) * cos((rayl % spl.light_direction) / (rayl.distance() * spl.light_direction.distance()));

            if(theta > spl.cutoff_angle)
            {
                return t_min;
            }

            //find if rayl is obscured by any object
            Ray *ray = new Ray(spl.point_light.light_pos, rayl);
            double *dummyColor = new double[3];
            double tMin = 1e50;

            double curT = this->intersect(ray, dummyColor, 0);

            for(auto object : objects)
            {
                double iterT = object->intersect(ray, dummyColor, 0);

                if(iterT >= 0)
                    if(tMin > iterT)
                    {
                        tMin = min(tMin, iterT);
                    }
            }

            //another object is closer to the ray i.e. current object obscured
            if(tMin < curT || curT < 0)
            {
                //only ambient components need to be calculated
                return t_min;
            }

            //diffuse component
            double lambertValue = (normal % rayl) / (normal.distance() * rayl.distance());
            lambertValue = max(lambertValue, 0.0);

            this->pixelColor = this->pixelColor + spl.point_light.color * coefficients[DIFFUSE] * lambertValue ^ intersectPointColor;

            //specular component
            Vector3D rayv = intersectPoint - cameraPosition;
            Vector3D rayr = 2.0 * (rayl % normal) * normal - rayl;

            double phongValue = (rayv % rayr) / (rayv.distance() * rayr.distance());
            phongValue = max(pow(phongValue,this->shine*1.0), 0.0);

            this->pixelColor = this->pixelColor + spl.point_light.color * coefficients[SPECULAR] * phongValue ^ intersectPointColor;

            color[0] = this->pixelColor.r;
            color[1] = this->pixelColor.g;
            color[2] = this->pixelColor.b;
        }

        return t_min;
    }
};

class Tile : public Object
{
    double tileX, tileY;
    double tileWidth;

public:
    Tile(double tileX, double tileY, double tileWidth)
    {
        //reference_point = Vector3D(-floorWidth/2.0, -floorWidth/2.0, 0);
        length = tileWidth;

        this->tileX = tileX;
        this->tileY = tileY;
        this->tileWidth = tileWidth;
    }

    void draw()
    {
        //glColor3f(color[0], color[1], color[2]);
        glColor3f(color.r, color.g, color.b);

        glBegin(GL_QUADS);
        {
            glVertex3f(tileX, tileY, 0);
            glVertex3f(tileX-tileWidth, tileY, 0);
            glVertex3f(tileX-tileWidth, tileY-tileWidth, 0);
            glVertex3f(tileX, tileY-tileWidth, 0);
        }
        glEnd();
    }

    double intersect(Ray *r, double *color, int level)
    {
        //plane: XY; point: (0,0,0), normal: z axis (0, 0, 1)
        //D = (0, 0, 1) dot (0, 0, 0) = 0
        //t = -(D + normal dot R0) / normal dot Rd
        if(r->dir.z == 0.0)
            return -1.0;

        double t = -(r->start.z) / (r->dir.z);

        //determine if the intersecting point is inside the tile/square
        //intersecting point = R0 + t*Rd
        Vector3D intersectPoint = r->start + r->dir * t;
        double x = intersectPoint.x;
        double y = intersectPoint.y;

        double xMin = tileX - tileWidth;
        double yMin = tileY - tileWidth;
        //since tile is on XY plane, compare x, y coordinates with those of the tile
        if((x >= xMin && x <= tileX) && (y >= yMin && y <= tileY))
            return t;

        return -1.0;
    }
};

class Floor : public Object
{
    double floorWidth, tileWidth;

public:
    Floor(double floorWidth, double tileWidth)
    {
        reference_point = Vector3D(-floorWidth/2.0, -floorWidth/2.0, 0);
        length = tileWidth;

        this->floorWidth = floorWidth;
        this->tileWidth = tileWidth;
    }

    void draw()
    {
        double floorOriginX = floorWidth/2.0;
        double floorOriginY = floorWidth/2.0;

        int tileNumbers = floorWidth / tileWidth;

        for(int i=0; i<tileNumbers; i++)
        {
            double x = floorOriginX - i * tileWidth;

            for(int j=0; j<tileNumbers; j++)
            {
                double y = floorOriginY - j * tileWidth;

                if(i % 2 == 0)
                {
                    if(j % 2 == 0)  glColor3f(0, 0, 0);
                    else    glColor3f(1, 1, 1);
                }
                else
                {
                    if(j % 2 == 1)  glColor3f(0, 0, 0);
                    else    glColor3f(1, 1, 1);
                }

                glBegin(GL_QUADS);
                {
                    glVertex3f(x, y, 0);
                    glVertex3f(x-tileWidth, y, 0);
                    glVertex3f(x-tileWidth, y-tileWidth, 0);
                    glVertex3f(x, y-tileWidth, 0);
                }
                glEnd();
            }
        }
    }

//    double intersect(Ray *r, double *color, int level)
//    {
//        //plane: XY; point: (0,0,0), normal: z axis (0, 0, 1)
//        //D = (0, 0, 1) dot (0, 0, 0) = 0
//        //t = -(D + normal dot R0) / normal dot Rd
//        if(r->dir.z == 0.0)
//            return -1.0;
//
//        double t = -(r->start.z) / (r->dir.z);
//
//        return t;
//    }

    double intersect(Ray *r, double *color, int level)
    {
        //plane: XY; point: (0,0,0), normal: z axis (0, 0, 1)
        //D = (0, 0, 1) dot (0, 0, 0) = 0
        //t = -(D + normal dot R0) / normal dot Rd
        if(r->dir.z == 0.0)
            return -1.0;

        double t = -(r->start.z) / (r->dir.z);

        //determine if the intersecting point is inside the tile/square
        //intersecting point = R0 + t*Rd
        Vector3D intersectPoint = r->start + r->dir * t;
        double x = intersectPoint.x;
        double y = intersectPoint.y;

        int i = int(floor((floorWidth / 2.0 - x) / tileWidth));
        int j = int(floor((floorWidth / 2.0 - y) / tileWidth));

//        if(i % 2 == 0)
//        {
//            if(j % 2 == 0)  this->setColor(0, 0, 0);
//            else    this->setColor(1, 1, 1);
//        }
//        else
//        {
//            if(j % 2 == 1)  this->setColor(0, 0, 0);
//            else    this->setColor(1, 1, 1);
//        }

        if(i % 2 == 0)
        {
            if(j % 2 == 0)
            {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            else
            {
                color[0] = 1;
                color[1] = 1;
                color[2] = 1;
            }
        }
        else
        {
            if(j % 2 == 1)
            {
                color[0] = 0;
                color[1] = 0;
                color[2] = 0;
            }
            else
            {
                color[0] = 1;
                color[1] = 1;
                color[2] = 1;
            }
        }

        return t;
    }
};


#endif // CLASSES_H_INCLUDED
