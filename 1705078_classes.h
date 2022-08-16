#ifndef CLASSES_H_INCLUDED
#define CLASSES_H_INCLUDED

#define pi (2*acos(0.0))

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

    friend Vector3D operator*(double a, const Vector3D& v)
    {
        return Vector3D(a*v.x, a*v.y, a*v.z);
    }
};

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
    double color[3];
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
        color[0] = red;
        color[1] = green;
        color[2] = blue;
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
        cout << "Object Color: r = " << color[0] << ", g = " << color[1] << ", b = " << color[2] << "\n";
    }
};

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

            glColor3f(color[0], color[1], color[2]);

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

        //if d^2 < 0, ray does not intersect the sphere
        if(d_square >= 0.0)
        {
            double d = sqrt(d_square);

            double t1 = (-b + d) / (2.0 * a);
            double t2 = (-b - d) / (2.0 * a);

            //take closest positive
            if(t1 >= 0 && t2 >= 0)
                return min(t1, t2);

            if(t2 < 0.0)
                return t1;

            if(t1 < 0.0)
                return t2;
        }

        return -1.0;
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
        glColor3f(color[0], color[1], color[2]);

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

        //beta+gamma < 1 and beta>0 and gamma>0 and t>0
        if( ! (beta+gamma < 1.0 && beta > 0 && gamma > 0 && detT > 0))
        {
            return -1.0;
        }

        return detT;
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

        //coefficient of t^2
        double aa = a*dx*dx + b*dy*dy + c*dz*dz + d*dx*dy + e*dy*dz + f*dz*dx;

        //coefficient of t
        double bb = 2.0*a*dx*rx + 2.0*b*dy*ry + 2.0*c*dz*rz + d*rx*dy + d*dx*ry + e*dy*rz + e*ry*dz + f*dx*rz + f*dz*rx + g*dx + h*dy + i*dz;

        //constant
        double cc = a*rx*rx + b*ry*ry + c*rz*rz + d*rx*ry + e*ry*rz + f*rx*rz + g*rx + h*ry + i*rz + j;

        double dd_square = bb * bb - 4.0 * aa * cc;

        //if dd^2 < 0, ray does not intersect the surface
        if(dd_square >= 0.0)
        {
            double dd = sqrt(dd_square);

            double t1 = (-bb + dd) / (2.0 * aa);
            double t2 = (-bb - dd) / (2.0 * aa);

            //take closest positive
            if(t1 >= 0 && t2 >= 0)
                return min(t1, t2);

            if(t2 < 0.0)
                return t1;

            if(t1 < 0.0)
                return t2;
        }

        return -1.0;
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
        glColor3f(color[0], color[1], color[2]);

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

        if(i % 2 == 0)
        {
            if(j % 2 == 0)  this->setColor(0, 0, 0);
            else    this->setColor(1, 1, 1);
        }
        else
        {
            if(j % 2 == 1)  this->setColor(0, 0, 0);
            else    this->setColor(1, 1, 1);
        }

        return t;
    }
};

class PointLight
{
    Vector3D light_pos;
    double color[3];

public:
    PointLight()    {}

    PointLight(Vector3D light_pos, double r, double g, double b)
    {
        this->light_pos = light_pos;
        this->color[0] = r;
        this->color[1] = g;
        this->color[2] = b;
    }
};

class SpotLight
{
    PointLight point_light;
    Vector3D light_direction;
    double cutoff_angle;

public:
    SpotLight(PointLight point_light, Vector3D direction, double cutoff_angle)
    {
        this->point_light = point_light;
        this->light_direction = direction;
        this->cutoff_angle = cutoff_angle;
    }
};



#endif // CLASSES_H_INCLUDED
