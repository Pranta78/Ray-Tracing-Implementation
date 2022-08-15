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
};

class Object
{
protected:
    Vector3D reference_point;
    double height, width, length;
    double color[3];
    double coefficients[4]; //ambient, diffuse, specular, reflection coefficient
    int shine;  //exponent term of specular component

public:
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
            int x = floorOriginX - i * tileWidth;

            for(int j=0; j<tileNumbers; j++)
            {
                int y = floorOriginY - j * tileWidth;

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
