#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

class Vector3D
{
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
    Vector3D reference_point;
    double height, width, length;
    double color[3];
    double coefficients[4]; //ambient, diffuse, specular, reflection coefficient
    int shine;  //exponent term of specular component

    Object(Vector3D reference_point)
    {
        this->reference_point = reference_point;
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

#endif // HEADER_H_INCLUDED
