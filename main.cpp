#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>

#define pi (2*acos(0.0))
#define root_2 sqrt(2.0)

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

struct point
{
	double x,y,z;
};

struct point cameraPosition = {100, 100, 0};
struct point u = {0, 0, 1};
struct point r = {-1.0/root_2, 1.0/root_2, 0};
struct point l = {-1.0/root_2, -1.0/root_2, 0};

struct point lookDirection = {cameraPosition.x + l.x, cameraPosition.y + l.y, cameraPosition.z + l.z};

double forwardMovementIncrement = 2.5;
double rightMovementIncrement = 0.5;
double upMovementIncrement = 1.5;

double lookupAngle = 0.1;
double lookRightAngle = 0.1;
double tiltAngle = 0.1;

//double sphereRadius = 15;
//double cylinderHeight = 20;
//how many steps or key presses to reach from sphere to cube and vice versa
const double sphereRadiusMax = 25;
const double cylinderHeightMax = 30;
const double steps = 25;

double sphereRadius = sphereRadiusMax * 15 / steps;
double cylinderHeight = cylinderHeightMax - cylinderHeightMax * 15 / steps;

double getModulus(struct point Point)
{
    return sqrt(Point.x*Point.x + Point.y*Point.y + Point.z*Point.z);
}

double getModulus(struct point Point1, struct point Point2)
{
    return sqrt((Point1.x-Point2.x) * (Point1.x-Point2.x)
                + (Point1.y-Point2.y) * (Point1.y-Point2.y)
                + (Point1.z-Point2.z) * (Point1.z-Point2.z));
}


void drawAxes()
{
    double lineLength = 200.0;

	if(drawaxes==1)
	{
		glColor3f(1.0, 1.0, 1.0);
		glBegin(GL_LINES);{
			glVertex3f( lineLength,0,0);
			glVertex3f(-lineLength,0,0);

			glVertex3f(0,-lineLength,0);
			glVertex3f(0, lineLength,0);

			glVertex3f(0,0, lineLength);
			glVertex3f(0,0,-lineLength);
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

//direction indicates along which axis values should be 0
//0 for x-axis, 1 for y, 2 for z
void drawSquare(double a, int direction)
{
    glColor3f(1.0,1.0,1.0);

    if(direction == 0)
    {
        glBegin(GL_QUADS);{
            glVertex3f(0, a/2, a/2);
            glVertex3f(0, a/2, -a/2);
            glVertex3f(0, -a/2, -a/2);
            glVertex3f(0, -a/2, a/2);
        }glEnd();
    }

	else if(direction == 1)
    {
        glBegin(GL_QUADS);{
            glVertex3f(a/2, 0, a/2);
            glVertex3f(a/2, 0, -a/2);
            glVertex3f(-a/2, 0, -a/2);
            glVertex3f(-a/2, 0, a/2);
        }glEnd();
    }

    else if(direction == 2)
    {
        glBegin(GL_QUADS);{
            glVertex3f(a/2, a/2, 0);
            glVertex3f(a/2, -a/2, 0);
            glVertex3f(-a/2, -a/2, 0);
            glVertex3f(-a/2, a/2, 0);
        }glEnd();
    }
}


void drawCircle(double radius,int segments)
{
    int i;
    struct point points[segments+1];
    glColor3f(0.7,0.7,0.7);
    //generate points
    for(i=0;i<=segments;i++)    // radius is drawn as well if only "<" is used
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
    struct point points[segments+1];
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
	// struct point points[100][100];
	struct point points[stacks+1][slices+1];
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


void drawOctetSphere(double radius,int slices,int stacks, int x, int y, int z)
{
    glColor3f(1, 0, 0);

	// struct point points[100][100];
	struct point points[stacks+1][slices+1];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=radius*sin(((double)i/(double)stacks)*(pi/2));
		r=radius*cos(((double)i/(double)stacks)*(pi/2));
		for(j=0;j<=slices;j++)
		{
			points[i][j].x=x * r*cos(((double)j/((double)slices * 2.0))*pi);
			points[i][j].y=y * r*sin(((double)j/((double)slices * 2.0))*pi);
			points[i][j].z=z * h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
			    //upper hemisphere
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}


void drawCylinder(double radius,int slices,int stacks)
{
	// struct point points[100][100];
	struct point points[stacks+1][slices+1];
	int i,j;
	double h,r;
	//generate points
	for(i=0;i<=stacks;i++)
	{
		h=i;
		r=radius;
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
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}


void drawOneFourthCylinder(double radius,int slices,int stacks, int x, int y, int z)
{
    glColor3f(0, 1, 0);

	// struct point points[100][100];
	struct point points[stacks+1][slices+1];
	int i,j;
	double h,r;
	//generate points
//	for(i=0;i<=stacks;i++)
//	{
//		h=i;
//		r=radius;
//		for(j=0;j<=slices;j++)
//		{
//			points[i][j].x = x * r*cos(((double)j/((double)slices*2.0))*pi);
//			points[i][j].y = y * r*sin(((double)j/((double)slices*2.0))*pi);
//			points[i][j].z = z * h;
//		}
//	}

    for(i=0;i<=stacks;i++)
	{
		h=i - stacks/2;
		r=radius;
		for(j=0;j<=slices;j++)
		{
			points[i][j].x = x * r*cos(((double)j/((double)slices*2.0))*pi);
			points[i][j].y = y * r*sin(((double)j/((double)slices*2.0))*pi);
			points[i][j].z = z * h;
		}
	}
	//draw quads using generated points
	for(i=0;i<stacks;i++)
	{
        //glColor3f((double)i/(double)stacks,(double)i/(double)stacks,(double)i/(double)stacks);
		for(j=0;j<slices;j++)
		{
			glBegin(GL_QUADS);{
				glVertex3f(points[i][j].x,points[i][j].y,points[i][j].z);
				glVertex3f(points[i][j+1].x,points[i][j+1].y,points[i][j+1].z);
				glVertex3f(points[i+1][j+1].x,points[i+1][j+1].y,points[i+1][j+1].z);
				glVertex3f(points[i+1][j].x,points[i+1][j].y,points[i+1][j].z);
			}glEnd();
		}
	}
}

void drawSphereCube()
{
    //Translate the sphere slices by sphereTranslation along each axis

    //draw the one-eight spheres in 8 places around the origin
    //translation co-ordinates are of form (+/- st, +/- st, +/- st)
    for(int i=0; i<8; i++)
    {
        int x = (i & 0b100) ? 1 : -1;
        int y = (i & 0b10) ? 1 : -1;
        int z = (i & 0b1) ? 1 : -1;

        glTranslatef(x*cylinderHeight, y*cylinderHeight, z*cylinderHeight);
        drawOctetSphere(sphereRadius, 100, 100, x, y, z);
        glTranslatef(-x*cylinderHeight, -y*cylinderHeight, -z*cylinderHeight);
    }

    //draw the one-fourth cylinder in 12 places around the origin
    //translation co-ordinates are of form (+/- st, +/- st, 0) -> 0 in exactly one of the 3 places
    for(int i=0; i<16; i++)
    {
        int x, y, z;

        if(i/4 == 0)
        {
            x = 0;
            y = (i & 0b10) ? 1 : -1;
            z = (i & 0b1) ? 1 : -1;

            if(y == -1)
            {
                glRotatef(90.0, 1, 0, 0);
                glTranslatef(cylinderHeight, -cylinderHeight, 0);
            }
            else
            {
                glRotatef(90.0, -1, 0, 0);
                glRotatef(180.0, 0, 1, 0);
                glTranslatef(cylinderHeight, cylinderHeight, 0);
            }
        }

        else if(i/4 == 1)
        {
            x = (i & 0b10) ? 1 : -1;
            y = 0;
            z = (i & 0b1) ? 1 : -1;

            glRotatef(90.0, 1, 0, 0);

            if(x == -1)
                glTranslatef(-cylinderHeight, cylinderHeight, 0);
            else
                glTranslatef(cylinderHeight, cylinderHeight, 0);
        }

        else if(i/4 == 2)
        {
            x = (i & 0b10) ? 1 : -1;
            y = (i & 0b1) ? 1 : -1;
            z = 0;

            glTranslatef(x*cylinderHeight, y*cylinderHeight, z*cylinderHeight);
        }

        else if(i/4 == 3)
        {
            x = (i & 0b10) ? 1 : -1;
            y = (i & 0b1) ? 1 : -1;
            z = 0;

            glRotatef(90.0, 0, 1, 0);
            glTranslatef(x*cylinderHeight, y*cylinderHeight, 0);
        }

        //printf("x=%d, y=%d, z=%d\n", x, y, z);

        //glTranslatef(x*sphereTranslation, y*sphereTranslation, z*sphereTranslation);

        int newX = (x == 0) ? 1 : x;
        int newY = (y == 0) ? 1 : y;
        int newZ = (z == 0) ? 1 : z;

        drawOneFourthCylinder(sphereRadius, 100, 2*cylinderHeight, newX, newY, newZ);

        //glTranslatef(-x*sphereTranslation, -y*sphereTranslation, -z*sphereTranslation);

        if(i/4 == 0)
        {
            if(y == -1)
            {
                glTranslatef(-cylinderHeight, cylinderHeight, 0);
                glRotatef(90.0, -1, 0, 0);
            }
            else
            {
                glTranslatef(-cylinderHeight, -cylinderHeight, 0);
                glRotatef(180.0, 0, -1, 0);
                glRotatef(90.0, 1, 0, 0);
                //glTranslatef(0, -sphereTranslation, 0);
            }

            //glRotatef(-90.0, 1, 0, 0);
        }

        else if(i/4 == 1)
        {
            if(x == -1)
                glTranslatef(cylinderHeight, -cylinderHeight, 0);
            else
                glTranslatef(-cylinderHeight, -cylinderHeight, 0);

            glRotatef(-90.0, 1, 0, 0);
        }

        else if(i/4 == 2)
        {
            glTranslatef(-x*cylinderHeight, -y*cylinderHeight, -z*cylinderHeight);
        }

        else if(i/4 == 3)
        {
            glTranslatef(-x*cylinderHeight, -y*cylinderHeight, 0);
            glRotatef(90.0, 0, -1, 0);
        }
    }

    //draw the squares in 6 places around the origin
    //translation co-ordinates are of form (+/- st, 0, 0) with exactly 2 axis as value zero
    for(int i=0; i<6; i++)
    {
        int x, y, z;

        if(i/2 == 0)
        {
            x = (i & 1) ? 1 : -1;
            y = 0;
            z = 0;
        }
        else if(i/2 == 1)
        {
            x = 0;
            y = (i & 1) ? 1 : -1;
            z = 0;
        }
        else if(i/2 == 2)
        {
            x = 0;
            y = 0;
            z = (i & 1) ? 1 : -1;
        }

        glTranslatef(x*(cylinderHeight+sphereRadius), y*(cylinderHeight+sphereRadius), z*(cylinderHeight+sphereRadius));

        //side length of the square = height of the cylinder + radius
        drawSquare(2*cylinderHeight , i/2);

        glTranslatef(-x*(cylinderHeight+sphereRadius), -y*(cylinderHeight+sphereRadius), -z*(cylinderHeight+sphereRadius));
    }
}

void drawSS()
{
    glColor3f(1,0,0);
    //drawSquare(20);

    glRotatef(angle,0,0,1);
    glTranslatef(110,0,0);
    glRotatef(2*angle,0,0,1);
    glColor3f(0,1,0);
    //drawSquare(15);
/*
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
*/
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

        //look left
        //rotate r and l with respect to u
		case '1':
			{
			    struct point newL;
			    // u cross l is -r
			    newL.x = l.x * cos(lookRightAngle) - r.x * sin(lookRightAngle);
			    newL.y = l.y * cos(lookRightAngle) - r.y * sin(lookRightAngle);
			    newL.z = l.z * cos(lookRightAngle) - r.z * sin(lookRightAngle);

			    struct point newR;
			    //u cross r is l
			    newR.x = r.x * cos(lookRightAngle) + l.x * sin(lookRightAngle);
			    newR.y = r.y * cos(lookRightAngle) + l.y * sin(lookRightAngle);
			    newR.z = r.z * cos(lookRightAngle) + l.z * sin(lookRightAngle);

			    l.x = newL.x;
			    l.y = newL.y;
			    l.z = newL.z;

			    r.x = newR.x;
			    r.y = newR.y;
			    r.z = newR.z;

			    break;
			}

        //look right
		case '2':
			{
			    struct point newL;
			    // u cross l is -r
			    newL.x = l.x * cos(-lookRightAngle) - r.x * sin(-lookRightAngle);
			    newL.y = l.y * cos(-lookRightAngle) - r.y * sin(-lookRightAngle);
			    newL.z = l.z * cos(-lookRightAngle) - r.z * sin(-lookRightAngle);

			    struct point newR;
			    //u cross r is l
			    newR.x = r.x * cos(-lookRightAngle) + l.x * sin(-lookRightAngle);
			    newR.y = r.y * cos(-lookRightAngle) + l.y * sin(-lookRightAngle);
			    newR.z = r.z * cos(-lookRightAngle) + l.z * sin(-lookRightAngle);

			    l.x = newL.x;
			    l.y = newL.y;
			    l.z = newL.z;

			    r.x = newR.x;
			    r.y = newR.y;
			    r.z = newR.z;

			    break;
			}

        //look up
        //rotate u and l with respect to r
        case '3':
			{
			    struct point newL;
			    // u is r cross l
			    newL.x = l.x * cos(lookupAngle) + u.x * sin(lookupAngle);
			    newL.y = l.y * cos(lookupAngle) + u.y * sin(lookupAngle);
			    newL.z = l.z * cos(lookupAngle) + u.z * sin(lookupAngle);

			    struct point newU;
			    //-l is r cross u
			    newU.x = u.x * cos(lookupAngle) - l.x * sin(lookupAngle);
			    newU.y = u.y * cos(lookupAngle) - l.y * sin(lookupAngle);
			    newU.z = u.z * cos(lookupAngle) - l.z * sin(lookupAngle);

			    l.x = newL.x;
			    l.y = newL.y;
			    l.z = newL.z;

			    u.x = newU.x;
			    u.y = newU.y;
			    u.z = newU.z;

//			    lookDirection.x = cameraPosition.x + l.x;
//			    lookDirection.y = cameraPosition.y + l.y;
//			    lookDirection.z = cameraPosition.z + l.z;
			    break;
			}

        //look down
        case '4':
			{
			    struct point newL;
			    newL.x = l.x * cos(-lookupAngle) + u.x * sin(-lookupAngle);
			    newL.y = l.y * cos(-lookupAngle) + u.y * sin(-lookupAngle);
			    newL.z = l.z * cos(-lookupAngle) + u.z * sin(-lookupAngle);

			    struct point newU;
			    //-l is r cross u
			    newU.x = u.x * cos(-lookupAngle) - l.x * sin(-lookupAngle);
			    newU.y = u.y * cos(-lookupAngle) - l.y * sin(-lookupAngle);
			    newU.z = u.z * cos(-lookupAngle) - l.z * sin(-lookupAngle);

			    l.x = newL.x;
			    l.y = newL.y;
			    l.z = newL.z;

			    u.x = newU.x;
			    u.y = newU.y;
			    u.z = newU.z;

//			    lookDirection.x = cameraPosition.x + l.x;
//			    lookDirection.y = cameraPosition.y + l.y;
//			    lookDirection.z = cameraPosition.z + l.z;
			    break;
			}

        //tilt clockwise
        //rotate u and r with respect to l
        case '5':
			{
			    struct point newR;
			    //l cross r is -u
			    newR.x = r.x * cos(tiltAngle) - u.x * sin(tiltAngle);
			    newR.y = r.y * cos(tiltAngle) - u.y * sin(tiltAngle);
			    newR.z = r.z * cos(tiltAngle) - u.z * sin(tiltAngle);

			    struct point newU;
			    //l cross u is r
			    newU.x = u.x * cos(tiltAngle) + r.x * sin(tiltAngle);
			    newU.y = u.y * cos(tiltAngle) + r.y * sin(tiltAngle);
			    newU.z = u.z * cos(tiltAngle) + r.z * sin(tiltAngle);

			    r.x = newR.x;
			    r.y = newR.y;
			    r.z = newR.z;

			    u.x = newU.x;
			    u.y = newU.y;
			    u.z = newU.z;

			    break;
			}

        //tilt anticlockwise
        case '6':
			{
			    struct point newR;
			    //l cross r is -u
			    newR.x = r.x * cos(-tiltAngle) - u.x * sin(-tiltAngle);
			    newR.y = r.y * cos(-tiltAngle) - u.y * sin(-tiltAngle);
			    newR.z = r.z * cos(-tiltAngle) - u.z * sin(-tiltAngle);

			    struct point newU;
			    //l cross u is r
			    newU.x = u.x * cos(-tiltAngle) + r.x * sin(-tiltAngle);
			    newU.y = u.y * cos(-tiltAngle) + r.y * sin(-tiltAngle);
			    newU.z = u.z * cos(-tiltAngle) + r.z * sin(-tiltAngle);

			    r.x = newR.x;
			    r.y = newR.y;
			    r.z = newR.z;

			    u.x = newU.x;
			    u.y = newU.y;
			    u.z = newU.z;

			    break;
			}

		default:
			break;
	}
}


void specialKeyListener(int key, int x,int y){
	switch(key){
		case GLUT_KEY_DOWN:		//down arrow key
			{
		        //cameraHeight -= 3.0;
                double dist = getModulus(l);
                //double dist = getModulus(l, cameraPosition);
                cameraPosition.x -= (forwardMovementIncrement * l.x / dist);
                cameraPosition.y -= (forwardMovementIncrement * l.y / dist);
                cameraPosition.z -= (forwardMovementIncrement * l.z / dist);
                break;
		    }
		case GLUT_KEY_UP:		// up arrow key
		    {
		        //cameraHeight += 3.0;
                double dist = getModulus(l);
                //double dist = getModulus(l, cameraPosition);
                cameraPosition.x += (forwardMovementIncrement * l.x / dist);
                cameraPosition.y += (forwardMovementIncrement * l.y / dist);
                cameraPosition.z += (forwardMovementIncrement * l.z / dist);
                break;
		    }

		case GLUT_KEY_RIGHT:
			{
		        //cameraAngle += 0.03;
                double dist = getModulus(r);
                //double dist = getModulus(r, cameraPosition);
                cameraPosition.x += (rightMovementIncrement * r.x / dist);
                cameraPosition.y += (rightMovementIncrement * r.y / dist);
                cameraPosition.z += (rightMovementIncrement * r.z / dist);
                break;
		    }
		case GLUT_KEY_LEFT:
			{
		        //cameraAngle -= 0.03;
                double dist = getModulus(r);
                //double dist = getModulus(r, cameraPosition);
                cameraPosition.x -= (rightMovementIncrement * r.x / dist);
                cameraPosition.y -= (rightMovementIncrement * r.y / dist);
                cameraPosition.z -= (rightMovementIncrement * r.z / dist);
                break;
		    }

		case GLUT_KEY_PAGE_UP:
			{
                double dist = getModulus(u);
                //double dist = getModulus(u, cameraPosition);
                cameraPosition.x += (upMovementIncrement * u.x / dist);
                cameraPosition.y += (upMovementIncrement * u.y / dist);
                cameraPosition.z += (upMovementIncrement * u.z / dist);
                break;
		    }
		case GLUT_KEY_PAGE_DOWN:
			{
                double dist = getModulus(u);
                //double dist = getModulus(u, cameraPosition);
                cameraPosition.x -= (upMovementIncrement * u.x / dist);
                cameraPosition.y -= (upMovementIncrement * u.y / dist);
                cameraPosition.z -= (upMovementIncrement * u.z / dist);
                break;
		    }

		case GLUT_KEY_INSERT:
			break;

		case GLUT_KEY_HOME:
			{
			    //increase sphere radius and decrease cylinder height
			    if(sphereRadius + sphereRadiusMax / steps < sphereRadiusMax)
                    //sphereRadius += 1;
                    sphereRadius += sphereRadiusMax / steps;
                else
                    sphereRadius = sphereRadiusMax;

                if(cylinderHeight > 1)
                    //cylinderHeight -= 2;
                    cylinderHeight -= cylinderHeightMax / steps;
                else
                    cylinderHeight = 0;

			    break;
			}
		case GLUT_KEY_END:
			{
			    //decrease sphere radius and increase cylinder height
			    if(sphereRadius > 0)
                    //sphereRadius -= 1;
                    sphereRadius -= sphereRadiusMax / steps;
                else
                    sphereRadius = 0;

                if(cylinderHeight + cylinderHeightMax / steps < cylinderHeightMax)
                    //cylinderHeight += 2;
                    cylinderHeight += cylinderHeightMax / steps;
                else
                    cylinderHeight = cylinderHeightMax;

			    break;
			}

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
	//gluLookAt(0,200,0,	0,0,0,	1,0,0);

	gluLookAt(cameraPosition.x, cameraPosition.y, cameraPosition.z,
           cameraPosition.x + l.x, cameraPosition.y + l.y, cameraPosition.z + l.z,
           //lookDirection.x, lookDirection.y, lookDirection.z,
           u.x, u.y, u.z);


	//again select MODEL-VIEW
	glMatrixMode(GL_MODELVIEW);


	/****************************
	/ Add your objects from here
	****************************/
	//add objects

	drawAxes();
	drawGrid();

    //glColor3f(1,0,0);
    //drawSquare(10);

    //drawSS();

    //drawCircle(30,99);

    //drawCone(20,100,24);

	//drawSphere(50,100,100);

	//drawOctetSphere(50,100,100);

    //drawCylinder(30,100,70);

    //drawOneFourthCylinder(30,100,70, 0, 1, 1);

    drawSphereCube();

	//ADD this line in the end --- if you use double buffer (i.e. GL_DOUBLE)
	glutSwapBuffers();
}


void animate(){
	angle+=(0.05 * 180 / pi);
	//codes for any changes in Models, Camera
	glutPostRedisplay();
}

void init(){
	//codes for initialization
	drawgrid=0;
	drawaxes=1;
	cameraHeight=150.0;
	cameraAngle=1.0;
	angle=0;

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
	gluPerspective(80,	1,	1,	1000.0);
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
