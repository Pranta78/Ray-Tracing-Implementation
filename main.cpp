#include<stdio.h>
#include<stdlib.h>
#include<math.h>

#include <windows.h>
#include <GL/glut.h>

#include <string>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>

using namespace std;

#include "1705078_header.h"
#include "bitmap_image.hpp"

#define root_2 sqrt(2.0)

#define windowHeight 500
#define windowWidth 500
#define viewAngle 90

double cameraHeight;
double cameraAngle;
int drawgrid;
int drawaxes;
double angle;

struct point
{
	double x,y,z;
};

//Vector3D cameraPosition = {100, 100, 0};
Vector3D cameraPosition = {128.482, 96.9838, 52.5};
//Vector3D cameraPosition = {63.4918, -92.1987, 39};
//Vector3D cameraPosition = {-37.3664, 107.725, 52.5};
//Vector3D cameraPosition = {151.476, 85.2762, 52.5};
//Vector3D cameraPosition = {-9.43898, -89.0654, 52.5};
//Vector3D cameraPosition = {8.73062, -72.1185, 45.0012};
Vector3D u = {0, 0, 1};
Vector3D r = {-1.0/root_2, 1.0/root_2, 0};
Vector3D l = {-1.0/root_2, -1.0/root_2, 0};

//Vector3D cameraPosition = {0, 0, 595.8767963};
//Vector3D u = {0, 1, 0};
//Vector3D r = {1, 0, 0};
//Vector3D l = {0, 0, -1};

Vector3D lookDirection = {cameraPosition.x + l.x, cameraPosition.y + l.y, cameraPosition.z + l.z};

double forwardMovementIncrement = 2.5;
double rightMovementIncrement = 1.5;
double upMovementIncrement = 1.5;

double lookupAngle = 0.1;
double lookRightAngle = 0.1;
double tiltAngle = 0.1;

double cameraRadius;

extern vector <Object*> objects;
extern vector <PointLight> pointLights;
extern vector <SpotLight> spotLights;

extern int recursion_level;
extern int pixelNumbers;

int imageWidth;
int imageHeight;

int counter = 11;

double getModulus(struct point Point)
{
    return sqrt(Point.x*Point.x + Point.y*Point.y + Point.z*Point.z);
}

double getModulus(Vector3D Point)
{
    return sqrt(Point.x*Point.x + Point.y*Point.y + Point.z*Point.z);
}

void capture()
{
    bitmap_image image;
    image.setwidth_height(pixelNumbers, pixelNumbers, true);

    double planeDistance = (windowHeight / 2.0) / tan(viewAngle * pi / (2.0 * 180));

    cout << "Plane Distance = " << planeDistance << "!\n";

    Vector3D ll(cameraPosition.x + l.x, cameraPosition.y + l.y, cameraPosition.z + l.z);
    Vector3D uu(u.x, u.y, u.z);
    Vector3D rr = ll ^ uu;

    Vector3D topleft = cameraPosition + l * planeDistance - r * windowWidth / 2.0 + u * windowHeight / 2.0;
    //Vector3D topleft = cameraPosition + lookDirection * planeDistance - rr * windowWidth / 2.0 + uu * windowHeight / 2.0;

    cout << "here -1!\n";

    imageWidth = pixelNumbers;
    imageHeight = pixelNumbers;

    double du = windowWidth * 1.0 / imageWidth;
    double dv = windowHeight * 1.0 / imageHeight;

    cout << "here asadsd!\n";

    //Choose middle of the grid cell
    topleft = topleft + r * (0.5 * du) - u * (0.5 * dv);
    //topleft = topleft + rr * (0.5 * du) - uu * (0.5 * dv);

    cout << "Initial topleft: x = " << topleft.x << ", y = " << topleft.y << ", z = " << topleft.z << "\n";

//    draw_circle = true;
//    cur_pixel = topleft;

    for(int i=0; i<imageWidth; i++)
    {
        for(int j=0; j<imageHeight; j++)
        {
            //topleft = topleft + r * (0.5 * du) * i - u * (0.5 * dv) * j;
            Vector3D curPixel = topleft + r * du * j - u * dv * i;
            //topleft = topleft + rr * (0.5 * du) * i - uu * (0.5 * dv) * j;

            //cout << "Curpixel: x = " << curPixel.x << ", y = " << curPixel.y << ", z = " << curPixel.z << "\n";
            //cur_pixel = curPixel;

            Ray *ray = new Ray(cameraPosition, curPixel - cameraPosition);
            double *dummyColor = new double[3];
            double nearestColor[3];

            Object *nearest = NULL;
            double t;
            double tMin = 1e50;

            //cout << "here 2!\n";

            for(auto object : objects)
            {
                double t = object->intersect(ray, dummyColor, 0);
                //cout << "here 2.1!\n";

//                if(t != -1)
//                    cout << "bingo! t = " << t << ", tMin = " << tMin << "\n";

                //save nearest object
                if(t >= 0)
                    if(tMin > t)
                    {
                        nearest = object;
                        nearestColor[0] = object->color[0];
                        nearestColor[1] = object->color[1];
                        nearestColor[2] = object->color[2];
                        tMin = min(tMin, t);
                        //cout << "Updated! t = " << t << ", tMin = " << tMin << "\n";
                        //image.set_pixel(j, i, 1, 0, 1);
                    }
            }

            if(nearest != NULL)
            {
                tMin = nearest -> intersect(ray, dummyColor, 1);
                nearestColor[0] = nearest->color[0];
                nearestColor[1] = nearest->color[1];
                nearestColor[2] = nearest->color[2];
                //cout << "Updating i = " << i << ", j = " << j << "\tr = " << nearest->color[0] << ", g = " << nearest->color[1] << ", b = " << nearest->color[2] << "\n";
                //image.set_pixel(j, i, nearest->color[0]*255.0, nearest->color[1]*255.0, nearest->color[2]*255.0);
                image.set_pixel(j, i, nearestColor[0]*255.0, nearestColor[1]*255.0, nearestColor[2]*255.0);
            }

            //image.set_pixel(j, i, 1, 0, 1);
        }

        cout << "here i = " << i << "!\n";
    }

    cout << "here 3!\n";

    //string filename = "output" + to_string(counter) + ".bmp";
    ostringstream ss;
    ss << counter;
    image.save_image(FOLDER+"output" + ss.str() + ".bmp");
    counter++;
}

void keyboardListener(unsigned char key, int x,int y){
	switch(key){

        case '0':
            {
                capture();
                break;
            }

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
            break;

		case GLUT_KEY_END:
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
				cout << "Camera Position: " << cameraPosition.x << ", " << cameraPosition.y << ", " << cameraPosition.z << "\n";
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

//	gluLookAt(cameraRadius*cos(cameraAngle), cameraRadius*sin(cameraAngle), cameraHeight,
//           0,0,0,		0,0,1);

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

	//cout << "Total objects: " << objects.size() << "\n";

	for(auto object : objects)
    {
        object -> draw();
    }

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
	drawgrid=1;
	drawaxes=0;
	cameraHeight=100.0;
	cameraAngle=1.0;
	cameraRadius=50;

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
	gluPerspective(viewAngle,	1,	1,	1000.0);
	//field of view in the Y (vertically)
	//aspect ratio that determines the field of view in the X direction (horizontally)
	//near distance
	//far distance
}

int main(int argc, char **argv){
	glutInit(&argc,argv);
	glutInitWindowSize(windowWidth, windowHeight);
	glutInitWindowPosition(0, 0);
	glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGB);	//Depth, Double buffer, RGB color

	glutCreateWindow("My OpenGL Program");

	loadData("scene.txt");

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
