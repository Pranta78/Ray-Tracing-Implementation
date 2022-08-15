#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

#include "1705078_classes.h"

const string FOLDER = "";
int recursion_level;
int pixelNumbers;

vector <Object> objects;
vector <PointLight> pointLights;
vector <SpotLight> spotLights;

void loadData(string filename)
{
    string input;

    ifstream sceneFile;
    sceneFile.open(FOLDER + "scene.txt");

    if( !sceneFile.is_open() )
    {
        cout << "Cannot open file: scene.txt\n";
        exit(-1);
    }

    //get line 1: level of recursion
    getline(sceneFile, input)
    cout << "input:" << input << "\n";
    istringstream iss(input);
    iss >> recursion_level;

    //get line 2: pixels along both dimensions
    getline(sceneFile, input)
    cout << "input:" << input << "\n";
    istringstream iss(input);
    iss >> pixelNumbers;

    //get line 3: number of objects
    int numObjects;
    getline(sceneFile, input)
    cout << "input:" << input << "\n";
    istringstream iss(input);
    iss >> numObjects;

    for(int k=0; k<numObjects; k++)
    {
        //get object type
        getline(sceneFile, input)
        cout << "Object type:" << input << "\n";

        if(input == "sphere")
        {
            //get sphere center
            getline(sceneFile, input);
            istringstream iss(input);
            Vector3D center;
            iss >> center.x >> center.y >> center.z;

            //get sphere radius
            getline(sceneFile, input);
            istringstream iss(input);
            double radius;
            iss >> radius;

            //get sphere color
            getline(sceneFile, input);
            istringstream iss(input);
            double r, g, b;
            iss >> r >> g >> b;

            //get sphere ambient, diffuse, specular, recursive reflection coefficient
            getline(sceneFile, input);
            istringstream iss(input);
            double ambient, diffuse, specular, recursion_reflection;
            iss >> ambient >> diffuse >> specular >> recursion_reflection;

            //get sphere shininess
            getline(sceneFile, input);
            istringstream iss(input);
            double shininess;
            iss >> shininess;

            Object *temp;
            temp = new Sphere(center, radius);
            temp->setColor(r, g, b);
            temp->setShine(shininess);
            temp->setCoEfficients(ambient, diffuse, specular, recursion_reflection);

            objects.push_back(temp);
        }

        else if(input == "triangle")
        {
            //get triangle points
            getline(sceneFile, input);
            istringstream iss(input);
            Vector3D point1;
            iss >> point1.x >> point1.y >> point1.z;

            getline(sceneFile, input);
            istringstream iss(input);
            Vector3D point2;
            iss >> point2.x >> point2.y >> point2.z;

            getline(sceneFile, input);
            istringstream iss(input);
            Vector3D point3;
            iss >> point3.x >> point3.y >> point3.z;

            //get triangle color
            getline(sceneFile, input);
            istringstream iss(input);
            double r, g, b;
            iss >> r >> g >> b;

            //get triangle ambient, diffuse, specular, recursive reflection coefficient
            getline(sceneFile, input);
            istringstream iss(input);
            double ambient, diffuse, specular, recursion_reflection;
            iss >> ambient >> diffuse >> specular >> recursion_reflection;

            //get triangle shininess
            getline(sceneFile, input);
            istringstream iss(input);
            double shininess;
            iss >> shininess;

            Object *temp;
            temp = new Triangle(point1, point2, point3);
            temp->setColor(r, g, b);
            temp->setShine(shininess);
            temp->setCoEfficients(ambient, diffuse, specular, recursion_reflection);

            objects.push_back(temp);
        }

        else if(input == "general")
        {
            //get general quadric surface coefficients
            getline(sceneFile, input);
            istringstream iss(input);
            double a, b, c, d, e, f, g, h, i, j;
            iss >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j;

            //get cube reference point, length, width, height
            getline(sceneFile, input);
            istringstream iss(input);
            Vector3D reference_point;
            double length, width, height;
            iss >> reference_point.x >> reference_point.y >> reference_point.z >> length >> width >> height;

            //get general quadric surface color
            getline(sceneFile, input);
            istringstream iss(input);
            double r, g, b;
            iss >> r >> g >> b;

            //get general quadric surface ambient, diffuse, specular, recursive reflection coefficient
            getline(sceneFile, input);
            istringstream iss(input);
            double ambient, diffuse, specular, recursion_reflection;
            iss >> ambient >> diffuse >> specular >> recursion_reflection;

            //get general quadric surface shininess
            getline(sceneFile, input);
            istringstream iss(input);
            double shininess;
            iss >> shininess;

            Object *temp;
            temp = new GeneralQuadric(reference_point, length, width, height, a, b, c, d, e, f, g, h, i, j);
            temp->setColor(r, g, b);
            temp->setShine(shininess);
            temp->setCoEfficients(ambient, diffuse, specular, recursion_reflection);

            objects.push_back(temp);
        }
    }

    while (getline(sceneFile, input))
	{
	    cout << "input:" << input << "\n";
		istringstream iss(input);
		iss >> up.x >> up.y >> up.z;
		break;
	}

    //get line 4: fovY aspectRatio near far
    while (getline(sceneFile, input))
	{
	    cout << "input:" << input << "\n";
		istringstream iss(input);
		iss >> fovYValue >> aspectRatioValue >> nearValue >> farValue;
		break;
	}
}

#endif // HEADER_H_INCLUDED
