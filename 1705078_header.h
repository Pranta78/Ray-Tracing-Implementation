#ifndef HEADER_H_INCLUDED
#define HEADER_H_INCLUDED

#include "1705078_classes.h"

const string FOLDER = "";
extern int recursion_level;
int pixelNumbers;

//vector <Object*> objects;
//vector <PointLight> pointLights;
//vector <SpotLight> spotLights;

extern vector <Object*> objects;
extern vector <PointLight> pointLights;
extern vector <SpotLight> spotLights;

void loadData(string filename)
{
    string input;

    ifstream sceneFile;
    sceneFile.open(FOLDER + filename);

    if( !sceneFile.is_open() )
    {
        cout << "Cannot open file: " << filename << "\n";
        exit(-1);
    }

    //get line 1: level of recursion
    getline(sceneFile, input);
    cout << "input:" << input << "\n";
    istringstream iss1(input);
    iss1 >> recursion_level;

    //get line 2: pixels along both dimensions
    getline(sceneFile, input);
    cout << "input:" << input << "\n";
    istringstream iss2(input);
    iss2 >> pixelNumbers;

    getline(sceneFile, input);  //empty line

    //get line 3: number of objects
    int numObjects;
    getline(sceneFile, input);
    cout << "input:" << input << "\n";
    istringstream iss3(input);
    iss3 >> numObjects;

    for(int k=0; k<numObjects; k++)
    {
        //get object type
        getline(sceneFile, input);
        cout << "Object type:" << input << "\n";

        if( input.find("sphere") != string::npos )
        {
            //get sphere center
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss1(input);
            Vector3D center;
            iss1 >> center.x >> center.y >> center.z;

            //get sphere radius
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss2(input);
            double radius;
            iss2 >> radius;

            //get sphere color
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss3(input);
            double r, g, b;
            iss3 >> r >> g >> b;

            //get sphere ambient, diffuse, specular, recursive reflection coefficient
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss4(input);
            double ambient, diffuse, specular, recursion_reflection;
            iss4 >> ambient >> diffuse >> specular >> recursion_reflection;

            //get sphere shininess
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss5(input);
            double shininess;
            iss5 >> shininess;

            Object *temp;
            temp = new Sphere(center, radius);
            temp->setColor(r, g, b);
            temp->setShine(shininess);
            temp->setCoEfficients(ambient, diffuse, specular, recursion_reflection);
            //temp->draw();

            objects.push_back(temp);
        }

        else if( input.find("triangle") != string::npos )
        {
            //get triangle points
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss1(input);
            Vector3D point1;
            iss1 >> point1.x >> point1.y >> point1.z;

            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss2(input);
            Vector3D point2;
            iss2 >> point2.x >> point2.y >> point2.z;

            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss3(input);
            Vector3D point3;
            iss3 >> point3.x >> point3.y >> point3.z;

            //get triangle color
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss4(input);
            double r, g, b;
            iss4 >> r >> g >> b;

            //get triangle ambient, diffuse, specular, recursive reflection coefficient
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss5(input);
            double ambient, diffuse, specular, recursion_reflection;
            iss5 >> ambient >> diffuse >> specular >> recursion_reflection;

            //get triangle shininess
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss6(input);
            double shininess;
            iss6 >> shininess;

            Object *temp;
            temp = new Triangle(point1, point2, point3);
            temp->setColor(r, g, b);
            temp->setShine(shininess);
            temp->setCoEfficients(ambient, diffuse, specular, recursion_reflection);
            //temp->draw();

            objects.push_back(temp);
        }

        else if( input.find("general") != string::npos )
        {
            //get general quadric surface coefficients
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss1(input);
            double a, b, c, d, e, f, g, h, i, j;
            iss1 >> a >> b >> c >> d >> e >> f >> g >> h >> i >> j;

            //get cube reference point, length, width, height
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss2(input);
            Vector3D reference_point;
            double length, width, height;
            iss2 >> reference_point.x >> reference_point.y >> reference_point.z >> length >> width >> height;

            //get general quadric surface color
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss3(input);
            double rr, gg, bb;
            iss3 >> rr >> gg >> bb;

            //get general quadric surface ambient, diffuse, specular, recursive reflection coefficient
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss4(input);
            double ambient, diffuse, specular, recursion_reflection;
            iss4 >> ambient >> diffuse >> specular >> recursion_reflection;

            //get general quadric surface shininess
            getline(sceneFile, input);
            cout << "input:" << input << "\n";
            istringstream iss5(input);
            double shininess;
            iss5 >> shininess;

            Object *temp;
            temp = new GeneralQuadric(reference_point, length, width, height, a, b, c, d, e, f, g, h, i, j);
            temp->setColor(rr, gg, bb);
            temp->setShine(shininess);
            temp->setCoEfficients(ambient, diffuse, specular, recursion_reflection);

            objects.push_back(temp);
        }

        getline(sceneFile, input);  // empty line
    }

    //get point light sources
    int numPointLightSources;
    getline(sceneFile, input);
    cout << "numPointLightSources:" << input << "\n";
    istringstream iss4(input);
    iss4 >> numPointLightSources;

    for(int k=0; k<numPointLightSources; k++)
    {
        //get position of the light source
        getline(sceneFile, input);
        cout << "input:" << input << "\n";
        istringstream iss1(input);
        Vector3D position;
        iss1 >> position.x >> position.y >> position.z;

        //get color of the light source
        getline(sceneFile, input);
        cout << "input:" << input << "\n";
        istringstream iss2(input);
        double r, g, b;
        iss2 >> r >> g >> b;

        PointLight pl(position, r, g, b);

        pointLights.push_back(pl);
    }

    getline(sceneFile, input);  // empty line

    //get spot light sources
    int numSpotLightSources;
    getline(sceneFile, input);
    cout << "numSpotLightSources:" << input << "\n";
    istringstream iss5(input);
    iss5 >> numSpotLightSources;

    for(int k=0; k<numSpotLightSources; k++)
    {
        //get position of the spot light source
        getline(sceneFile, input);
        cout << "input:" << input << "\n";
        istringstream iss1(input);
        Vector3D position;
        iss1 >> position.x >> position.y >> position.z;

        //get color of the spot light source
        getline(sceneFile, input);
        cout << "input:" << input << "\n";
        istringstream iss2(input);
        double r, g, b;
        iss2 >> r >> g >> b;

        //get direction of the spot light source
        getline(sceneFile, input);
        cout << "input:" << input << "\n";
        istringstream iss3(input);
        Vector3D direction;
        iss3 >> direction.x >> direction.y >> direction.z;

        //get cutoff angle of the spot light source
        getline(sceneFile, input);
        cout << "input:" << input << "\n";
        istringstream iss4(input);
        double cutoff_angle;
        iss4 >> cutoff_angle;

        PointLight pl(position, r, g, b);
        SpotLight sl(pl, direction, cutoff_angle);

        spotLights.push_back(sl);
    }

    //the floor
    Object *temp;
    temp = new Floor(1000, 20);
    //temp->setCoEfficients(0.4, 0.2, 0.1, 0.1);
    temp->setCoEfficients(0.3, 0.4, 0.0, 0.3);
    //temp->setShine(20);
    temp->setShine(0.0);
    objects.push_back(temp);

    sceneFile.close();
}

#endif // HEADER_H_INCLUDED
