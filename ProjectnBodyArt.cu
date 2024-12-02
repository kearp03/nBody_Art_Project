// nvcc nBodyArtB.cu -o nBodyArt -lglut -lm -lGLU -lGL
//To stop hit "control c" in the window you launched it from.
#include <iostream>
#include <fstream>
#include <sstream>
#include <string.h>
#include <GL/glut.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <cuda.h>
using namespace std;

#define PI 3.141592653589793238462643383279502884197169399375105820974944592307816406286

FILE* ffmpeg;

// defines for terminal stuff.
#define BOLD_ON  "\e[1m"
#define BOLD_OFF   "\e[m"

typedef struct {
	int id;
	float mass;
	float4 pos;
	float4 vel;
	float4 force;
	float4 color;
} Body;

FILE* MovieFile;

// Globals
int NumberOfBodies;
float TotalRunTime;
float Dt;
float G;
float H;
float dForce;
float Epsilon;
float MassOfBody;
float DiameterOfBody;
float VelocityMax;
float Drag;
int DrawRate;
int PrintRate;
// int NumberOfBodies;
int Capacity = 100;

// Other Globals
int Pause;
int LClickOn = 0;
int DrawTimer, PrintTimer;
float RunTime;
int* Buffer;
int MovieOn;
int MovieFlag;
int Trace;
double MouseX, MouseY, MouseZ;
float4 NextColor = {0.0f,0.0f,0.0f,1.0f};
string NextColorString = "Random";
int RandomColor = 1;
int GToggle = 0;
int HToggle = 0;

// Window globals
static int Window;
int XWindowSize;
int YWindowSize; 
double Near;
double Far;
double EyeX;
double EyeY;
double EyeZ;
double CenterX;
double CenterY;
double CenterZ;
double UpX;
double UpY;
double UpZ;

// Prototyping functions
void addBody(Body);
void setSimulationParameters();
void allocateMemory();
void setInitailConditions();
void drawPicture();
void nBody();
void errorCheck(const char*);
void terminalPrint();
void setup();
void movieOn();
void movieOff();
void screenShot();
float4 centerOfMass();
float4 linearVelocity();
void zeroOutSystem();

//#include "./callBackFunctions.h"

Body* Bodies = NULL;

void addBody(Body newBody) 
{
    // Reallocate memory to accommodate the new body
    if (NumberOfBodies >= Capacity) //if the new body will exceed the current capacity
    {
        Capacity *= 2; //double the capacity
        Body* temp = (Body*)realloc(Bodies, Capacity*sizeof(Body)); //reallocate memory to accommodate the new body
        if (temp == NULL)  //if memory allocation fails
        {
            fprintf(stderr, "Memory allocation failed\n");
            exit(1);
        }
        Bodies = temp;//assign the new memory to the bodies array, so long as memory allocation was successful
        //printf("Reallocated memory to capacity: %d\n", capacity);
    }

    /// Add the new body to the array
    Bodies[NumberOfBodies] = newBody;

    // Increment the number of bodies
    NumberOfBodies++;

	drawPicture();
    //for debugging
    //printf("Body %d added at (%f, %f, %f) with velocity (%f, %f, %f)\n", newBody.id, newBody.pos.x, newBody.pos.y, newBody.pos.z, newBody.vel.x, newBody.vel.y, newBody.vel.z);
}

void freeBodies() 
{
    free(Bodies);
}

void Display()
{
	glClear(GL_COLOR_BUFFER_BIT);
	glClear(GL_DEPTH_BUFFER_BIT);
	drawPicture();
}

void idle()
{
	nBody();
}

void reshape(int w, int h)
{
	glViewport(0, 0, (GLsizei) w, (GLsizei) h);
}

// void reshape(int w, int h)
// {
//     glViewport(0, 0, (GLsizei) w, (GLsizei) h);
//     glMatrixMode(GL_PROJECTION);
//     glLoadIdentity();

//     // Calculate aspect ratio
//     float aspect = (float)w / (float)h;

//     // Adjust the projection matrix based on the aspect ratio
//     if (aspect >= 1.0f) {
//         // Wider than tall
//         glOrtho(-aspect, aspect, -1.0, 1.0, Near, Far);
//     } else {
//         // Taller than wide
//         glOrtho(-1.0, 1.0, -1.0 / aspect, 1.0 / aspect, Near, Far);
//     }

//     glMatrixMode(GL_MODELVIEW);
//     glLoadIdentity();
// }

void KeyPressed(unsigned char key, int x, int y)
{	
	if(key == 'q')
	{
		// pclose(ffmpeg);
		if (ffmpeg != NULL) {
            pclose(ffmpeg);
            ffmpeg = NULL; // Optionally set to NULL after closing
        } else {
            fprintf(stderr, "Warning: Attempted to close a NULL file pointer\n");
        }
		glutDestroyWindow(Window);
		printf("\nGood Bye\n");
		exit(0);
	}
	if(key == 'o')
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glOrtho(-1.0, 1.0, -1.0, 1.0, Near, Far);
		glMatrixMode(GL_MODELVIEW);
		drawPicture();
	}
	if(key == 'f')
	{
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glFrustum(-0.2, 0.2, -0.2, 0.2, Near, Far);
		glMatrixMode(GL_MODELVIEW);
		drawPicture();
	}
	if(key == 'p')
	{
		if(Pause == 1) Pause = 0;
		else Pause = 1;
		drawPicture();
	}
	if(key == 't') // Turns tracers on and off
	{
		if(Trace == 1) Trace = 0;
		else Trace = 1;
		drawPicture();
	}
	if(key == 'M')  // Movie on/off
	{
		if(MovieFlag == 0) 
		{
			MovieFlag = 1;
			movieOn();
		}
		else 
		{
			MovieFlag = 0;
			movieOff();
		}
	}
	
	if(key == 'S')  // Screenshot
	{	
		screenShot();
	}

	if(key == 'C') // Center out system
	{
		zeroOutSystem();
		drawPicture();
	}

	if(key == 'N') // Turn on or off adding bodies
	{
		if(LClickOn == 1) LClickOn = 0;
		else LClickOn = 1;
	}

	if(key== '0')
	{
		RandomColor = 1;
		NextColorString = "Random";
	}
	if(key== '1')
	{
		RandomColor = 0;
		NextColor={1.0f,0.4f,0.5f,1.0f};
		NextColorString = "Pink";
	}
	if(key== '2')
	{
		RandomColor = 0;
		NextColor={0.9f,1.0f,0.2f,1.0f};
		NextColorString = "Yellow";
	}
	if(key== '3')
	{
		RandomColor = 0;
		NextColor={0.9f,0.07f,0.07f,1.0f};
		NextColorString = "Red";
	}
	if(key== '4')
	{
		RandomColor = 0;
		NextColor={0.9f,0.45f,0.07f,1.0f};
		NextColorString = "Orange";
	}
	if(key== '5')
	{
		RandomColor = 0;
		NextColor={0.5f,0.92,0.4f,1.0f};
		NextColorString = "Green";
	}
	if(key== '6')
	{
		RandomColor = 0;
		NextColor={0.4f,0.69f,0.92f,1.0f};
		NextColorString = "Blue";
	}
	if(key== '7')
	{
		RandomColor = 0;
		NextColor={0.4f,0.4f,0.92f,1.0f};
		NextColorString = "Purple";
	}
	if(key== '8')
	{
		RandomColor = 0;
		NextColor={1.0f,1.0f,1.0f,1.0f};
		NextColorString = "White";
	}
	if(key== '9')
	{
		RandomColor = 0;
		NextColor={0.0f,0.0f,0.0f,1.0f};
		NextColorString = "Black";
	}
	if(key=='B'){
		NumberOfBodies--;
		drawPicture();
	}
	if(key=='R'){
		NumberOfBodies = 0;
		drawPicture();
	}
	if(key=='G'){
		if(GToggle == 0) 
		{
			GToggle = 1;
			HToggle = 0;
		}
		else GToggle = 0;
	}
	if(key=='H'){
		if(HToggle == 0)
		{
			HToggle = 1;
			GToggle = 0;
		}
		else HToggle = 0;
	}
	if(key=='['){
		if(GToggle == 1) G -= 0.1;
		if(HToggle == 1) H -= 0.01;
	}
	if(key==']'){
		if(GToggle == 1) G += 0.1;
		if(HToggle == 1) H += 0.01;
	}

	terminalPrint();
}

void mousePassiveMotionCallback(int x, int y) 
{
	// This function is called when the mouse moves without any button pressed
	// x and y are the current mouse coordinates
	
	// x and y come in as 0 to XWindowSize and 0 to YWindowSize. 
	// Use this if you choose to.
}

// This is called when you push a mouse button.
void mymouse(int button, int state, int x, int y)
{	
	if(state == GLUT_DOWN)
	{	
		if(button == GLUT_LEFT_BUTTON)
		{	
			// Do stuff in here if you choose to when the left mouse button is pressed.
			if(LClickOn == 1)
			{
				// float VelocityMax = 1.0;
				//generate random numbers for all the properties of the new body
                int index = NumberOfBodies; // Define and initialize index
				float xpos = (float)x/(float)XWindowSize*2.0 - 1.0;
				float ypos = -(float)y/(float)YWindowSize*2.0 + 1.0;
                float mass = MassOfBody;

                float colorx = ((float)rand()/(float)RAND_MAX);
                float colory = ((float)rand()/(float)RAND_MAX);
                float colorz = ((float)rand()/(float)RAND_MAX);

                Body newBody; //create a new body with the body struct

                //assign all the properties of the new body
                newBody.id = index;
				if(index == 0)
				{
					newBody.mass = 50.0;
					newBody.pos = {0.0, 0.0, 0.0, 0.0}; // Directly assign values to float4
					newBody.vel = {0.0, 0.0, 0.0, 0.0}; // Directly assign values to float4
					newBody.color = {1.0, 1.0, 1.0, 1.0f}; // Directly assign values to float4
				}
				else
				{
					newBody.mass = mass;
					newBody.pos = {xpos, ypos, 0.0f, 0.0f}; // Directly assign values to float4
					newBody.vel = {4*ypos, -4*xpos, 0.0f, 0.0f}; // Directly assign values to float4
					if(RandomColor != 1)
						newBody.color = NextColor; // Directly assign values to float4
					else
						newBody.color = {(float)rand()/(float)RAND_MAX, (float)rand()/(float)RAND_MAX, (float)rand()/(float)RAND_MAX, 1.0f}; // Directly assign values to float4

				}
                newBody.force = {0.0f, 0.0f, 0.0f, 0.0f}; // Directly assign values to float4

                addBody(newBody);
			}
		}
		else if(button == GLUT_RIGHT_BUTTON) // Right Mouse button down
		{
			// Do stuff in here if you choose to when the right mouse button is pressed.
		}
		else if(button == GLUT_MIDDLE_BUTTON)
		{
			// Do stuff in here if you choose to when the middle mouse button is pressed.
		}
	}
	
	// If no mouse button is down (state 0, they don't have a nice word like GLUT_NOT_DOWN) 
	// but you move the mouse wheel this is called.
	if(state == 0)
	{
		float dz = 0.05f;
		// When you turn the mouse wheel forward this is called.
		if(button == 3)
		{
			glTranslatef(0.0, 0.0, dz);
			drawPicture();
			terminalPrint();
		}
		
		// When you turn the mouse wheel backward this is called.
		else if(button == 4)
		{
			glTranslatef(0.0, 0.0, -dz);
			drawPicture();
			terminalPrint();
		}
	}
}

string getTimeStamp()
{
	// Want to get a time stamp string representing current date/time, so we have a
	// unique name for each video/screenshot taken.
	time_t t = time(0); 
	struct tm * now = localtime( & t );
	int month = now->tm_mon + 1, day = now->tm_mday, year = now->tm_year, 
				curTimeHour = now->tm_hour, curTimeMin = now->tm_min, curTimeSec = now->tm_sec;
	stringstream smonth, sday, syear, stimeHour, stimeMin, stimeSec;
	smonth << month;
	sday << day;
	syear << (year + 1900); // The computer starts counting from the year 1900, so 1900 is year 0. So we fix that.
	stimeHour << curTimeHour;
	stimeMin << curTimeMin;
	stimeSec << curTimeSec;
	string timeStamp;

	if (curTimeMin <= 9)	
		timeStamp = smonth.str() + "-" + sday.str() + "-" + syear.str() + '_' + stimeHour.str() + ".0" + stimeMin.str() + 
					"." + stimeSec.str();
	else			
		timeStamp = smonth.str() + "-" + sday.str() + '-' + syear.str() + "_" + stimeHour.str() + "." + stimeMin.str() +
					"." + stimeSec.str();
	return timeStamp;
}

void movieOn()
{
	string ts = getTimeStamp();
	ts.append(".mp4");

	// Setting up the movie buffer.
	/*const char* cmd = "ffmpeg -loglevel quiet -r 60 -f rawvideo -pix_fmt rgba -s 1000x1000 -i - "
		      "-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip output.mp4";*/

	// string baseCommand = "ffmpeg -loglevel quiet -r 60 -f rawvideo -pix_fmt rgba -s 1000x1000 -i - "
	// 			"-c:v libx264rgb -threads 0 -preset fast -y -pix_fmt yuv420p -crf 0 -vf vflip ";
	string baseCommand = "ffmpeg -loglevel quiet -r 60 -f rawvideo -pix_fmt rgba -s " + to_string(XWindowSize) + "x" + to_string(YWindowSize) + " -i - "
				"-c:v libx264rgb -threads 0 -preset fast -y -pix_fmt yuv420p -crf 0 -vf vflip ";

	string z = baseCommand + ts;

	const char *ccx = z.c_str();
	MovieFile = popen(ccx, "w");
	//Buffer = new int[XWindowSize*YWindowSize];
	Buffer = (int*)malloc(XWindowSize*YWindowSize*sizeof(int));
	MovieOn = 1;
}

void movieOff()
{
	if(MovieOn == 1) 
	{
		pclose(MovieFile);
	}
	free(Buffer);
	MovieOn = 0;
}

void screenShot()
{	
	int pauseFlag;
	FILE* ScreenShotFile;
	int* buffer;

	string tempString = "ffmpeg -loglevel quiet -framerate 60 -f rawvideo -pix_fmt rgba -s " + to_string(XWindowSize) + "x" + to_string(YWindowSize) + " -i - "
				"-c:v libx264rgb -threads 0 -preset fast -y -crf 0 -vf vflip output1.mp4";
	// const char* cmd = "ffmpeg -loglevel quiet -framerate 60 -f rawvideo -pix_fmt rgba -s 1000x1000 -i - "
	// 			"-c:v libx264rgb -threads 0 -preset fast -y -crf 0 -vf vflip output1.mp4";
	//const char* cmd = "ffmpeg -r 60 -f rawvideo -pix_fmt rgba -s 1000x1000 -i - "
	//              "-threads 0 -preset fast -y -pix_fmt yuv420p -crf 21 -vf vflip output1.mp4";
	const char* cmd = tempString.c_str();
	ScreenShotFile = popen(cmd, "w");
	buffer = (int*)malloc(XWindowSize*YWindowSize*sizeof(int));
	
	if(Pause == 0) 
	{
		Pause = 1;
		pauseFlag = 0;
	}
	else
	{
		pauseFlag = 1;
	}
	
	for(int i =0; i < 1; i++)
	{
		drawPicture();
		glReadPixels(5, 5, XWindowSize, YWindowSize, GL_RGBA, GL_UNSIGNED_BYTE, buffer);
		fwrite(buffer, sizeof(int)*XWindowSize*YWindowSize, 1, ScreenShotFile);
	}
	
	pclose(ScreenShotFile);
	free(buffer);

	string ts = getTimeStamp(); // Only storing in a separate variable for debugging purposes.
	string s = "ffmpeg -loglevel quiet -i output1.mp4 -qscale:v 1 -qmin 1 -qmax 1 " + ts + ".jpeg";
	// Convert back to a C-style string.
	const char *ccx = s.c_str();
	system(ccx);
	system("rm output1.mp4");
	printf("\nScreenshot Captured: \n");
	cout << "Saved as " << ts << ".jpeg" << endl;

	
	//system("ffmpeg -i output1.mp4 screenShot.jpeg");
	//system("rm output1.mp4");
	
	Pause = pauseFlag;
	//ffmpeg -i output1.mp4 output_%03d.jpeg
}

void setSimulationParameters()
{
	NumberOfBodies = 0;

	TotalRunTime = 10000.0;

	Dt = 0.002;

	// This is a lennard-Jones type force G*m1*m2/(r^2) - H*m1*m2/(r^4).
	// If you want a gravity type force just set G to your gravity and set H equal 0.
	G = 0.4;

	H = 0.02;

	dForce = 0.9;

	Epsilon = 0.01;

	MassOfBody = 1.0;

	DiameterOfBody = 0.1;

	VelocityMax = 0.0;

	Drag = 0.001;

	DrawRate = 8;
	
	PrintRate = 100;
}

void allocateMemory()
{
// Allocate initial memory for the bodies array
    Bodies = (Body*)malloc(Capacity*sizeof(Body));
    if (Bodies == NULL) 
    {
        fprintf(stderr, "Initial memory allocation failed\n");
        exit(1);
    }
    printf("Initial memory allocated with capacity: %d\n", Capacity);
}

void setInitailConditions()
{
    float dx, dy, dz, d, d2;
    int test;
	time_t t;
	
	srand((unsigned) time(&t));
	for(int i = 0; i < NumberOfBodies; i++)
	{
		test = 0;
		while(test == 0)
		{
			// 2D Box Shape
			// Get random number between -1 at 1.
			Bodies[i].pos.x = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
			Bodies[i].pos.y = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
			Bodies[i].pos.z = ((float)rand()/(float)RAND_MAX)*2.0 - 1.0; //0.0;
			// Heart Shape
			// float temp = 2*PI*((float)rand()/(float)RAND_MAX);
			// BodyPositionX[i] = 16*pow(sin(temp),3)/7.0;
			// BodyPositionY[i] = (13*cos(temp) - 5*cos(2*temp) - 2*cos(3*temp) - cos(4*temp))/7.0;
			// BodyPositionZ[i] = 0.0;
			test = 1;
			
			for(int j = 0; j < i; j++)
			{
				dx = Bodies[i].pos.x - Bodies[j].pos.x;
				dy = Bodies[i].pos.y - Bodies[j].pos.y;
				dz = Bodies[i].pos.z - Bodies[j].pos.z;
				d2  = dx*dx + dy*dy + dz*dz;
				d = sqrt(d2);
				if(d < DiameterOfBody)
				{
					test = 0;
					break;
				}
			}
			
			if(test == 1)
			{
				Bodies[i].vel.x = 0.0; //VelocityMax*((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
				Bodies[i].vel.y = 0.0; //VelocityMax*((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
				Bodies[i].vel.z = 0.0;  //VelocityMax*((float)rand()/(float)RAND_MAX)*2.0 - 1.0;
				
				Bodies[i].color.x = ((float)rand()/(float)RAND_MAX);
				Bodies[i].color.y = ((float)rand()/(float)RAND_MAX);
				Bodies[i].color.z = ((float)rand()/(float)RAND_MAX);
			}
		}
	}
}

float4 centerOfMass()
{
	float totalMass;
	float4 centerOfMass;
	
	centerOfMass.x = 0.0;
	centerOfMass.y = 0.0;
	centerOfMass.z = 0.0;
	totalMass = 0.0;
	
	for(int i = 0; i < NumberOfBodies; i++)
	{
    	centerOfMass.x += Bodies[i].pos.x*Bodies[i].mass;
		centerOfMass.y += Bodies[i].pos.y*Bodies[i].mass;
		centerOfMass.z += Bodies[i].pos.z*Bodies[i].mass;
		totalMass += Bodies[i].mass;
	}
	centerOfMass.x /= totalMass;
	centerOfMass.y /= totalMass;
	centerOfMass.z /= totalMass;
	
	return(centerOfMass);
}

float4 linearVelocity()
{
	float totalMass;
	float4 linearVelocity;
	
	linearVelocity.x = 0.0;
	linearVelocity.y = 0.0;
	linearVelocity.z = 0.0;
	totalMass = 0.0;
	
	for(int i = 0; i < NumberOfBodies; i++)
	{
    	linearVelocity.x += Bodies[i].vel.x*Bodies[i].mass;
		linearVelocity.y += Bodies[i].vel.y*Bodies[i].mass;
		linearVelocity.z += Bodies[i].vel.z*Bodies[i].mass;
		totalMass += Bodies[i].mass;
	}
	linearVelocity.x /= totalMass;
	linearVelocity.y /= totalMass;
	linearVelocity.z /= totalMass;
	
	return(linearVelocity);
}

void zeroOutSystem()
{
	float4 pos, vel;
	pos = centerOfMass();
	vel = linearVelocity();
		
	for(int i = 0; i < NumberOfBodies; i++)
	{
		Bodies[i].pos.x -= pos.x;
		Bodies[i].pos.y -= pos.y;
		Bodies[i].pos.z -= pos.z;
		
		Bodies[i].vel.x -= vel.x;
		Bodies[i].vel.y -= vel.y;
		Bodies[i].vel.z -= vel.z;
	}
}

void drawPicture()
{
	if(Trace == 0)
	{
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);
	}
		
	for(int i = 0; i < NumberOfBodies; i++)
	{
		// BodyColorX[i] = (float)rand()/(float)RAND_MAX;
		// BodyColorY[i] = (float)rand()/(float)RAND_MAX;
		// BodyColorZ[i] = (float)rand()/(float)RAND_MAX;
		glColor3d(Bodies[i].color.x, Bodies[i].color.y, Bodies[i].color.z);
		glPushMatrix();
			glTranslatef(Bodies[i].pos.x, Bodies[i].pos.y, Bodies[i].pos.z);
			glutSolidSphere(DiameterOfBody/2.0, 20, 20);
		glPopMatrix();
	}
	glutSwapBuffers();
	
	if(MovieOn == 1)
	{
		glReadPixels(5, 5, XWindowSize, YWindowSize, GL_RGBA, GL_UNSIGNED_BYTE, Buffer);
		fwrite(Buffer, sizeof(int)*XWindowSize*YWindowSize, 1, MovieFile);
	}
}

void getForces(Body* bodies, float G, float H, float Epsilon, float drag, float dt, int n)
{
	float dx, dy, dz, d2, d;
	float forceMag;
    
	for(int i = 0; i < n; i++)
	{
		bodies[i].force.x = 0.0;
		bodies[i].force.y = 0.0;
		bodies[i].force.z = 0.0;
	}
	
	// Getting force
	for(int i = 0; i < n; i++)
	{   
		for(int j = i + 1; j < n; j++)
		{
			dx = bodies[j].pos.x - bodies[i].pos.x;
			dy = bodies[j].pos.y - bodies[i].pos.y;
			dz = bodies[j].pos.z - bodies[i].pos.z;
		 	d2 = dx*dx + dy*dy + dz*dz + Epsilon;
		 	d = sqrt(d2);
			forceMag  = (G*bodies[i].mass*bodies[j].mass)/(d2) - (H*bodies[i].mass*bodies[j].mass)/(d2*d2);
			bodies[i].force.x += forceMag*dx/d;
			bodies[i].force.y += forceMag*dy/d;
			bodies[i].force.z += forceMag*dz/d;
			bodies[j].force.x -= forceMag*dx/d;
			bodies[j].force.y -= forceMag*dy/d;
			bodies[j].force.z -= forceMag*dz/d;
		}
    }
    
    // Updating positions
	for(int i = 0; i < n; i++)
	{
		bodies[i].vel.x += ((bodies[i].force.x - drag*bodies[i].vel.x)/bodies[i].mass)*dt;
		bodies[i].vel.y += ((bodies[i].force.y - drag*bodies[i].vel.y)/bodies[i].mass)*dt;
		bodies[i].vel.z += ((bodies[i].force.z - drag*bodies[i].vel.z)/bodies[i].mass)*dt;
		
		bodies[i].pos.x += bodies[i].vel.x*dt;
		bodies[i].pos.y += bodies[i].vel.y*dt;
		bodies[i].pos.z += bodies[i].vel.z*dt;
	}
	// Force changes over time
	// G *= (1 - dForce);
	// H *= (1 - dForce);

	// Move the system so that the first body is at the origin with no velocity
	for(int i = 0; i < n; i++)
	{
		bodies[i].pos.x -= bodies[0].pos.x;
		bodies[i].pos.y -= bodies[0].pos.y;
		bodies[i].pos.z -= bodies[0].pos.z;
		bodies[i].vel.x -= bodies[0].vel.x;
		bodies[i].vel.y -= bodies[0].vel.y;
		bodies[i].vel.z -= bodies[0].vel.z;
	}	
}

void nBody()
{
	if(Pause != 1)
	{	
		getForces(Bodies,  G,  H,  Epsilon,  Drag, Dt, NumberOfBodies);
        
        DrawTimer++;
		if(DrawTimer == DrawRate) 
		{
			drawPicture();
			DrawTimer = 0;
		}
		
		PrintTimer++;
		if(PrintTimer == PrintRate) 
		{
			terminalPrint();
			PrintTimer = 0;
		}
		
		RunTime += Dt; 
		if(TotalRunTime < RunTime)
		{
			printf("\n\n Done\n");
			exit(0);
		}
	}
}

void terminalPrint()
{
	/*
	default  \033[0m
	Black:   \033[0;30m
	Red:     \033[0;31m
	Green:   \033[0;32m
	Yellow:  \033[0;33m
	Blue:    \033[0;34m
	Magenta: \033[0;35m
	Cyan:    \033[0;36m
	White:   \033[0;37m
	printf("\033[0;30mThis text is black.\033[0m\n");
	
	BOLD_ON  "\e[1m"
	BOLD_OFF   "\e[m"
	*/
	
	system("clear");
	
	printf("\n\n");
	printf("\033[0m");
	printf("\n");
	printf(" o/f: Orthographic/Frustum Projection");
	printf("\n");
	printf("\n Mouse Wheel: Zoom in/out\n");
	printf("\n [] : Increase/Decrease G or H\n");
	printf("\n G: Edit G toggle --> ");
	printf(" G ");
	if (GToggle == 1) 
	{
		printf("\e[1m" " \033[0;32mCan be edited\n" "\e[m");
	}
	else 
	{
		printf("\e[1m" " \033[0;31mCannot be edited\n" "\e[m");
	}
	printf("\n G: %f", G);
	printf("\n");
	printf("\n H: Edit H toggle --> ");
	printf(" H ");
	if (HToggle == 1) 
	{
		printf("\e[1m" " \033[0;32mCan be edited\n" "\e[m");
	}
	else 
	{
		printf("\e[1m" " \033[0;31mCannot be edited\n" "\e[m");
	}
	printf("\n H: %f", H);
	printf("\n");
	printf("\n p: Pause on/off toggle --> ");
	printf(" The simulation is:");
	if (Pause == 1) 
	{
		printf("\e[1m" " \033[0;31mPaused\n" "\e[m");
	}
	else 
	{
		printf("\e[1m" " \033[0;32mRunning\n" "\e[m");
	}
	
	printf("\n");
	printf("\033[0m");
	printf(" t: Trace on/off toggle --> ");
	printf(" Trace is:");
	if (Trace == 1) 
	{
		printf("\e[1m" " \033[0;31mOn\n" "\e[m");
	}
	else 
	{
		printf("\e[1m" " \033[0;32mOff\n" "\e[m");
	}
	
	printf("\n M: Video On/Off toggle --> ");
	if (MovieFlag == 0) 
	{
		printf("\033[0;31m");
		printf(BOLD_ON "Video Recording Off" BOLD_OFF); 
	}
	else 
	{
		printf("\033[0;32m");
		printf(BOLD_ON "Video Recording On" BOLD_OFF);
	}

	printf("\n");
	printf("\n N: Add Bodies On/Off toggle --> ");
	if (LClickOn == 0) 
	{
		printf("\033[0;31m");
		printf(BOLD_ON "Adding Bodies Off" BOLD_OFF); 
	}
	else 
	{
		printf("\033[0;32m");
		printf(BOLD_ON "Adding Bodies On" BOLD_OFF);
	}
	
	printf("\n");
	printf("\n S: Screenshot");

	printf("\n");
	printf("\n C: Center out system");
	
	printf("\n");
	printf("\n q: Terminates the simulation");

	printf("\n");
	printf("\n 0: Random");
	printf("\n 1: Pink!!");
	printf("\n 2: Yellow");
	printf("\n 3: Red");
	printf("\n 4: Orange");
	printf("\n 5: Green");
	printf("\n 6: Blue");
	printf("\n 7: Purple");
	printf("\n 8: White");
	printf("\n 9: Black");
	printf("\n Next Color: %s", NextColorString.c_str());
	
	printf("\n");
	printf("\n B: Remove Last Body");
	
	printf("\n");
	printf("\n R: Remove All Bodies");
	printf("\n");
}

void setup()
{	
	setSimulationParameters();
	allocateMemory();
	setInitailConditions();
	zeroOutSystem();
    	DrawTimer = 0;
    	PrintRate = 0;
	RunTime = 0.0;
	Trace = 0;
	Pause = 1;
	MovieOn = 0;
	terminalPrint();
}

int main(int argc, char** argv)
{
	setup();
	
	XWindowSize = 1000;
	YWindowSize = 1000; 
	//Buffer = new int[XWindowSize*YWindowSize];

	// Clip plains
	Near = 0.2;
	Far = 30.0;

	//Direction here your eye is located location
	EyeX = 0.0;
	EyeY = 0.0;
	EyeZ = 2.0;

	//Where you are looking
	CenterX = 0.0;
	CenterY = 0.0;
	CenterZ = 0.0;

	//Up vector for viewing
	UpX = 0.0;
	UpY = 1.0;
	UpZ = 0.0;
	
	glutInit(&argc,argv);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutInitWindowSize(XWindowSize,YWindowSize);
	glutInitWindowPosition(5,5);
	Window = glutCreateWindow("N Body");
	
	gluLookAt(EyeX, EyeY, EyeZ, CenterX, CenterY, CenterZ, UpX, UpY, UpZ);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	glFrustum(-0.2, 0.2, -0.2, 0.2, Near, Far);
	glMatrixMode(GL_MODELVIEW);
	glClearColor(0.0, 0.0, 0.0, 0.0);
	
	GLfloat light_position[] = {1.0, 1.0, 1.0, 0.0};
	GLfloat light_ambient[]  = {0.0, 0.0, 0.0, 1.0};
	GLfloat light_diffuse[]  = {1.0, 1.0, 1.0, 1.0};
	GLfloat light_specular[] = {1.0, 1.0, 1.0, 1.0};
	GLfloat lmodel_ambient[] = {0.2, 0.2, 0.2, 1.0};
	GLfloat mat_specular[]   = {1.0, 1.0, 1.0, 1.0};
	GLfloat mat_shininess[]  = {10.0};
	glShadeModel(GL_SMOOTH);
	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	glLightfv(GL_LIGHT0, GL_POSITION, light_position);
	glLightfv(GL_LIGHT0, GL_AMBIENT, light_ambient);
	glLightfv(GL_LIGHT0, GL_DIFFUSE, light_diffuse);
	glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mat_specular);
	glMaterialfv(GL_FRONT, GL_SHININESS, mat_shininess);
	glEnable(GL_LIGHTING);
	glEnable(GL_LIGHT0);
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
	
	glutPassiveMotionFunc(mousePassiveMotionCallback);
	glutMouseFunc(mymouse);
	glutDisplayFunc(Display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(KeyPressed);
	glutIdleFunc(idle);
	terminalPrint();
	glutMainLoop();
	return 0;
}