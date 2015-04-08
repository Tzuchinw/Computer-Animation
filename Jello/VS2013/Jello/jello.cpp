/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1
  Tzu-Chin Wang (Rin)

*/

#include "jello.h"
#include "showCube.h"
#include "input.h"
#include "physics.h"
#include "performanceCounter.h"

// camera parameters
double Theta = pi / 6;
double Phi = pi / 6;
double R = 6;

// mouse control
int g_iMenuId;
int g_vMousePos[2];
int g_iLeftMouseButton,g_iMiddleMouseButton,g_iRightMouseButton;

// number of images saved to disk so far
int sprite=0;

// these variables control what is displayed on screen
int shear=0, bend=0, structural=1, pause=0, viewingMode=0, saveScreenToFile=0;

struct world jello;
struct ball sphere;

int windowWidth, windowHeight;

//Star Wars text scrolling speed
GLdouble textDownScroll = -10.0;

//Collision method switcher
bool positionBasedCollision = false;
bool switchCollision = false;

//FPS counter
PerformanceCounter theTimer;
double timePerFrame;

//	Star Wars style credit
void showCredit()
{
	int lenghOfQuote, i;
	
	int numberOfQuotes = 0;
	char quote[3][25];

	strcpy(quote[2],"Computer Animation");
    strcpy(quote[1],"Assignment 1: the jello");
    strcpy(quote[0],"Tzu-Chin Wang (Rin)");

    numberOfQuotes = 3;

	glPushMatrix();

    glTranslated(0.0, -100, textDownScroll);
	glRotated(180, 0.0, 1.0, 0.0);
	glRotated(-120, 1.0, 0.0, 0.0);
    glScalef(0.1, 0.1, 0.1);

    for(int i=0; i<numberOfQuotes; i++)
    {
		lenghOfQuote = (int)strlen(quote[i]);
		glPushMatrix();
		glTranslated(-(lenghOfQuote * 35), (i*200), 0.0);
        for (int j = 0; j < lenghOfQuote; j++)
        {
			glColor3f(1.0,1.0,1.0);
			glutStrokeCharacter(GLUT_STROKE_ROMAN, quote[i][j]);
		}
		glPopMatrix();
	}
	glPopMatrix();
	glFlush();
}

//	Show FPS and if it's on position based collision mode
void showState()
{
	int lenghOfState, i;
	
	int numberOfStates = 0;
	char State[5][50];
	
	//FPS
	int fps = 1/timePerFrame;
	sprintf (State[0], "FPS: %d", fps);
	//Collision method
	sprintf (State[1], "Position Based Collision: %s (a)", positionBasedCollision ? "true" : "false");

    numberOfStates = 2;

	glPushMatrix();
	glTranslated(1.75, -2.0, 1.75);
	glRotated(90, 1.0, 0.0, 0.0);
	glRotated(180, 0.0, 1.0, 0.0);
	glScalef(0.0015, 0.0015, 0.0015);
	
	for(int i=0; i<numberOfStates; i++)
    {
		int lenghOfState = (int)strlen(State[i]);
		glTranslated(0, -(i*150), 0.0);
		glPushMatrix();
		for (int j = 0; j < lenghOfState; j++)
		{
			glColor3f(1.0,1.0,0.0);
			glutStrokeCharacter(GLUT_STROKE_ROMAN, State[i][j]);
		}
		glPopMatrix();
	}
	glPopMatrix();
	
	//glWindowPos(10, 20);
	//glutBitmapString(GLUT.BITMAP_HELVETICA_12, WHT);
}

//show the ball
void showBall(struct ball * sphere)
{
	glPushMatrix();
	glColor3f(0.8,0,0.75);
	GLUquadric *quad;
	quad = gluNewQuadric();
	glTranslatef(sphere->p.x,sphere->p.y,sphere->p.z);
	gluSphere(quad, sphere->r, 10,10);
	glPopMatrix();
}

//moving the ball
void moveBall(struct ball * sphere, const double speed, const double &time)
{
	static bool dirX = true;

	if(sphere->p.x > 3)	dirX =false;
	else if(sphere->p.x < -3)	dirX =true;

	if(dirX)	sphere->p.x += speed * time;
	else	sphere->p.x -= speed * time;

}

//check if point already in the list
int check(point * target, int n, point list)
{
	for(int i = 0; i<n; i++)
	{
		if(target[i]== list)	return n;
	}
	return -1;
}

//show the inclined plane, which I hate the most
void showInclinedPlane(struct world * jello)
{
	point intersect[12];
	point todraw[12];//not yet sorted!

	double a = jello->a;
	double b = jello->b;
	double c = jello->c;
	double d = jello->d;

	//find all intersect points
	intersect[0].y = 2;
	intersect[0].z = 2;
	intersect[0].x = -((b * intersect[0].y + c * intersect[0].z + d) / a);

	intersect[1].y = 2;
	intersect[1].z = -2;
	intersect[1].x = -((b * intersect[1].y + c * intersect[1].z + d) / a);

	intersect[2].y = -2;
	intersect[2].z = -2;
	intersect[2].x = -((b * intersect[2].y + c * intersect[2].z + d) / a);

	intersect[3].y = -2;
	intersect[3].z = 2;
	intersect[3].x = -((b * intersect[3].y + c * intersect[3].z + d) / a);
	
	intersect[4].x = 2;
	intersect[4].z = 2;
	intersect[4].y = -((a * intersect[4].x + c * intersect[4].z + d) / b);

	intersect[5].x = -2;
	intersect[5].z = 2;
	intersect[5].y = -((a * intersect[5].x + c * intersect[5].z + d) / b);

	intersect[6].x = -2;
	intersect[6].z = -2;
	intersect[6].y = -((a * intersect[6].x + c * intersect[6].z + d) / b);

	intersect[7].x = 2;
	intersect[7].z = -2;
	intersect[7].y = -((a * intersect[7].x + c * intersect[7].z + d) / b);

	intersect[8].x = 2;
	intersect[8].y = 2;
	intersect[8].z = -((a * intersect[8].x + b * intersect[8].y + d) / c);

	intersect[9].x = -2;
	intersect[9].y = 2;
	intersect[9].z = -((a * intersect[9].x + b * intersect[9].y + d) / c);

	intersect[10].x = -2;
	intersect[10].y = -2;
	intersect[10].z = -((a * intersect[10].x + b * intersect[10].y + d) / c);

	intersect[11].x = 2;
	intersect[11].y = -2;
	intersect[11].z = -((a * intersect[11].x + b * intersect[11].y + d) / c);
	
	int i = 0;
	int counter = 0;
	for(i; i<4; i++)
	{
		if(intersect[i].x >= -2 && intersect[i].x <= 2)
		{
			if(check(todraw, counter, intersect[i])==-1)
			{
				todraw[counter] = intersect[i];
				counter++;
			}
		}
	}
	for(i; i<8; i++)
	{
		if(intersect[i].y >= -2 && intersect[i].y <= 2)
		{
			if(check(todraw, counter, intersect[i])==-1)
			{
				todraw[counter] = intersect[i];
				counter++;
			}
		}
	}
	for(i; i<12; i++)
	{
		if(intersect[i].z >= -2 && intersect[i].z <= 2)
		{
			if(check(todraw, counter, intersect[i])==-1)
			{
				todraw[counter] = intersect[i];
				counter++;
			}
		}
	}

	glColor4d(1.0, 1.0, 1.0, 0.1);//alpha doesn't work

	if(counter == 3)//3 points no need to sort points
	{
		glBegin(GL_TRIANGLES);
		glVertex3d(todraw[0].x,todraw[0].y,todraw[0].z);
		glVertex3d(todraw[1].x,todraw[1].y,todraw[1].z);
		glVertex3d(todraw[2].x,todraw[2].y,todraw[2].z);
		glEnd();
		
		glBegin(GL_TRIANGLES);
		glVertex3d(todraw[2].x,todraw[2].y,todraw[2].z);
		glVertex3d(todraw[1].x,todraw[1].y,todraw[1].z);
		glVertex3d(todraw[0].x,todraw[0].y,todraw[0].z);
		glEnd();
		//two side
	}
	else
	{
		//sort point in order!!!!!!
		point new_z(a,b,c);
		new_z.normalize();

		point c;//find center point
		for(int i = 0; i < counter; i++)	c +=todraw[i];
		c/= counter;
		//find new basis axis
		point new_y = (todraw[0] - c);
		new_y.normalize();
		point new_x = new_y.cross(new_z);

		double angle[12];//calculate arc tan
		for(int i = 0; i < counter; i++)
		{
			double xv = todraw[i].dot(new_x);
			double yv = todraw[i].dot(new_y);
			angle[i] = atan2(xv,yv);
		}

		point sortedraw[12];
		for(int i = 0; i < counter; i++)
		{
			double min = angle[0] ;
			int loc = 0;
			for(int j = 0; j < counter; j++)
			{
				if(angle[j] < min)
				{
					min = angle[j];
					loc = j;
				}
			}
			sortedraw[i] = todraw[loc];
			angle[loc] = 99999;//set as used
		}//sorting

		//Draw
		glBegin(GL_POLYGON);
		for(int i = 0; i < counter; i ++)	glVertex3d(sortedraw[i].x,sortedraw[i].y,sortedraw[i].z);
		glEnd();
		//two side
		glBegin(GL_POLYGON);
		for(int i = counter-1; i >= 0; i --)	glVertex3d(sortedraw[i].x,sortedraw[i].y,sortedraw[i].z);
		glEnd();
	}
	
}


void myinit()
{
	sphere.p = point(-2,0,0);
	sphere.r = 0.3;

	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(90.0,1.0,0.01,1000.0);

	// set background color to grey, no! dark green!
	glClearColor(0.05, 0.15, 0.05, 0.0);

	glCullFace(GL_BACK);
	glEnable(GL_CULL_FACE);

	glShadeModel(GL_SMOOTH);
	glEnable(GL_POLYGON_SMOOTH);
	glEnable(GL_LINE_SMOOTH);
	glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);	//doesn't work
	glEnable (GL_BLEND);
	
  
	theTimer.StartCounter();//initialize counter
	return; 
}



void reshape(int w, int h) 
{
	// Prevent a divide by zero, when h is zero.
	// You can't make a window of zero height.
	if(h == 0)
		h = 1;

	glViewport(0, 0, w, h);

	// Reset the coordinate system before modifying
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	// Set the perspective
	double aspectRatio = 1.0 * w / h;
	gluPerspective(60.0f, aspectRatio, 0.01f, 1000.0f);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity(); 

	windowWidth = w;
	windowHeight = h;

	glutPostRedisplay();
}



void display()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();

	// camera parameters are Phi, Theta, R
	gluLookAt(R * cos(Phi) * cos (Theta), R * sin(Phi) * cos (Theta), R * sin (Theta),
				0.0,0.0,0.0, 0.0,0.0,1.0);


	/* Lighting */
	/* You are encouraged to change lighting parameters or make improvements/modifications
		to the lighting model . 
		This way, you will personalize your assignment and your assignment will stick out. 
	*/

	// global ambient light
	GLfloat aGa[] = { 0.1, 0.05, 0.1, 0.0 };	//purple!
  
	// light 's ambient, diffuse, specular
	GLfloat lKa0[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd0[] = { 0.0, 0.0, 1.0, 1.0 };
	GLfloat lKs0[] = { 0.0, 0.0, 1.0, 1.0 };

	GLfloat lKa1[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd1[] = { 0.5, 0.0, 0.0, 1.0 };
	GLfloat lKs1[] = { 0.5, 0.0, 0.0, 1.0 };

	GLfloat lKa2[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd2[] = { 0.5, 0.5, 0.0, 1.0 };
	GLfloat lKs2[] = { 0.5, 0.5, 0.0, 1.0 };

	GLfloat lKa3[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd3[] = { 0.0, 0.5, 0.5, 1.0 };
	GLfloat lKs3[] = { 0.0, 0.5, 0.5, 1.0 };

	GLfloat lKa4[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd4[] = { 0.0, 0.0, 0.5, 1.0 };
	GLfloat lKs4[] = { 0.0, 0.0, 0.5, 1.0 };

	GLfloat lKa5[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd5[] = { 0.5, 0.0, 0.5, 1.0 };
	GLfloat lKs5[] = { 0.5, 0.0, 0.5, 1.0 };

	GLfloat lKa6[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd6[] = { 1.0, 0.0, 0.0, 1.0 };
	GLfloat lKs6[] = { 1.0, 0.0, 0.0, 1.0 };

	GLfloat lKa7[] = { 0.0, 0.0, 0.0, 1.0 };
	GLfloat lKd7[] = { 0.0, 0.5, 0.5, 1.0 };
	GLfloat lKs7[] = { 0.0, 0.5, 0.5, 1.0 };

	// light positions and directions
	static double lp0x = -1.999; 
	static double lp0y = -1.999;
	static double lp0z = -1.999;
	static bool inverted = false;
	double lightmove = (inverted)? -timePerFrame*0.1: timePerFrame*0.1;

	lp0x += timePerFrame*0.2;
	lp0y += timePerFrame*0.2;
	if(lp0x > 2.0)	inverted = true;
	else if(lp0x < -2.0)	inverted = false;

	GLfloat lP0[] = { lp0x, lp0y, 2, 1.0 };
	GLfloat lP1[] = { 1.999, -1.999, -1.999, 1.0 };
	GLfloat lP2[] = { 1.999, 1.999, -1.999, 1.0 };
	GLfloat lP3[] = { -1.999, 1.999, -1.999, 1.0 };
	GLfloat lP4[] = { -1.999, -1.999, 1.999, 1.0 };
	GLfloat lP5[] = { 1.999, -1.999, 1.999, 1.0 };
	GLfloat lP6[] = { 2- lp0x, 2 - lp0y, -2, 1.0 };
	GLfloat lP7[] = { -1.999, 1.999, 1.999, 1.0 };
  
	// jelly material color
	GLfloat mKa[] = { 0.1, 0.1, 0.1, 1.0 };
	GLfloat mKd[] = { 0.75, 0.75, 0.75, 1.0 };
	GLfloat mKs[] = { 1.0, 1.0, 1.0, 1.0 };
	GLfloat mKe[] = { 0.1, 0.0, 0.1, 1.0 };

	/* set up lighting */
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, aGa);
	glLightModelf(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);

	// set up cube color
	glMaterialfv(GL_FRONT, GL_AMBIENT, mKa);
	glMaterialfv(GL_FRONT, GL_DIFFUSE, mKd);
	glMaterialfv(GL_FRONT, GL_SPECULAR, mKs);
	glMaterialfv(GL_FRONT, GL_EMISSION, mKe);
	glMaterialf(GL_FRONT, GL_SHININESS, 120);
    
	// macro to set up light i
	#define LIGHTSETUP(i)\
	glLightfv(GL_LIGHT##i, GL_POSITION, lP##i);\
	glLightfv(GL_LIGHT##i, GL_AMBIENT, lKa##i);\
	glLightfv(GL_LIGHT##i, GL_DIFFUSE, lKd##i);\
	glLightfv(GL_LIGHT##i, GL_SPECULAR, lKs##i);\
	glEnable(GL_LIGHT##i)
  
	LIGHTSETUP (0);
	LIGHTSETUP (1);
	LIGHTSETUP (2);
	LIGHTSETUP (3);
	LIGHTSETUP (4);
	LIGHTSETUP (5);
	LIGHTSETUP (6);
	LIGHTSETUP (7);

	// enable lighting
	glEnable(GL_LIGHTING);    
	glEnable(GL_DEPTH_TEST);

	// show the cube
	showCube(&jello);
	

	showBall(&sphere);
	
	// show the bounding box
	showBoundingBox();
	showInclinedPlane(&jello);
	glDisable(GL_LIGHTING);//make the bounding box shinning

	//Show text
	showCredit();
	showState();
	
 
	glutSwapBuffers();
}

void doIdle()
{
	theTimer.StopCounter();
	timePerFrame = theTimer.GetElapsedTime();
	theTimer.StartCounter();

	char s[20]="picxxxx.ppm";

	// save screen to file
	s[3] = 48 + (sprite / 1000);
	s[4] = 48 + (sprite % 1000) / 100;
	s[5] = 48 + (sprite % 100 ) / 10;
	s[6] = 48 + sprite % 10;

	static int prtscreencounter = 0;
	if (saveScreenToFile == 1)
	{
		if(prtscreencounter==5)
		{
			saveScreenshot(windowWidth, windowHeight, s);
			//saveScreenToFile=0; // save only once, change this if you want continuos image generation (i.e. animation)
			sprite++;
			prtscreencounter = 0;
		}
		prtscreencounter++;
	}

	if (sprite >= 300) // allow only 300 snapshots
	{
		exit(0);	
	}

	if (pause == 0)
	{
    // insert code which appropriately performs one step of the cube simulation:
		for (int i=1; i<=jello.n; i++)
		{
			if (jello.integrator[0] == 'E')		Euler(&jello);
			else								RK4(&jello);

				for(int i = 0; i < 8 ; i++)
					for(int j = 0; j < 8 ; j++)
						for(int k = 0; k < 8 ; k++)
						{
							if(positionBasedCollision)
								collision(&jello, i, j, k);
							collisionBall(&jello, i, j, k, &sphere);
							collisionPlane(&jello, i, j, k);
						}

		}
		if(switchCollision)	//switching collision mode
		{	
			positionBasedCollision = !positionBasedCollision;
			switchCollision = false;
		}
	}

	if (textDownScroll > -100)
		textDownScroll -= 10*timePerFrame;
	else
		textDownScroll = 60;

	double ballspeed = 0.25;
	moveBall(&sphere,ballspeed,timePerFrame);

	glutPostRedisplay();
}

int main (int argc, char ** argv)
{
  if (argc<2)
  {  
    printf ("Oops! You didn't say the jello world file!\n");
    printf ("Usage: %s [worldfile]\n", argv[0]);
    exit(0);
  }

  readWorld(argv[1],&jello);

  glutInit(&argc,argv);
  
  /* double buffered window, use depth testing, 640x480 */
  glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
  
  windowWidth = 640;
  windowHeight = 480;
  glutInitWindowSize (windowWidth, windowHeight);
  glutInitWindowPosition (0,0);
  glutCreateWindow ("Jello cube");

  /* tells glut to use a particular display function to redraw */
  glutDisplayFunc(display);

  /* replace with any animate code */
  glutIdleFunc(doIdle);

  /* callback for mouse drags */
  glutMotionFunc(mouseMotionDrag);

  /* callback for window size changes */
  glutReshapeFunc(reshape);

  /* callback for mouse movement */
  glutPassiveMotionFunc(mouseMotion);

  /* callback for mouse button changes */
  glutMouseFunc(mouseButton);

  /* register for keyboard events */
  glutKeyboardFunc(keyboardFunc);

  /* do initialization */
  myinit();

  /* forever sink in the black hole */
  glutMainLoop();

  return(0);
}

