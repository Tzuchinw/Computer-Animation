
Readme 2/12
-----------------------
Computer Animation Assignment 1: the jello THE UPDATE VERSION
Tzu-Chin Wang(Rin)

Dev Environment: Microsoft Visual Studio 2010 on Windows 7 x64
Animation: the screen shot is made once per 5 frame

For the basic requirement, I think what I did is quite normal, iterate all mass point and calculate corresponding forces.
I didn't use those define function given in the starter code, because first I like to write my own code, and there are a few bugs with them(pMAKE and pDIFFERENCE), instead I extended the point class and overloaded serveral operators.
While I'm not familiar with OpenGL (really had bad time trying to figure out how to project things to screen space), I did try to implement some extra features.

Extra features
Moving lights:
There are two light sources moving around.

States Display:
FPS and collision mode is showed on the back of bounding cube

Star Wars Styled Credit:
I love Star Wars!
My neighbors were partying while I was working on this assignment. Their noise made me upset, so that I decided to do some funny stuff to comfort myself…
And surprisingly they played Star Wars title music when I finished implementing this feature, which was kind of creepy…
The Credit is showed right behind the bounding cube, saying "Computer Animation ..."

Position Based Collision:
Position Based Collision is implemented as another collision choice.
User may press 'a' at anytime to switch between two methods.

Poke force:
Left click mouse to apply unified force toward the screen to the jello.

Moving Ball:
Collision(position based) with a moving ball.

Collision with inclined plane:
Collision with inclined plane, which is correctly displayed in the screen.
Note that I did sort the intersect points.

Thank you!


