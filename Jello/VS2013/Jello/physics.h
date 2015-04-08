/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1
  Tzu-Chin Wang (Rin)

*/

#ifndef _PHYSICS_H_
#define _PHYSICS_H_

void computeAcceleration(struct world * jello, struct point a[8][8][8]);

// perform one step of Euler and Runge-Kutta-4th-order integrators
// updates the jello structure accordingly
void Euler(struct world * jello);
void RK4(struct world * jello);
void collision(struct world * jello, int &i, int &j, int &k);
void collisionBall(struct world * jello, int &i, int &j, int &k, struct ball * sphere);
void collisionPlane(struct world * jello, int &i, int &j, int &k);
#endif

