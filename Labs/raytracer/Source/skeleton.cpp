#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>
#include <limits.h>
#include <math.h> //M_PI

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;

SDL_Event event;

#define SCREEN_WIDTH 100
#define SCREEN_HEIGHT 100
#define FULLSCREEN_MODE false

/* STRUCTURES */

struct Intersection
{
  vec4 position;
  float distance;
  int triangleIndex;
};
struct Camera{
  vec4 position;
  mat4 basis;

};
struct Light{
    vec4 position;
    vec3 color;
};

vector<Triangle> triangles;
/* ASSIGNING STRUCTURES */
Camera camera = {
    .position = vec4(0,0,-2, 1), // AS we want our image plane at 0,0,-1 (because the cornel box is 2x2 with centre 0,0,0) camera is moved 1 unit back
    .basis = mat4(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1))
};
Light dLight{
  .position=  vec4( 0, -0.5, -0.7, 1.0 ),
  .color=14.f * vec3( 1, 1, 1 )
};
/* ASSIGNING VARIABLES/CONSTANTS */
float m = std::numeric_limits<float>::max();
int focalLength = SCREEN_WIDTH/2;
float yaw = 0;
mat4 translation;
vec3 indirectLight = 0.5f*vec3( 1, 1, 1 );




/*------------------*/

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */
bool closestIntersection(vec4 start, vec4 d, const vector<Triangle>& triangles, Intersection& closestIntersectionm, Camera camera);
bool Update();
mat4 RotMatrix(vec3 thetas);
mat4 LookAt(vec3 from,vec3 to);
vec3 DirectLight( const Intersection& i,Light dLight, const vector<Triangle>& triangles );
void Draw(screen* screen);


int main( int argc, char* argv[] )
{

  screen *screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  LoadTestModel(triangles);
  while ( Update())
    {
      Draw(screen);
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
  /* Clear buffer */
  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));

  vec3 black(0.0,0.0,0.0);
  int W=screen->width;
  int H=screen->height;


  for(int i=0; i<W;i++){
    for(int j=0; j<H;j++){
      vec4 d(i - W/2,j-H/2,focalLength,1);//from centre of image to pixels
      Intersection intersect;
      if (closestIntersection(camera.position, d, triangles, intersect,camera)){
        vec3 ro=triangles[intersect.triangleIndex].color;
        vec3 D=DirectLight(intersect,dLight,triangles);
        vec3 total_light= D+indirectLight;
        vec3 light_color= ro*total_light;

        PutPixelSDL(screen, i, j, light_color);

        //vec3 only_direct= ro*D;
        //PutPixelSDL(screen, i, j, only_direct);
        //PutPixelSDL(screen, i, j, ro);


      }
      else{
        PutPixelSDL(screen, i, j, black);
      }
    }
  }
}

/*Place updates of parameters here*/
bool Update()
{
  /*static int t = SDL_GetTicks();
  Compute frame time */
  /*int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;*/
 mat4 translation(vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(0,0,0,1));

 /* TO REMEMBER:
 each vec4 is column.
 trans matrix is vec4(1,0,0,0), vec4(0,1,0,0), vec4(0,0,1,0), vec4(DX,DY,DZ,1)
 meaning: (1,0,0,DX
           0,1,0,DY
           0,0,1,DZ
           0,0,0,1)
*/
 //camera.basis=LookAt(vec3(camera.position[0],camera.position[1],camera.position[2]),vec3(0,0,0));
  SDL_Event e;
  while(SDL_PollEvent(&e))
    {
      if (e.type == SDL_QUIT)
	{
	  return false;
	}
      else
	if (e.type == SDL_KEYDOWN)
	  {
	    int key_code = e.key.keysym.sym;
	    switch(key_code)
	      {

        case SDLK_UP:
        /* zoom in */
        translation[3][2] = 0.1;
        camera.position = translation*camera.position;

    break;
        case SDLK_DOWN:
        /* zoom out */
        translation[3][2] =-0.1;
        camera.position = translation*camera.position;
    break;
        case SDLK_LEFT:
        yaw+=10*M_PI/180;
        camera.basis=RotMatrix(vec3(0, yaw, 0)) * mat4(vec4(1,0,0,(camera.position)[0]),vec4(0,1,0,camera.position[1]),vec4(0,0,1,(camera.position)[2]),vec4(0,0,0,1));
        //dLight.position= glm::inverse(camera.basis)*dLight.position;

    /* Move camera left */
    break;
        case SDLK_RIGHT:
        yaw-=10*M_PI/180;
        camera.basis=RotMatrix(vec3(0, yaw, 0)) * mat4(vec4(1,0,0,(camera.position)[0]),vec4(0,1,0,camera.position[1]),vec4(0,0,1,(camera.position)[2]),vec4(0,0,0,1));
        //dLight.position= glm::inverse(camera.basis)*dLight.position;
    /* Move camera right */
		break;
    case SDLK_w:
        /* zoom in */
        translation[3][2] = 0.1;
        dLight.position = translation*dLight.position;

    break;
        case SDLK_s:
        /* zoom out */
        translation[3][2] =-0.1;
        dLight.position = translation*dLight.position;
    break;
        case SDLK_a:
        translation[3][0] =-0.1;
        dLight.position = translation*dLight.position;
    break;
        case SDLK_d:
        translation[3][0] =0.1;
        dLight.position = translation*dLight.position;
    break;
    case SDLK_q:
        translation[3][1] =0.1;
        dLight.position = translation*dLight.position;
    break;
    case SDLK_e:
        translation[3][1] =-0.1;
        dLight.position = translation*dLight.position;
    break;
	      case SDLK_ESCAPE:
		/* Move camera quit */
		return false;
	      }
	  }
    }
  return true;
}
bool closestIntersection(vec4 start, vec4 d,const vector<Triangle>& triangles,Intersection& closestIntersection, Camera camera){
  closestIntersection.distance = m;
  for( uint i=0; i<triangles.size(); i++ ){
    Triangle triangle = triangles[i];
    vec4 v0 = glm::inverse(camera.basis)*triangle.v0;
    vec4 v1 = glm::inverse(camera.basis)*triangle.v1;
    vec4 v2 = glm::inverse(camera.basis)*triangle.v2;
    vec3 e1 = vec3(v1[0]-v0[0],v1[1]-v0[1],v1[2]-v0[2]);
    vec3 e2 = vec3(v2[0]-v0[0],v2[1]-v0[1],v2[2]-v0[2]);
    vec3 b = vec3(start[0]-v0[0],start[1]-v0[1],start[2]-v0[2]);
    vec3 dir(d[0], d[1], d[2]);
    mat3 A( -dir, e1, e2 );
    vec3 x = glm::inverse( A ) * b;
    if(x[0]>0 && x[1]>=0 && x[2]>=0 && x[1]+x[2]<=1 && x[0] < closestIntersection.distance){
      closestIntersection.distance= x[0];
      closestIntersection.position= start + x[0]*d;
      closestIntersection.triangleIndex = i;
    }
  }
  if (closestIntersection.distance == m) return false;
  return true;
}
/* creates the rotation matrix given a vec3 of angles(thetax,thetay,thetaz)*/
mat4 RotMatrix(vec3 thetas){
  vec4 col1 = vec4(cos(thetas[1])*cos(thetas[2]), cos(thetas[1])*sin(thetas[2]), -sin(thetas[1]), 0);
  vec4 col2 = vec4(-cos(thetas[0])*sin(thetas[2]) + sin(thetas[0]*sin(thetas[1])*cos(thetas[2])), cos(thetas[0])*cos(thetas[2])+sin(thetas[0])*sin(thetas[1])*sin(thetas[2]), sin(thetas[0])*cos(thetas[1]), 0);
  vec4 col3 = vec4(sin(thetas[0])*sin(thetas[2])+cos(thetas[0])*sin(thetas[1])*cos(thetas[2]), -sin(thetas[0])*cos(thetas[2])+cos(thetas[0])*sin(thetas[1])*sin(thetas[2]), cos(thetas[0])*cos(thetas[1]), 0);
  vec4 col4 = vec4(0,0,0,1);
  return (mat4(col1,col2,col3,col4));
}

mat4 LookAt(vec3 from,vec3 to){
  vec3 forward = normalize(from - to);// Z axis
  vec3 tmp(0, 1, 0); // we cross_product with this ""random"" vector to produce a perpendicular vector
  // Think about this: if the forward vector is (0,0,1), then the right vector ought to be (1,0,0). This can only be done if we choose as our arbitrary vector, the vector (0,1,0).
  vec3 right = cross(normalize(tmp), forward);// X axis
  vec3 up = cross(forward, right); // Y axis
  //convert to vec4
  vec4 forward_4 (forward[0],forward[1],forward[2],0);
  vec4 right_4(right[0],right[1],right[2],0);
  vec4 up_4(up[0],up[1],up[2],0);
  mat4 cam_coordinates(right_4,up_4,forward_4,vec4(from[0],from[1],from[2],1));
  return(cam_coordinates);
  ////niiice
}
vec3 DirectLight( const Intersection& i, Light dLight,const vector<Triangle>& triangles ){

  float r = glm::distance(dLight.position, i.position);
  float A= 4*M_PI*(r*r);
  vec3 B= dLight.color/A;
  vec4 normal = normalize(triangles[i.triangleIndex].normal);//normal pointing out of the surface/triangle
  vec4 direction = normalize(dLight.position - i.position );//direction from triangle to light source
  float dot_d_n = dot(direction, normal);
  vec3 D= B*max(dot_d_n,0.f);
  float n_mag=sqrt(pow(normal[0],2)+pow(normal[1],2)+pow(normal[2],2)+pow(normal[3],2));
  float d_mag=sqrt(pow(direction[0],2)+pow(direction[1],2)+pow(direction[2],2)+pow(direction[3],2));
  float angle_D_N =acos(dot_d_n/(d_mag*n_mag)) *180/ M_PI;

  //CHECK IF LIGHT IS THE CLOSEST INTERSECTION. IF NOT, DISPLAY SHADOW
  Intersection inter;
  if (closestIntersection(i.position,direction,triangles,inter,camera)){
    if((glm::distance(i.position,inter.position)<r &&glm::distance(i.position,inter.position)>0.001) || angle_D_N>90) return(vec3(0,0,0)); // used the second part as we might intersect with the same triangle(maybe), or a triangle very close. used 0.001 for now, need to figure out if there is any exact number
  }// NEED TO FIND THE CORRECT DISTANCE AND NOT 0.001

  return (D);
}



/* QUESTIONS: WHY DO WE NORMALIZE THE VECTORS?
              EXPLANATION OF LOOKAT AND HOW TO USE IT*/
