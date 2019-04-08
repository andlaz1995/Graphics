#include <iostream>
#include <glm/glm.hpp>
#include <SDL.h>
#include "SDLauxiliary.h"
#include "TestModelH.h"
#include <stdint.h>

using namespace std;
using glm::vec3;
using glm::mat3;
using glm::vec4;
using glm::mat4;
using glm::ivec2;
using glm::vec2;

SDL_Event event;

#define SCREEN_WIDTH 320
#define SCREEN_HEIGHT 256
#define FULLSCREEN_MODE false

/* GLOBAL VARIABLE*/
std::vector<Triangle> triangles;
vec4 camera_pos(0, 0, -3.001,1);
mat4 R;
float f = SCREEN_HEIGHT;
float yaw = (0.0f * 3.1415926 / 180); // Yaw angle controlling camera rotation around y-axis

//SDL_Surface* screen;

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, ivec2& p );
void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result );
void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color );
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen );


int main( int argc, char* argv[] )
{
  
  screen* screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );

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
/*  memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
  
  vec3 colour(1.0,0.0,0.0);
  for(int i=0; i<1000; i++)
    {
      uint32_t x = rand() % screen->width;
      uint32_t y = rand() % screen->height;
      PutPixelSDL(screen, x, y, colour);
    }
    */
    
    memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
    LoadTestModel(triangles);
    R = mat4(cos(yaw), 0, sin(yaw), 0, 0, 1, 0, 0, -sin(yaw), 0, cos(yaw),0 ,0 ,0, 0, -1);
    for( uint32_t i=0; i<triangles.size(); ++i )
    {
        vector<vec4> vertices(3);
        vertices[0] = triangles[i].v0;
        vertices[1] = triangles[i].v1;
        vertices[2] = triangles[i].v2;
        
         DrawPolygonEdges (vertices, screen);

        for(int v=0; v<3; ++v)
        {
            /*ivec2 projPos;
           
            VertexShader( vertices[v], projPos );
            vec3 color(1,1,1);*/
            //PutPixelSDL( screen, projPos.x, projPos.y, color );
            //DrawPolygonEdges( vertices [v] );
           // projPos1 = projPos;
        }
    }
}
void VertexShader(const vec4& v, ivec2& projPos)
{
    vec4 P_; //P_ (Origin of the camera)
    P_ = (v - camera_pos) * R;

    float x_ = (P_[0]) / (P_[2]) * f + SCREEN_WIDTH / 2;
    float y_ = (P_[1]) / (P_[2]) * f + SCREEN_HEIGHT / 2;
    projPos.x= x_;
    projPos.y= y_;

   /* p.x = x_;
    p.y = y_;
    p.zinv = z_;
    p.pos3d = v;
    p.triangle_index = triangle_index;*/
}

void Interpolate( ivec2 a, ivec2 b, vector<ivec2>& result )
{
int N = result.size();
vec2 step = vec2(b-a) / float(max(N-1,1));
vec2 current( a );
for( int i=0; i<N; ++i )
{
result[i] = round(current);
current += step;
}
}

void DrawLineSDL( screen* screen, ivec2 a, ivec2 b, vec3 color )
{

    ivec2 delta = glm::abs( a - b );
    int pixels = glm::max( delta.x, delta.y ) + 1;
    
    
    vector<ivec2> line( pixels );
    Interpolate( a, b, line );
    
    for (int i=0; i<pixels; i++)
    {
    PutPixelSDL(screen, line[i].x, line[i].y, color);
    }

    
    /*vec3 light_area;
    Intersection inter, shadow_inter;
    // line = vector<Pixel> (pixels);*/


   /* for (int i = 0; i < pixels; i++)
    {
        // PutPixelSDL(surface, line[i].x, line[i].y, color);
        if (line[i].x >= 0 && line[i].x < SCREEN_WIDTH  && line[i].y >= 0 && line[i].y < SCREEN_HEIGHT)
        {
            if (line[i].zinv > depthBuffer[line[i].x][line[i].y])
            {
                depthBuffer[line[i].x][line[i].y] = line[i].zinv;
                vec3 dis = light_pos - line[i].pos3d;
                float r = glm::length(dis);

                float result = dis[0] * triangles[line[i].triangle_index].normal[0] + dis[1] * triangles[line[i].triangle_index].normal[1] + dis[2] * triangles[line[i].triangle_index].normal[2];
                float camera_pos1 = 4.0 * 3.1415926 * r * r;

                if (result > 0.0)
                    light_area = result / camera_pos1 * lightPower;
                else
                    light_area = vec3(0.0, 0.0, 0.0);

                if (light_pos[1] < -0.20f)
                {
                    if (line[i].pos3d[1] >= -0.20f)
                    {
                        if (closest_intersection(line[i].pos3d, dis, triangles, inter))
                        {
                            vec3 dis_ = inter.position - line[i].pos3d;
                            if (r > glm::length(dis_) && result > 0.0 && line[i].triangle_index != inter.triangle_index)
                                light_area = vec3(0.0, 0.0, 0.0);
                        }

                        light_area = 0.5f * (indirectLightPowerPerArea + light_area);
                        original_img[line[i].x][line[i].y] = color * light_area;
                    }
                    else
                    {

                        light_area = 0.5f * (indirectLightPowerPerArea + light_area);
                        original_img[line[i].x][line[i].y] = color * light_area;
                    }
                }
                else
                {

                    if (closest_intersection(line[i].pos3d, dis, triangles, inter))
                    {
                        vec3 dis_ = inter.position - line[i].pos3d;
                        if (r > glm::length(dis_) && result > 0.0 && line[i].triangle_index != inter.triangle_index)
                            light_area = vec3(0.0, 0.0, 0.0);
                    }

                    light_area = 0.5f * (indirectLightPowerPerArea + light_area);
                    original_img[line[i].x][line[i].y] = color * light_area;
                }
            }
        }
    }*/
}

void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen )

{
int V = vertices.size();
// Transform each vertex from 3D world position to 2D image position:
vector<ivec2> projectedVertices( V );
for( int i=0; i<V; ++i )
{
VertexShader( vertices[i], projectedVertices[i] );
}
// Loop over all vertices and draw the edge from it to the next vertex:
for( int i=0; i<V; ++i )
{
int j = (i+1)%V; // The next vertex
vec3 color( 1, 1, 1 );

DrawLineSDL (screen, projectedVertices[i], projectedVertices[j], color );
}
}


/*Place updates of parameters here*/
bool Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  float dt = float(t2-t);
  t = t2;

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
		/* Move camera forward */
		break;
	      case SDLK_DOWN:
		/* Move camera backwards */
		break;
	      case SDLK_LEFT:
		/* Move camera left */
		break;
	      case SDLK_RIGHT:
		/* Move camera right */
		break;
	      case SDLK_ESCAPE:
		/* Move camera quit */
		return false;
	      }
	  }  
    }
  return true;
}
