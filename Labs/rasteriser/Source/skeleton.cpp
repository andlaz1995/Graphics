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
#define PI 3.14159265

/* GLOBAL VARIABLE*/
std::vector<Triangle> triangles;
vec4 camera_pos(0, 0, -3.001,1);
mat4 R;
float f = SCREEN_HEIGHT;
float yaw = (0.0f * 3.1415926 / 180); // Yaw angle controlling camera rotation around y-axis
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH];
static vec3 anti_aliasing[SCREEN_WIDTH / 2][SCREEN_HEIGHT / 2];
static vec3 original_img[SCREEN_WIDTH][SCREEN_HEIGHT];
vec3 lightPos(0,-0.5,-3.5);

vec3 lightPower = 14.f*vec3( 1, 1, 1 );
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );
Triangle *currentTriangle;

struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 illumination;
    vec4 position;
    vec2 texturePosition;
    //int triangle_index;

};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS                                                                   */

bool Update();
void Draw(screen* screen);
void VertexShader( const vec4& v, Pixel& p, int vertexnum, vec4 currentNormal);
void Interpolate( Pixel a, Pixel b, vector<Pixel>& result );
void DrawLineSDL( screen* screen, Pixel a, Pixel b, vec3 color );
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen );
void DrawPolygon( const vector<vec4>& vertices, vec3 color, screen* screen, vec4 currentNormal);
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels );
void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, vec3 color, screen* screen, vec4 currentNormal);
void PixelShader( const Pixel& p, vec3 color, screen* screen, vec4 currentNormal );

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
    
    //memset(screen->buffer, 0, screen->height*screen->width*sizeof(uint32_t));
   memset(depthBuffer, 0, sizeof(depthBuffer));
   
   //initial depth buffer
   
    for( int y=0; y<SCREEN_HEIGHT; ++y )
    for( int x=0; x<SCREEN_WIDTH; ++x )
    depthBuffer[y][x] = 0;
    
    LoadTestModel(triangles);
    R = mat4(cos(yaw), 0, sin(yaw), 0, 0, 1, 0, 0, -sin(yaw), 0, cos(yaw),0 ,0 ,0, 0, -1);
    for( uint32_t i=0; i<triangles.size(); ++i )
    {
        vector<vec4> vertices(3);
        vertices[0]= triangles[i].v0;
        vertices[1]= triangles[i].v1;
        vertices[2]= triangles[i].v2;
        vec3 color  = triangles[i].color;
        
        vec4 currentNormal = triangles[i].normal;
        
        currentTriangle = &triangles[i];
        

        
          //DrawPolygonEdges (vertices, screen);
          DrawPolygon(vertices, color, screen, currentNormal);

        //for(int v=0; v<3; ++v)
        //{
            /*Pixel projPos;
           
            VertexShader( vertices[v], projPos );
            vec3 color(1,1,1);*/
            //PutPixelSDL( screen, projPos.x, projPos.y, color );
            //DrawPolygonEdges( vertices [v] );
           // projPos1 = projPos;
       //}
    }
}
void VertexShader(const vec4& v, Pixel& projPos, int vertexnum, vec4 currentNormal) //this function is to covert 3D to 2D vertices
{
    vec4 P_; //P_ (Origin of the camera)
    P_ = (v - camera_pos) * R;

    float x_ = (P_[0]) / (P_[2]) * f + SCREEN_WIDTH / 2;
    float y_ = (P_[1]) / (P_[2]) * f + SCREEN_HEIGHT / 2;
    float z_ = 1 / (P_[2]); //inverse function for depth
    
    //projPos.x= x_;
    //projPos.y= y_;

    projPos.x = x_;
    projPos.y = y_;
    projPos.zinv = z_;
    projPos.position.x=P_.x;
    projPos.position.y=P_.y;
    projPos.position.z=P_.z;
    projPos.position.w=P_.w;
     //projPos.position = P_;
    //p.triangle_index = triangle_index;
  
           /*vec4 light;
  
    
           vec3 D;
           float k;
           
           light = vec4(lightPos.x - camera_pos.x, lightPos.y - camera_pos.y, lightPos.z - camera_pos.z, 1);
           light = R * light; //R is rotation
           
           vec4 r = vec4(light.x - projPos.position.x, light.y - projPos.position.y, light.z - projPos.position.z, 1);
           float r_2 = r.x * r.x + r.y * r.y + r.z * r.z;
           
           k = currentNormal.x * r.x + currentNormal.y * r.y + currentNormal.z * r.z;
           k = (k > 0)?k : 0;
           k = k / (4 * PI * r_2);
           
           D = k * lightPower;
           projPos.illumination = (D + indirectLightPowerPerArea);*/
           
        

   
}

void Interpolate( Pixel a, Pixel b, vector<Pixel>& result )
{
int N = result.size();
//vec2 step = vec2(b-a) / float(max(N-1,1));
//vec2 current( a );

    float step_x = float(b.x - a.x) / float(max(N - 1, 1));
    float step_y = float(b.y - a.y) / float(max(N - 1, 1));
    float step_zinv = float(b.zinv - a.zinv) / float(max(N - 1, 1));
        
    for (int i = 0; i < N; i++)
    {
    
        result[i].x = a.x + step_x * i;
        result[i].y = a.y + step_y * i;
        result[i].zinv = a.zinv + step_zinv * i;

        result[i].position[2] = 1 / result[i].zinv + camera_pos[2];
        result[i].position[1] = (result[i].y - SCREEN_HEIGHT / 2) / f * (result[i].position[2] - camera_pos[2]) + camera_pos[1];
        result[i].position[0] = (result[i].x - SCREEN_WIDTH / 2) / f * (result[i].position[2] - camera_pos[2]) + camera_pos[0];
        result[i].illumination = a.illumination;
                   
    }
    
/*for( int i=0; i<N; ++i )
{
result[i] = round(current);
current += step;
} */

}

void DrawLineSDL( screen* screen, Pixel a, Pixel b, vec3 color )
{
    int delta_x = glm::abs(a.x - b.x);
    int delta_y = glm::abs(a.y - b.y);
    //Pixel delta = glm::abs( a - b );
    int pixels = glm::max( delta_x, delta_y ) + 1;
    
    
    vector<Pixel> line(pixels);
    Interpolate( a, b, line );
    
    for (int i=0; i<pixels; i++)
    {
    PutPixelSDL(screen, line[i].x, line[i].y, color);
    }

    
    /*vec3 light_area;
    Intersection inter, shadow_inter;
    // line = vector<Pixel> (pixels);*/


   /*for (int i = 0; i < pixels; i++)
    {
        PutPixelSDL(screen, line[i].x, line[i].y, color);
        if (line[i].x >= 0 && line[i].x < SCREEN_WIDTH  && line[i].y >= 0 && line[i].y < SCREEN_HEIGHT)
        {
            if (line[i].zinv > depthBuffer[line[i].x][line[i].y])
            {
                depthBuffer[line[i].x][line[i].y] = line[i].zinv;
                }
        }
    }
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
vector<Pixel> projectedVertices( V );
for( int i=0; i<V; ++i )
{
//VertexShader( vertices[i], projectedVertices[i], i ); //to be uncommented
}
// Loop over all vertices and draw the edge from it to the next vertex:
for( int i=0; i<V; ++i )
{
int j = (i+1)%V; // The next vertex
vec3 color( 1, 1, 1 );

DrawLineSDL (screen, projectedVertices[i], projectedVertices[j], color );
}
}

void DrawPolygon( const vector<vec4>& vertices, vec3 color, screen* screen, vec4 currentNormal )
   {
       int V = vertices.size();
       vector<Pixel> vertexPixels( V );
       for( int i=0; i<V; ++i )
       {
           VertexShader( vertices[i], vertexPixels[i], i, currentNormal );
           
          
        }
       vector<Pixel> leftPixels;
       vector<Pixel> rightPixels;
        vec3 acu_, acu;

       ComputePolygonRows( vertexPixels, leftPixels, rightPixels );
       DrawPolygonRows( leftPixels, rightPixels, color, screen, currentNormal);
       
    /*int height = 0;
    for (int i = 0; i < SCREEN_HEIGHT; i = i + 2)
    {
        int width = 0;
        for (int j = 0; j < SCREEN_WIDTH; j = j + 2)
        {
            acu_ = vec3(0.0f, 0.0f, 0.0f);

            if (depthBuffer[j][i] > 0 && (1.0f / depthBuffer[j][i]) > 2.5)
            {
                for (int k = 0; k < 4; k++)
                {
                    int depth = (((1.0f / depthBuffer[j][i]) - 2.5) * 4.0f);
                    depth = depth - depth / 2 + depth / 2 * (rand() / float(RAND_MAX));
                    if (j + depth < SCREEN_WIDTH && j - depth >= 0 && i + depth < SCREEN_HEIGHT && i - depth >= 0) {

                        acu = (original_img[j][i] +
                               original_img[j + depth][i] +
                               original_img[j][i + depth] +
                               original_img[j + depth][i + depth] +
                               original_img[j - depth][i] +
                               original_img[j][i - depth] +
                               original_img[j - depth][i - depth] +
                               original_img[j + depth][i - depth] +
                               original_img[j - depth][i + depth]) / vec3(9.0f, 9.0f, 9.0f);
                    }
                    acu_ = acu_ + acu;
                }

                anti_aliasing[width][height] = acu_ / 4.0f;
                PutPixelSDL(screen, width, height, anti_aliasing[width][height]);
                anti_aliasing[width][height] = vec3(0, 0, 0);
                width++;
            }
            else
            {
                anti_aliasing[width][height] = (original_img[j][i] + original_img[j + 1][i] + original_img[j][i + 1] + original_img[j + 1][i + 1]) / vec3(4.0f, 4.0f, 4.0f);
                PutPixelSDL(screen, width, height, anti_aliasing[width][height]);
                anti_aliasing[width][height] = vec3(0, 0, 0);
                width++;
            }
        }
        height++;
    } */
}

void ComputePolygonRows(
       const vector<Pixel>& vertexPixels,
       vector<Pixel>& leftPixels,
       vector<Pixel>& rightPixels ) // Trying to compute the start and end postion of array 
{

   // 1. Find max and min y-value of the polygon
   //    and compute the number of rows it occupies.
    int min = numeric_limits<int>::max();
    int max = 0;
    vector<Pixel> edge;

    for (size_t i = 0; i < vertexPixels.size(); i++)
    {
        if (min > vertexPixels[i].y)
            min = vertexPixels[i].y;
        if (max < vertexPixels[i].y)
            max = vertexPixels[i].y;
    }
    int ROWS= max - min + 1; //To calculate the rows (40 − 10 + 1 = 31 rows) 40=max value and 10=min value in notes
   
   // 2. Resize leftPixels and rightPixels so that they have an element for each row.
   // for computing left and right vertices
    
    leftPixels = vector<Pixel>(ROWS);
    rightPixels = vector<Pixel>(ROWS);
    
  // 3. Initialize the x-coordinates in leftPixels to some really large value and the x-coordinates in rightPixels to some really small value.
    
    for (int i = 0; i < ROWS; i++)
    {
        leftPixels[i].x = +numeric_limits<int>::max();
        rightPixels[i].x = -numeric_limits<int>::max();
        leftPixels[i].y = min + i;
        rightPixels[i].y = min + i;
    }
    
   // 4. Loop through all edges of the polygon and use
   //    linear interpolation to find the x-coordinate for
   //    each row it occupies. Update the corresponding
// values in rightPixels and leftPixels. 
    for (size_t i = 0; i < vertexPixels.size(); i++)
    {
        int j = (i + 1) % vertexPixels.size(); // The next vertex
        // DrawLineSDL(screen, vertexPixels[i], vertexPixels[j], color, edge);
        int delta_x = glm::abs(vertexPixels[i].x - vertexPixels[j].x);
        int delta_y = glm::abs(vertexPixels[i].y - vertexPixels[j].y);
        int pixels = glm::max(delta_x, delta_y) + 1;
        edge = vector<Pixel> (pixels);
        
        
        Interpolate(vertexPixels[i], vertexPixels[j], edge);
        
        for (int a = 0; a < ROWS; a++)
        {
            for (size_t b = 0; b < edge.size(); b++)
            {
                if (edge[b].y == min + a)
                {
                    if (edge[b].x < leftPixels[a].x)
                    {
                        leftPixels[a].x = edge[b].x;
                        leftPixels[a].zinv = edge[b].zinv;
                        leftPixels[a].position= edge[b].position;
                        leftPixels[a].texturePosition=edge[b].texturePosition;
                        leftPixels[a].illumination=edge[b].illumination;
                        
                        
                    }
                    if (edge[b].x > rightPixels[a].x)
                    {
                        rightPixels[a].x = edge[b].x;
                        rightPixels[a].zinv = edge[b].zinv;
                        rightPixels[a].position= edge[b].position;
                        rightPixels[a].illumination=edge[b].illumination;
                    }
                }
            }
        }
    }
}

void DrawPolygonRows(vector<Pixel>& leftPixels, vector<Pixel>& rightPixels, vec3 color, screen* screen,vec4 currentNormal)
{

int size = leftPixels.size();
  vector<Pixel> line;
 
  for(int i = 0; i < size; i++){
    ivec2 delta = glm::abs(ivec2(leftPixels[i].x - rightPixels[i].x, leftPixels[i].y - rightPixels[i].y));
    
    int pixels = glm::max(delta.x, delta.y) + 1;

    line.clear();
    line.resize(pixels);
        

    Interpolate(leftPixels[i], rightPixels[i], line);

    for(int j = 0; j < pixels; j++){
      PixelShader(line[j], color, screen, currentNormal);
    }
  }
  
  
    /*for (size_t i = 0; i < leftPixels.size(); i++) {

        DrawLineSDL(screen, leftPixels[i], rightPixels[i], color);
        PixelShader(leftPixels[i], color, screen );
    }*/
}

void PixelShader( const Pixel& p, vec3 color, screen* screen, vec4 currentNormal)
{
  
  int x = p.x;
  int y = p.y;
  vec4 light;
  vec3 illumination;
  
       if( p.zinv > depthBuffer[y][x] )
       {
  
    
          vec3 D;
           float k;
           
           light = vec4(lightPos.x - camera_pos.x, lightPos.y - camera_pos.y, lightPos.z - camera_pos.z, 1);
           light = R * light; //R is rotation
           
           vec4 r = vec4(light.x - p.position.x, light.y - p.position.y, light.z - p.position.z, 1);
           float r_2 = r.x * r.x + r.y * r.y + r.z * r.z;
           
           k = currentNormal.x * r.x + currentNormal.y * r.y + currentNormal.z * r.z;
           k = (k > 0)?k : 0;
           k = k / (4 * PI * r_2);
           
           D = k * lightPower;
           illumination = (D + indirectLightPowerPerArea);
          
           
           //p.illumination = illumination;

           depthBuffer[y][x] = p.zinv;
           color.x =  illumination.x * color.x;
           color.y =  illumination.y * color.y;
           color.z  = illumination.z * color.z;
           
           PutPixelSDL( screen, x, y, color);
           
       }
}


/*Place updates of parameters here*/
bool Update()
{
  static int t = SDL_GetTicks();
  /* Compute frame time */
  int t2 = SDL_GetTicks();
  //float dt = float(t2-t);
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
