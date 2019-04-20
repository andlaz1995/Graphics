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
float depthBuffer[SCREEN_HEIGHT][SCREEN_WIDTH]; //to provide defect to objects
static vec3 anti_aliasing[SCREEN_WIDTH / 2][SCREEN_HEIGHT / 2];
static vec3 original_img[SCREEN_WIDTH][SCREEN_HEIGHT];
vec3 lightPos(0,-0.5,-0.7); // determines the position of light 
vec3 lightPower = 14.f*vec3( 1, 1, 1 );// power of light
vec3 indirectLightPowerPerArea = 0.5f*vec3( 1, 1, 1 );// used to calculate illumination
Triangle *currentTriangle;
SDL_Surface *textureSurface;

//struct pixel - stores the vetrex, pixel, illumination of each vertex of every triangle
struct Pixel
{
    int x;
    int y;
    float zinv;
    vec3 illumination;
    vec4 position;
    vec2 texturePosition;
    //int triangle_index

};

/* ----------------------------------------------------------------------------*/
/* FUNCTIONS DECLRATION                                                             */

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
  //initialise screen with width and height
  screen* screen = InitializeSDL( SCREEN_WIDTH, SCREEN_HEIGHT, FULLSCREEN_MODE );
  
  const char file[] = "./space.bmp";
  textureSurface = SDL_LoadBMP(file);
  if(textureSurface == NULL){
    cout << "Fail to load texture" << endl;
  } else {
    cout << "image loaded" << endl;
  }
 
//update() and draw() is called until user close the screen
  while ( Update())
    {
      Draw(screen); 
      SDL_Renderframe(screen);
    }

  SDL_SaveImage( screen, "screenshot.bmp" );// output image saved with this name

  KillSDL(screen);
  return 0;
}

/*Place your drawing here*/
void Draw(screen* screen)
{
	//to spot dot on the screen
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
    //clear the depth buffer
   	memset(depthBuffer, 0, sizeof(depthBuffer));
   
   //initialize depth buffer   
    for( int y=0; y<SCREEN_HEIGHT; ++y )
    for( int x=0; x<SCREEN_WIDTH; ++x )
    depthBuffer[y][x] = 0;
    
    //Load the triangles - defined in TestModelH.h
    LoadTestModel(triangles);
    //initialise Rotation angle 
    R = mat4(cos(yaw), 0, sin(yaw), 0, 0, 1, 0, 0, -sin(yaw), 0, cos(yaw),0 ,0 ,0, 0, -1);
    
    //traverse every traingle and store the vertex, color, normal informations to Draw 
    for( uint32_t i=0; i<triangles.size(); ++i )
    {
        vector<vec4> vertices(3);//vertor to store 3 vertexes of triangle
        vertices[0]= triangles[i].v0;
        vertices[1]= triangles[i].v1;
        vertices[2]= triangles[i].v2;
        vec3 color  = triangles[i].color;
        vec4 currentNormal = triangles[i].normal;
        currentTriangle = &triangles[i];
        
        //DrawPolygonEdges (vertices, screen); // this function is called when raw polygon is draw with triangle information
        DrawPolygon(vertices, color, screen, currentNormal);// draw polygon with color and shapes

    }
}

/***********       Function to convert 3D vertices to 2D vertices     *******************/
void VertexShader(const vec4& v, Pixel& projPos, int vertexnum, vec4 currentNormal) 
{
    vec4 P_; //P_ (Origin of the camera)
    P_ = (v - camera_pos) * R; //(Origin of the camera = postion of the vertices - position of camera)* rotation matrix (all are 4*4 matrices)

    float x_ = (P_[0]) / (P_[2]) * f + SCREEN_WIDTH / 2; //(X/Z * focallength +(W/2))
    float y_ = (P_[1]) / (P_[2]) * f + SCREEN_HEIGHT / 2; //(Y/Z * focallength +(H/2))
    float z_ = 1 / (P_[2]); //inverse function for depth

//For texture mapping

 // vec4 e1(R[0][0], R[0][1], R[0][2], 0);
  //vec4 e2(R[1][0], R[1][1], R[1][2], 0);
  //vec4 e3(R[2][0], R[2][1], R[2][2], 0);
  //vec4 e4(0, 0, 0, 1);
  //mat4 R_homo(e1, e2, e3, e4);
  //P_ = R_homo * P_;
   
 	//add the information to the pixel struct
    projPos.x = x_;
    projPos.y = y_;
    projPos.zinv = z_;
    projPos.position.x=P_.x;
    projPos.position.y=P_.y;
    projPos.position.z=P_.z;
    projPos.position.w=P_.w;
    //projPos.position = P_;
    //p.triangle_index = triangle_index;
    
    //texture
  if(vertexnum == 0){
    projPos.texturePosition = currentTriangle->v0_tp;
    //std::cout << "apa3" << std::endl;
  } else if(vertexnum == 1){
    projPos.texturePosition = currentTriangle->v1_tp;
    //std::cout << "apa4" << std::endl;
  } else if(vertexnum == 2){
    projPos.texturePosition = currentTriangle->v2_tp;
    //std::cout << "apa5" << std::endl;
  }
 
    /* Per vertex illumination
    vec4 light;
	vec3 D;
	float k;
       
    light = vec4(lightPos.x - camera_pos.x, lightPos.y - camera_pos.y, lightPos.z - camera_pos.z, 1);
    light = R * light; //R is rotation
    vec4 r = vec4(light.x - projPos.position.x, light.y - projPos.position.y, light.z - projPos.position.z, 1);
    float r_2 = r.x * r.x + r.y * r.y + r.z * r.z;
	//formula to calculate illumination
    k = currentNormal.x * r.x + currentNormal.y * r.y + currentNormal.z * r.z;
	k = (k > 0)?k : 0;
    k = k / (4 * PI * r_2);
           
    D = k * lightPower;
	projPos.illumination = (D + indirectLightPowerPerArea);*/
}


void Interpolate(Pixel a, Pixel b, vector<Pixel> &result){
  int N = result.size();
  vec2 step = vec2(b.x - a.x, b.y - a.y) / float(max(N-1, 1));
  float depth = (b.zinv - a.zinv) / float(max(N-1, 1));

  vec4 pos = vec4(b.position * b.zinv - a.position * a.zinv) / float(max(N - 1, 1));
  vec2 texture_step = vec2(b.texturePosition - a.texturePosition) / float(max(N - 1, 1));
  Pixel current = a;

  vec4 q = current.position * current.zinv;

  for(int i = 0; i < N; i++){
    result[i] = current;
    current.x += step.x;
    current.y += step.y;
    current.zinv += depth;
    q += pos;
    current.position = q / current.zinv;
    current.texturePosition += texture_step;
    current.illumination = a.illumination;
    result[i] = current;
    // printf("result.zinv %f\n", result[i].zinv);
  }
}

/**********       Function to find interpolation to draw line   ************/
/*void Interpolate( Pixel a, Pixel b, vector<Pixel>& result ) //interpolate vector is stored in vector<Pixel>& result
{
	int N = result.size();//before finding interpolate - identify no of pixels between the vector

    float step_x = float(b.x - a.x) / float(max(N - 1, 1));
    float step_y = float(b.y - a.y) / float(max(N - 1, 1));
    float step_zinv = float(b.zinv - a.zinv) / float(max(N - 1, 1));
    
    vec4 pos = vec4(b.position * b.zinv - a.position * a.zinv) / float(max(N - 1, 1));
    vec2 texture_step = vec2(b.texturePosition - a.texturePosition) / float(max(N - 1, 1));
    Pixel current = a;

    vec4 q = current.position * current.zinv;
   
    //stores pixel information to draw line    
   
    for (int i = 0; i < N; i++)
    {
    
        result[i].x += step_x;
        result[i].y += step_y;
        result[i].zinv += step_zinv;
        
        q+= pos;
        
      //  result[i].position = q / current.zinv;
        result[i].texturePosition += texture_step;
        
          

        result[i].position[2] = 1 / result[i].zinv + camera_pos[2];
        result[i].position[1] = (result[i].y - SCREEN_HEIGHT / 2) / f * (result[i].position[2] - camera_pos[2]) + camera_pos[1];
        result[i].position[0] = (result[i].x - SCREEN_WIDTH / 2) / f * (result[i].position[2] - camera_pos[2]) + camera_pos[0];
        result[i].illumination = a.illumination;
                   
    }
}*/


/********* Function to daw line with the information from interpolate function         **************/
void DrawLineSDL( screen* screen, Pixel a, Pixel b, vec3 color )
{
    int delta_x = glm::abs(a.x - b.x);//used when pixel information is used to draw polygon
    int delta_y = glm::abs(a.y - b.y);
    
    //Pixel delta = glm::abs( a - b );//used when vertex information is used draw polygon
    int pixels = glm::max( delta_x, delta_y ) + 1; //before finding interplote - identify no of pixels between the vector
     
    vector<Pixel> line(pixels);
    Interpolate( a, b, line ); // find 'line' no of points between 'a' and 'b'
    
    for (int i=0; i<pixels; i++)
    {
   // PutPixelSDL(screen, line[i].x, line[i].y, color); // draw line with interpolate points
    }

}

/***********  Function to draw just the edges of polygon                  ************/
void DrawPolygonEdges( const vector<vec4>& vertices, screen* screen )
{
	int V = vertices.size();
	// Transform each vertex from 3D world position to 2D image position:
	vector<Pixel> projectedVertices( V );
	//for( int i=0; i<V; ++i )
	//{
	//	VertexShader( vertices[i], projectedVertices[i], i, currentNormal ); //to be uncommented
	//}
	// Loop over all vertices and draw the edge from it to the next vertex:
	for( int i=0; i<V; ++i )
	{
		int j = (i+1)%V; // The next vertex
		vec3 color( 1, 1, 1 );
		DrawLineSDL (screen, projectedVertices[i], projectedVertices[j], color );//draw line
	}
}


/****************vFunction to draw polygon with color    ****************/
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

    ComputePolygonRows( vertexPixels, leftPixels, rightPixels ); // compute rows(y-axis) and  left pixels & right pixels(y-axis)
    DrawPolygonRows( leftPixels, rightPixels, color, screen, currentNormal);    // draw polygon with computed rows in above function
   
}

/*******************       to compute no of rows           *****************/
void ComputePolygonRows(const vector<Pixel>& vertexPixels, vector<Pixel>& leftPixels, vector<Pixel>& rightPixels ) // To compute the start and end postion of array 
{
/* // 1. Find max and min y-value of the polygon and compute the number of rows it occupies.
    int V= vertexPixels.size();
    int min = numeric_limits<int>::max();
    int max = 0;
    vector<Pixel> edge;
    for (size_t i = 0; i < vertexPixels.size(); i++)
    {
        if (min > vertexPixels[i].y)
            min = vertexPixels[i].y; //10
        if (max < vertexPixels[i].y)
            max = vertexPixels[i].y;//40
    }
    int ROWS= max - min + 1; //To calculate the rows (40 − 10 + 1 = 31 rows) 40=max value and 10=min value in notes
   
// 2. Resize leftPixels and rightPixels so that they have an element for each row for computing left and right vertices
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
    
// 4. Loop through all edges of the polygon and use linear interpolation to find the x-coordinate for each row it occupies. Update the corresponding
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
                       // leftPixels[a].texturePosition = a.texturePosition + (b.texturePosition - a.texturePosition);
                        leftPixels[a].illumination=edge[b].illumination;
                        
                        
                    }
                    if (edge[b].x > rightPixels[a].x)
                    {
                        rightPixels[a].x = edge[b].x;
                        rightPixels[a].zinv = edge[b].zinv;
                        rightPixels[a].position= edge[b].position;
                       // rightPixels[a].texturePosition = a.texturePosition + (b.texturePosition - a.texturePosition);
                        rightPixels[a].illumination=edge[b].illumination;
                    }
                }
            }
        }
    }*/
    
    
   int V = vertexPixels.size();

   //find the max and min y
   int max = std::max(vertexPixels[0].y, vertexPixels[1].y);
   max = std::max(max, vertexPixels[2].y);
   int min = std::min(vertexPixels[0].y, vertexPixels[1].y);
   min = std::min(min, vertexPixels[2].y);
   int rows = max - min + 1;

   //reserve size
   leftPixels.resize(rows);
   rightPixels.resize(rows);

   //initialize x in left to largest, in right to smallest
   for(int i = 0; i < rows; i++){
    leftPixels[i].x = numeric_limits<int>::max();
    rightPixels[i].x = -numeric_limits<int>::max();

  }
    for(int i = 0; i < V; i++){
    Pixel a = vertexPixels[i];
    Pixel b = vertexPixels[(i + 1) % V];
    float step = 1.f / (glm::length(vec2(a.x - b.x, a.y - b.y)) + 1);

    for(float t = 0; t < 1; t += step){
      Pixel pixel;

      pixel.x = static_cast<int>(a.x + t * (b.x - a.x));
      pixel.y = static_cast<int>(a.y + t * (b.y - a.y));
      pixel.zinv = a.zinv + t * (b.zinv - a.zinv);
      pixel.position = ((a.position * a.zinv) + t * (b.position * b.zinv - a.position * a.zinv)) / pixel.zinv;
      pixel.texturePosition = a.texturePosition + t * (b.texturePosition - a.texturePosition);
      

      if(pixel.x < 0){
        pixel.x = 0;
      } else if (pixel.x > SCREEN_WIDTH) {
        pixel.x = SCREEN_WIDTH;
      }

      if(pixel.y < 0){
        pixel.y = 0;
      } else if (pixel.y > SCREEN_HEIGHT) {
        pixel.y = SCREEN_HEIGHT;
      }

      //update left and right
      int y = pixel.y - min;
      if(pixel.x <= leftPixels[y].x){
        leftPixels[y] = pixel;
      }
      if(pixel.x >= rightPixels[y].x){
        rightPixels[y] = pixel;
      }
    }
  }
}


// draw the polygon with the above computed methods
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
      		PixelShader(line[j], color, screen, currentNormal); // pixel shader is introduced to give depth effect and per pixel illumination
      		
    }
  }
  
}


// to give depth effect and per pixel illumination
void PixelShader( const Pixel& p, vec3 color, screen* screen, vec4 currentNormal)
{
  
  int x = p.x;
  int y = p.y;
  vec4 light;
  vec3 illumination;
  vec3 texturecolor;
  
       if( p.zinv > depthBuffer[y][x] )
       {
       		//per pixel illumination
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
          
           depthBuffer[y][x] = p.zinv;//depth
           
           //calculate color with illumination
           color.x =  illumination.x * color.x;
           color.y =  illumination.y * color.y;
           color.z  = illumination.z * color.z;
           
          texturecolor = color;
          if(currentTriangle->texture){
          texturecolor = GetPixelSDL(textureSurface, (int)(p.texturePosition.x), (int)(p.texturePosition.y)); 
         // std::cout << p.texturePosition.x << std::endl; 
          //std::cout << p.texturePosition.y << std::endl;
          
    }
           
           PutPixelSDL( screen, x, y, illumination * texturecolor);
           
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
