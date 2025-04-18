/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "showCube.h"
#include <cmath>

int pointMap(int side, int i, int j)
{
  int r;

  switch (side)
  {
  case 1: //[i][j][0] bottom face
    r = 64 * i + 8 * j;
    break;
  case 6: //[i][j][7] top face
    r = 64 * i + 8 * j + 7;
    break;
  case 2: //[i][0][j] front face
    r = 64 * i + j;
    break;
  case 5: //[i][7][j] back face
    r = 64 * i + 56 + j;
    break;
  case 3: //[0][i][j] left face
    r = 8 * i + j;
    break;
  case 4: //[7][i][j] right face
    r = 448 + 8 * i + j;
    break;
  }

  return r;
}

// display inclined plane
void findOrthoVectors(const struct point & normal, struct point & u, struct point & v)
{
  struct point nabs;
  nabs.x = fabs(normal.x);
  nabs.y = fabs(normal.y);
  nabs.z = fabs(normal.z);

  struct point temp = {0,0,0};
  if (nabs.x <= nabs.y && nabs.x <= nabs.z)
    temp.x = 1.0;
  else if (nabs.y <= nabs.x && nabs.y <= nabs.z)
    temp.y = 1.0;
  else
    temp.z = 1.0;

  // u = normal x temp
  CROSSPRODUCTp(normal, temp, u);
  double length = sqrt(u.x*u.x + u.y*u.y + u.z*u.z);
  if (length < 1e-9) {
    u.x=0;u.y=1;u.z=0;
  } else {
    u.x /= length;
    u.y /= length;
    u.z /= length;
  }

  // v = normal cross u
  struct point vtemp;
  CROSSPRODUCTp(normal, u, vtemp);
  // normalize
  length = sqrt(vtemp.x*vtemp.x + vtemp.y*vtemp.y + vtemp.z*vtemp.z);
  v.x = vtemp.x/length;
  v.y = vtemp.y/length;
  v.z = vtemp.z/length;
}

void showInclinedPlane(const struct world * jello)
{
  if (jello->incPlanePresent != 1) return;

  // plane: a*x + b*y + c*z + d=0
  double a = jello->a;
  double b = jello->b;
  double c = jello->c;
  double d = jello->d;

  double length = sqrt(a * a + b * b + c * c);
  struct point normal = { a/length, b/length, c/length };

  // find orthonormal vectors u,v
  struct point u,v_;
  findOrthoVectors(normal, u, v_);

  double halfSize = 3.0; 
  double denom = length*length;
  double dist = -d/denom;
  struct point base;
  base.x = dist * a;
  base.y = dist * b;
  base.z = dist * c;

  // compute corners
  struct point corner[4];
  double coords[4][2] = {
    {-halfSize, -halfSize},
    {+halfSize, -halfSize},
    {+halfSize, +halfSize},
    {-halfSize, +halfSize}
  };

  for(int idx=0; idx<4; idx++){
    double s = coords[idx][0];
    double t = coords[idx][1];
    corner[idx].x = base.x + s*u.x + t*v_.x;
    corner[idx].y = base.y + s*u.y + t*v_.y;
    corner[idx].z = base.z + s*u.z + t*v_.z;
  }

  // set color or material
  glDisable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(1.0, 0.8, 0.8, 0.3);
  

  glBegin(GL_QUADS);
    glNormal3f(normal.x, normal.y, normal.z);
    glVertex3f(corner[0].x, corner[0].y, corner[0].z);
    glVertex3f(corner[1].x, corner[1].y, corner[1].z);
    glVertex3f(corner[2].x, corner[2].y, corner[2].z);
    glVertex3f(corner[3].x, corner[3].y, corner[3].z);
  glEnd();

  glDisable(GL_BLEND);
  glEnable(GL_LIGHTING);
}

void showCoordinates() {
  glDisable(GL_LIGHTING);
  glLineWidth(2.0f);

  // x
  glColor3f(0.5,0.5,1.0);
  glBegin(GL_LINES);
  glVertex3f(-0.5, 0.0, 0.0);
  glVertex3f(0.5, 0.0, 0.0);
  glVertex3f(0.5, 0.0, 0.0);
  glVertex3f(0.4, 0.1, 0.0);
  glVertex3f(0.5, 0.0, 0.0);
  glVertex3f(0.4, -0.1, 0.0);
  glEnd();

  // y
  glColor3f(1.0,0.5,0.5);
  glBegin(GL_LINES);
  glVertex3f(0.0, -0.5, 0.0);
  glVertex3f(0.0, 0.5, 0.0);
  glVertex3f(0.0, 0.5, 0.0);
  glVertex3f(0.1, 0.4, 0.0);
  glVertex3f(0.0, 0.5, 0.0);
  glVertex3f(-0.1, 0.4, 0.0);
  glEnd();

  // z
  glColor3f(0.5,1.0,0.5);
  glBegin(GL_LINES);
  glVertex3f(0.0, 0.0, -0.5);
  glVertex3f(0.0, 0.0, 0.5);
  glVertex3f(0.0, 0.0, 0.5);
  glVertex3f(0.0, 0.1, 0.4);
  glVertex3f(0.0, 0.0, 0.5);
  glVertex3f(0.0, -0.1, 0.4);
  glEnd();
  
  glEnable(GL_LIGHTING);
}

void showCube(struct world * jello)
{
  int i,j,k,ip,jp,kp;
  point r1,r2,r3; // aux variables
  
  /* normals buffer and counter for Gourad shading*/
  struct point normal[8][8];
  int counter[8][8];

  int face;
  double faceFactor, length;

  if (fabs(jello->p[0][0][0].x) > 10)
  {
    printf ("Your cube somehow escaped way out of the box.\n");
    exit(0);
  }

  
  #define NODE(face,i,j) (*((struct point * )(jello->p) + pointMap((face),(i),(j))))

  
  #define PROCESS_NEIGHBOUR(di,dj,dk) \
    ip=i+(di);\
    jp=j+(dj);\
    kp=k+(dk);\
    if\
    (!( (ip>7) || (ip<0) ||\
      (jp>7) || (jp<0) ||\
    (kp>7) || (kp<0) ) && ((i==0) || (i==7) || (j==0) || (j==7) || (k==0) || (k==7))\
       && ((ip==0) || (ip==7) || (jp==0) || (jp==7) || (kp==0) || (kp==7))) \
    {\
      glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);\
      glVertex3f(jello->p[ip][jp][kp].x,jello->p[ip][jp][kp].y,jello->p[ip][jp][kp].z);\
    }\

 
  if (viewingMode==0) // render wireframe
  {
    glLineWidth(1);
    glPointSize(5);
    glDisable(GL_LIGHTING);
    for (i=0; i<=7; i++)
      for (j=0; j<=7; j++)
        for (k=0; k<=7; k++)
        {
          if (i*j*k*(7-i)*(7-j)*(7-k) != 0) // not surface point
            continue;

          glBegin(GL_POINTS); // draw point
            glColor4f(0,0,0,0);  
            glVertex3f(jello->p[i][j][k].x,jello->p[i][j][k].y,jello->p[i][j][k].z);        
          glEnd();

          //
          //if ((i!=7) || (j!=7) || (k!=7))
          //  continue;

          glBegin(GL_LINES);      
          // structural
          if (structural == 1)
          {
            glColor4f(0,0,1,1);
            PROCESS_NEIGHBOUR(1,0,0);
            PROCESS_NEIGHBOUR(0,1,0);
            PROCESS_NEIGHBOUR(0,0,1);
            PROCESS_NEIGHBOUR(-1,0,0);
            PROCESS_NEIGHBOUR(0,-1,0);
            PROCESS_NEIGHBOUR(0,0,-1);
          }
          
          // shear
          if (shear == 1)
          {
            glColor4f(0,1,0,1);
            PROCESS_NEIGHBOUR(1,1,0);
            PROCESS_NEIGHBOUR(-1,1,0);
            PROCESS_NEIGHBOUR(-1,-1,0);
            PROCESS_NEIGHBOUR(1,-1,0);
            PROCESS_NEIGHBOUR(0,1,1);
            PROCESS_NEIGHBOUR(0,-1,1);
            PROCESS_NEIGHBOUR(0,-1,-1);
            PROCESS_NEIGHBOUR(0,1,-1);
            PROCESS_NEIGHBOUR(1,0,1);
            PROCESS_NEIGHBOUR(-1,0,1);
            PROCESS_NEIGHBOUR(-1,0,-1);
            PROCESS_NEIGHBOUR(1,0,-1);

            PROCESS_NEIGHBOUR(1,1,1)
            PROCESS_NEIGHBOUR(-1,1,1)
            PROCESS_NEIGHBOUR(-1,-1,1)
            PROCESS_NEIGHBOUR(1,-1,1)
            PROCESS_NEIGHBOUR(1,1,-1)
            PROCESS_NEIGHBOUR(-1,1,-1)
            PROCESS_NEIGHBOUR(-1,-1,-1)
            PROCESS_NEIGHBOUR(1,-1,-1)
          }
          
          // bend
          if (bend == 1)
          {
            glColor4f(1,0,0,1);
            PROCESS_NEIGHBOUR(2,0,0);
            PROCESS_NEIGHBOUR(0,2,0);
            PROCESS_NEIGHBOUR(0,0,2);
            PROCESS_NEIGHBOUR(-2,0,0);
            PROCESS_NEIGHBOUR(0,-2,0);
            PROCESS_NEIGHBOUR(0,0,-2);
          }           
          glEnd();
        }
    glEnable(GL_LIGHTING);
  }
  
  else
  {
    glPolygonMode(GL_FRONT, GL_FILL); 
    
    for (face=1; face <= 6; face++) 
      // face == face of a cube
      // 1 = bottom, 2 = front, 3 = left, 4 = right, 5 = far, 6 = top
    {
      
      if ((face==1) || (face==3) || (face==5))
        faceFactor=-1; // flip orientation
      else
        faceFactor=1;
      

      for (i=0; i <= 7; i++) // reset buffers
        for (j=0; j <= 7; j++)
        {
          normal[i][j].x=0;normal[i][j].y=0;normal[i][j].z=0;
          counter[i][j]=0;
        }

      /* process triangles, accumulate normals for Gourad shading */
  
      for (i=0; i <= 6; i++)
        for (j=0; j <= 6; j++) // process block (i,j)
        {
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i,j),r1); // first triangle
          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i,j),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i][j],r3,normal[i][j]);
          counter[i][j]++;

          pDIFFERENCE(NODE(face,i,j+1),NODE(face,i+1,j+1),r1); // second triangle
          pDIFFERENCE(NODE(face,i+1,j),NODE(face,i+1,j+1),r2);
          CROSSPRODUCTp(r1,r2,r3); pMULTIPLY(r3,faceFactor,r3);
          pNORMALIZE(r3);
          pSUM(normal[i+1][j],r3,normal[i+1][j]);
          counter[i+1][j]++;
          pSUM(normal[i][j+1],r3,normal[i][j+1]);
          counter[i][j+1]++;
          pSUM(normal[i+1][j+1],r3,normal[i+1][j+1]);
          counter[i+1][j+1]++;
        }

      
        /* the actual rendering */
        for (j=1; j<=7; j++) 
        {

          if (faceFactor  > 0)
            glFrontFace(GL_CCW); // the usual definition of front face
          else
            glFrontFace(GL_CW); // flip definition of orientation
          glEnable(GL_TEXTURE_2D);
          glBindTexture(GL_TEXTURE_2D, jelloTextureID);
          glBegin(GL_TRIANGLE_STRIP);
          for (i=0; i<=7; i++)
          { 
            float u = i/7.0f; float v = j/7.0f;
            glNormal3f(normal[i][j].x / counter[i][j],normal[i][j].y / counter[i][j], normal[i][j].z / counter[i][j]);
            glTexCoord2f(u, v);
            glVertex3f(NODE(face,i,j).x, NODE(face,i,j).y, NODE(face,i,j).z);
            
            float v2 = (j-1)/7.0f;
            glNormal3f(normal[i][j-1].x / counter[i][j-1],normal[i][j-1].y/ counter[i][j-1], normal[i][j-1].z / counter[i][j-1]);
            glTexCoord2f(u, v2);
            glVertex3f(NODE(face,i,j-1).x, NODE(face,i,j-1).y, NODE(face,i,j-1).z);
          }
          glEnd();
        }
        
        
    }  
  } // end for loop over faces
  glFrontFace(GL_CCW);
  showInclinedPlane(jello);
  showCoordinates();
}

void showBoundingBox()
{
  int i,j;

  glColor4f(0.6,0.6,0.6,0);

  glBegin(GL_LINES);

  // front face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,-2,-2);
    glVertex3f(i,-2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(2,-2,j);
  }

  // back face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(i,2,-2);
    glVertex3f(i,2,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,2,j);
    glVertex3f(2,2,j);
  }

  // left face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(-2,i,-2);
    glVertex3f(-2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(-2,-2,j);
    glVertex3f(-2,2,j);
  }

  // right face
  for(i=-2; i<=2; i++)
  {
    glVertex3f(2,i,-2);
    glVertex3f(2,i,2);
  }
  for(j=-2; j<=2; j++)
  {
    glVertex3f(2,-2,j);
    glVertex3f(2,2,j);
  }
  
  glEnd();

  return;
}

