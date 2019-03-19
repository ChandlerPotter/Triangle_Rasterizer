/*
  Chandler Poter CIS 441 Winter 2019.
  This program is a triangle rasterizer and renders pixel correct images.
  It uses techniques such as scanline algorithm, linear interpolation, z-buffer, and phong shading.
  It takes a VKT geometry file as input and produces a .png image file as output.
*/
#include <iostream>
#include <vtkDataSet.h>
#include <vtkImageData.h>
#include <vtkPNGWriter.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataReader.h>
#include <vtkPoints.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <math.h>
#include <stdio.h> 

using std::cerr;
using std::endl;

double ceil_fun(double f)
{
    return ceil(f-0.00001);
}

double floor_fun(double f)
{
    return floor(f+0.00001);
}

struct LightingParameters
{
    LightingParameters(void)
    {
         lightDir[0] = -0.6;
         lightDir[1] = 0;
         lightDir[2] = -0.8;
         Ka = 0.3;
         Kd = 0.7;
         Ks = 2.3;
         alpha = 2.5;
    };
  

    double lightDir[3]; // The direction of the light source
    double Ka;           // The coefficient for ambient lighting.
    double Kd;           // The coefficient for  lighting.
    double Ks;           // The coefficient for specular lighting.
    double alpha;        // The exponent term for specular lighting.
};

LightingParameters lp;

class Matrix
{
  public:
    double          A[4][4];

    void            TransformPoint(const double *ptIn, double *ptOut);
    static Matrix   ComposeMatrices(const Matrix &, const Matrix &);
    void            Print(ostream &o);
};

void Matrix::Print(ostream &o)
{
    for (int i = 0 ; i < 4 ; i++)
    {
        char str[256];
        sprintf(str, "(%.7f %.7f %.7f %.7f)\n", A[i][0], A[i][1], A[i][2], A[i][3]);
        o << str;
    }
}

Matrix Matrix::ComposeMatrices(const Matrix &M1, const Matrix &M2)
{
    Matrix rv;
    for (int i = 0 ; i < 4 ; i++)
        for (int j = 0 ; j < 4 ; j++)
        {
            rv.A[i][j] = 0;
            for (int k = 0 ; k < 4 ; k++)
                rv.A[i][j] += M1.A[i][k]*M2.A[k][j];
        }

    return rv;
}

void Matrix::TransformPoint(const double *ptIn, double *ptOut)
{
    ptOut[0] = ptIn[0]*A[0][0]
             + ptIn[1]*A[1][0]
             + ptIn[2]*A[2][0]
             + ptIn[3]*A[3][0];
    ptOut[1] = ptIn[0]*A[0][1]
             + ptIn[1]*A[1][1]
             + ptIn[2]*A[2][1]
             + ptIn[3]*A[3][1];
    ptOut[2] = ptIn[0]*A[0][2]
             + ptIn[1]*A[1][2]
             + ptIn[2]*A[2][2]
             + ptIn[3]*A[3][2];
    ptOut[3] = ptIn[0]*A[0][3]
             + ptIn[1]*A[1][3]
             + ptIn[2]*A[2][3]
             + ptIn[3]*A[3][3];
}

class Camera
{
  public:
    double          near, far;
    double          angle;
    double          position[3];
    double          focus[3];
    double          up[3];
};


double SineParameterize(int curFrame, int nFrames, int ramp)
{
    int nNonRamp = nFrames-2*ramp;
    double height = 1./(nNonRamp + 4*ramp/M_PI);
    if (curFrame < ramp)
    {
        double factor = 2*height*ramp/M_PI;
        double eval = cos(M_PI/2*((double)curFrame)/ramp);
        return (1.-eval)*factor;
    }
    else if (curFrame > nFrames-ramp)
    {
        int amount_left = nFrames-curFrame;
        double factor = 2*height*ramp/M_PI;
        double eval =cos(M_PI/2*((double)amount_left/ramp));
        return 1. - (1-eval)*factor;
    }
    double amount_in_quad = ((double)curFrame-ramp);
    double quad_part = amount_in_quad*height;
    double curve_part = height*(2*ramp)/M_PI;
    return quad_part+curve_part;
}

Camera GetCamera(int frame, int nframes)
{
    double t = SineParameterize(frame, nframes, nframes/10);
    Camera c;
    c.near = 5;
    c.far = 200;
    c.angle = M_PI/6;
    c.position[0] = 40*sin(2*M_PI*t);
    c.position[1] = 40*cos(2*M_PI*t);
    c.position[2] = 40;
    c.focus[0] = 0;
    c.focus[1] = 0;
    c.focus[2] = 0;
    c.up[0] = 0;
    c.up[1] = 1;
    c.up[2] = 0;
    return c;
}



vtkImageData *NewImage(int width, int height)
{
    vtkImageData *img = vtkImageData::New();
    img->SetDimensions(width, height, 1);
    img->AllocateScalars(VTK_UNSIGNED_CHAR, 3);

    return img;
}

void WriteImage(vtkImageData *img, const char *filename)
{
   std::string full_filename = filename;
   full_filename += ".png";
   vtkPNGWriter *writer = vtkPNGWriter::New();
   writer->SetInputData(img);
   writer->SetFileName(full_filename.c_str());
   writer->Write();
   writer->Delete();
}

class Triangle
{
  public:
      double         X[3];
      double         Y[3];
      double         Z[3];
      double         colors[3][3];
	    double 		     normals[3][3];
      double         shade[3];
  // would some methods for transforming the triangle in place be helpful?
};

class Screen
{
  public:
      unsigned char   *buffer;
      double          *zee_buffer;
      int width, height;
      int *triangleID;

      void SetPixel(double r, double c, double color, int ind);
      void Set_Z_Pixel(int r, int c, double z, int red, int green, int blue, int ind, int Zind, int triangleID);
      void Print();
  // would some methods for accessing and setting pixels be helpful?
};

void Screen::SetPixel(double r, double c, double color, int ind)
{
  if(r < 0 || r >= this->height || c < 0 || c >= this->width)
  {
    return;
  }
  else
  {
    buffer[ind] = color;
  }
}
/*
  Set_Z_Pixel() puts the color values into the buffer, drawing them onto that pixel. It checks the screen bounds and the 
  Z depth value to decide if the color should be drawn. Does the same as SetPixel but checks the zeebuffer.
*/
void Screen::Set_Z_Pixel(int r, int c, double z, int red, int green, int blue, int ind, int Zind, int triangleID)
{
  if(r < 0 || r >= this->height || c < 0 || c >= this->width || z <= this->zee_buffer[Zind])
  {
    return;
  }
  else
  {
  this->triangleID[Zind] = triangleID;
 	this->zee_buffer[Zind] = z;
	buffer[ind] = red;
	buffer[ind+1] = green;
	buffer[ind+2] = blue;
  }
}

/*
  A print method for debugging.
*/
void Screen::Print()
{
  int cnt = 0;
  for (int i = 0; i < this->height; i++)
  {
    for (int j = 0; j < this->width; j++)
    {
      printf("Index: %d (col=%d, row=%d) comes from triangle %d and has color %d, %d, %d and z-value %f\n", cnt, j, i, this->triangleID[((i * this->width) + j)], this->buffer[((i * this->width) + j) * 3], this->buffer[(((i * this->width) + j) * 3)+1], buffer[(((i * this->width) + j) * 3)+2], this->zee_buffer[((i * this->width) + j)]);
      cnt += 1;
    }
  }
}

#define NORMALS
std::vector<Triangle>
GetTriangles(void)
{
    vtkPolyDataReader *rdr = vtkPolyDataReader::New();
    rdr->SetFileName("proj1e_geometry.vtk");
    cerr << "Reading" << endl;
    rdr->Update();
    cerr << "Done reading" << endl;
    if (rdr->GetOutput()->GetNumberOfCells() == 0)
    {
        cerr << "Unable to open file!!" << endl;
        exit(EXIT_FAILURE);
    }
    vtkPolyData *pd = rdr->GetOutput();

    int numTris = pd->GetNumberOfCells();
    vtkPoints *pts = pd->GetPoints();
    vtkCellArray *cells = pd->GetPolys();
    vtkDoubleArray *var = (vtkDoubleArray *) pd->GetPointData()->GetArray("hardyglobal");
    double *color_ptr = var->GetPointer(0);
    //vtkFloatArray *var = (vtkFloatArray *) pd->GetPointData()->GetArray("hardyglobal");
    //float *color_ptr = var->GetPointer(0);
    vtkFloatArray *n = (vtkFloatArray *) pd->GetPointData()->GetNormals();
    float *normals = n->GetPointer(0);
    std::vector<Triangle> tris(numTris);
    vtkIdType npts;
    vtkIdType *ptIds;
    int idx;
    for (idx = 0, cells->InitTraversal() ; cells->GetNextCell(npts, ptIds) ; idx++)
    {
        if (npts != 3)
        {
            cerr << "Non-triangles!! ???" << endl;
            exit(EXIT_FAILURE);
        }
        double *pt = NULL;
        pt = pts->GetPoint(ptIds[0]);
        tris[idx].X[0] = pt[0];
        tris[idx].Y[0] = pt[1];
        tris[idx].Z[0] = pt[2];
#ifdef NORMALS
        tris[idx].normals[0][0] = normals[3*ptIds[0]+0];
        tris[idx].normals[0][1] = normals[3*ptIds[0]+1];
        tris[idx].normals[0][2] = normals[3*ptIds[0]+2];
#endif
        pt = pts->GetPoint(ptIds[1]);
        tris[idx].X[1] = pt[0];
        tris[idx].Y[1] = pt[1];
        tris[idx].Z[1] = pt[2];
#ifdef NORMALS
        tris[idx].normals[1][0] = normals[3*ptIds[1]+0];
        tris[idx].normals[1][1] = normals[3*ptIds[1]+1];
        tris[idx].normals[1][2] = normals[3*ptIds[1]+2];
#endif
        pt = pts->GetPoint(ptIds[2]);
        tris[idx].X[2] = pt[0];
        tris[idx].Y[2] = pt[1];
        tris[idx].Z[2] = pt[2];
#ifdef NORMALS
        tris[idx].normals[2][0] = normals[3*ptIds[2]+0];
        tris[idx].normals[2][1] = normals[3*ptIds[2]+1];
        tris[idx].normals[2][2] = normals[3*ptIds[2]+2];
#endif

        // 1->2 interpolate between light blue, dark blue
        // 2->2.5 interpolate between dark blue, cyan
        // 2.5->3 interpolate between cyan, green
        // 3->3.5 interpolate between green, yellow
        // 3.5->4 interpolate between yellow, orange
        // 4->5 interpolate between orange, brick
        // 5->6 interpolate between brick, salmon
        double mins[7] = { 1, 2, 2.5, 3, 3.5, 4, 5 };
        double maxs[7] = { 2, 2.5, 3, 3.5, 4, 5, 6 };
        unsigned char RGB[8][3] = { { 71, 71, 219 },
                                    { 0, 0, 91 },
                                    { 0, 255, 255 },
                                    { 0, 128, 0 },
                                    { 255, 255, 0 },
                                    { 255, 96, 0 },
                                    { 107, 0, 0 },
                                    { 224, 76, 76 }
                                  };
        for (int j = 0 ; j < 3 ; j++)
        {
            float val = color_ptr[ptIds[j]];
            int r;
            for (r = 0 ; r < 7 ; r++)
            {
                if (mins[r] <= val && val < maxs[r])
                    break;
            }
            if (r == 7)
            {
                cerr << "Could not interpolate color for " << val << endl;
                exit(EXIT_FAILURE);
            }
            double proportion = (val-mins[r]) / (maxs[r]-mins[r]);
            tris[idx].colors[j][0] = (RGB[r][0]+proportion*(RGB[r+1][0]-RGB[r][0]))/255.0;
            tris[idx].colors[j][1] = (RGB[r][1]+proportion*(RGB[r+1][1]-RGB[r][1]))/255.0;
            tris[idx].colors[j][2] = (RGB[r][2]+proportion*(RGB[r+1][2]-RGB[r][2]))/255.0;
        }
    }

    return tris;
}



bool double_equals(double a, double b, double epsilon = 0.0001)
{
    return fabs(a - b) < epsilon;
}

/*
  For lineal interpolation, a function to find the t value of 3 given doubles.
*/
double t_value(double x, double a, double b)
{
  return ((x-a)/(b-a));
}

/*
  This is the linear interpolation function and it uses the t value returned from the t_value function.
*/
double LERP(double t, double a, double b)
{
  return (a + (t * (b - a)));
}

/*
  A Rasterization routine for going down trianges. Uses scanline to draw pixels. 
*/
void RasterizeGoingDownTriangle(Triangle T, Screen screen, int k)
{
  double minY, topLeftX, topRightX;
  double minX, topLeftY, topRightY;
  double minZ, topLeftZ, topRightZ;
  double min_red, topLeft_red, topRight_red;
  double min_blue, topLeft_blue, topRight_blue;
  double min_green, topLeft_green, topRight_green;
  double t, fieldVal_z_left, fieldVal_red_left, fieldVal_green_left, fieldVal_blue_left;
  double fieldVal_z_right, fieldVal_red_right, fieldVal_green_right, fieldVal_blue_right;
  double minShade, topLeftShade, topRightShade;


  minY = std::min(T.Y[0], T.Y[1]);              //Finds the minimun Y value using the std::min method.
  minY = std::min(minY, T.Y[2]);

  if (minY == T.Y[0])                           //Checks if the first Y value is the minimum Y value (lowest Y point on a graph).
  {
    minX = T.X[0];                              //if that is the case, assigns all corresponding X, Z, and color values.
    minZ = T.Z[0];
    min_red = T.colors[0][0];
    min_green = T.colors[0][1];
    min_blue = T.colors[0][2];
    minShade = T.shade[0];
    if(T.X[1] < T.X[2])                         //Compares the other 2 X values to find which is the right or left, and assigns corresponding Y, Z, and color values.
    {
      topLeftX = T.X[1];
      topRightX = T.X[2];
      topLeftY = T.Y[1];
      topRightY = T.Y[2];
      topLeftZ = T.Z[1];
      topRightZ = T.Z[2];
      topLeft_red = T.colors[1][0];
      topLeft_green = T.colors[1][1];
      topLeft_blue = T.colors[1][2];
      topRight_red = T.colors[2][0];
      topRight_green = T.colors[2][1];
      topRight_blue = T.colors[2][2];
      topLeftShade = T.shade[1];
      topRightShade = T.shade[2];
    }
    else                                        //Assigns the correct X and Y values if X[2] is > X[1].
    {
      topLeftX = T.X[2];
      topRightX = T.X[1];
      topLeftY = T.Y[2];
      topRightY = T.Y[1];
      topLeftZ = T.Z[2];
      topRightZ = T.Z[1];
      topLeft_red = T.colors[2][0];
      topLeft_green = T.colors[2][1];
      topLeft_blue = T.colors[2][2];
      topRight_red = T.colors[1][0];
      topRight_green = T.colors[1][1];
      topRight_blue = T.colors[1][2];
      topLeftShade = T.shade[2];
      topRightShade = T.shade[1];
    }
  }
  if (minY == T.Y[1])                         //same as previos min check but for the next Y value and the corresponding if checks do the same as before.
  {
    minX = T.X[1];
    minZ = T.Z[1];
    min_red = T.colors[1][0];
    min_green = T.colors[1][1];
    min_blue = T.colors[1][2];
    minShade = T.shade[1];
    if(T.X[0] < T.X[2])
    {
      topLeftX = T.X[0];
      topRightX = T.X[2];
      topLeftY = T.Y[0];
      topRightY = T.Y[2];
      topLeftZ = T.Z[0];
      topRightZ = T.Z[2];
      topLeft_red = T.colors[0][0];
      topLeft_green = T.colors[0][1];
      topLeft_blue = T.colors[0][2];
      topRight_red = T.colors[2][0];
      topRight_green = T.colors[2][1];
      topRight_blue = T.colors[2][2];
      topLeftShade = T.shade[0];
      topRightShade = T.shade[2];
    }
    else
    {
      topLeftX = T.X[2];
      topRightX = T.X[0];
      topLeftY = T.Y[2];
      topRightY = T.Y[0];
      topLeftZ = T.Z[2];
      topRightZ = T.Z[0];
      topLeft_red = T.colors[2][0];
      topLeft_green = T.colors[2][1];
      topLeft_blue = T.colors[2][2];
      topRight_red = T.colors[0][0];
      topRight_green = T.colors[0][1];
      topRight_blue = T.colors[0][2];
      topLeftShade = T.shade[2];
      topRightShade = T.shade[0];
    }
  }
  if (minY == T.Y[2])
  {
    minX = T.X[2];
    minZ = T.Z[2];
    min_red = T.colors[2][0];
    min_green = T.colors[2][1];
    min_blue = T.colors[2][2];
    minShade = T.shade[2];
    if(T.X[1] < T.X[0])
    {
      topLeftX = T.X[1];
      topRightX = T.X[0];
      topLeftY = T.Y[1];
      topRightY = T.Y[0];
      topLeftZ = T.Z[1];
      topRightZ = T.Z[0];
      topLeft_red = T.colors[1][0];
      topLeft_green = T.colors[1][1];
      topLeft_blue = T.colors[1][2];
      topRight_red = T.colors[0][0];
      topRight_green = T.colors[0][1];
      topRight_blue = T.colors[0][2];
      topLeftShade = T.shade[1];
      topRightShade = T.shade[0];
    }
    else
    {
      topLeftX = T.X[0];
      topRightX = T.X[1];
      topLeftY = T.Y[0];
      topRightY = T.Y[1];
      topLeftZ = T.Z[0];
      topRightZ = T.Z[1];
      topLeft_red = T.colors[0][0];
      topLeft_green = T.colors[0][1];
      topLeft_blue = T.colors[0][2];
      topRight_red = T.colors[1][0];
      topRight_green = T.colors[1][1];
      topRight_blue = T.colors[1][2];
      topLeftShade = T.shade[0];
      topRightShade = T.shade[1];
    }
  }

  double mLeft, mRight; //y = mx + b
  double bLeft, bRight;

  mLeft = ((topLeftY - minY)/(topLeftX - minX));                  //Slope of the left side of the triangle.
  mRight = ((topRightY - minY)/(topRightX - minX));               //Slope of the right side of the triangle.
  bLeft = minY - (mLeft * minX);                                  //y intercept of left and right sides.
  bRight = minY - (mRight * minX);

  double row = ceil_fun(minY);
  double maxRow = floor_fun(topLeftY);

  for(double minRow = row; minRow <= maxRow; minRow++)
  {
    double xLeft;
    double xRight;
    double t_left, t_right;

    xLeft = ((minRow - bLeft)/mLeft);
    xRight = ((minRow - bRight)/mRight);

    if (topLeftX == minX)
    {
      xLeft = topLeftX;
    }
    if (topRightX == minX)
    {
      xRight = topRightX;
    }
    t_left = t_value(xLeft, minX, topLeftX);
    t_right = t_value(xRight, minX, topRightX);

    if(xRight == minX)
    {
      t_right = t_value(minRow, minY, topRightY);
    }
    if(xLeft == minX)
    {
      t_left = t_value(minRow, minY, topLeftY);
    }

    fieldVal_z_left = LERP(t_left, minZ, topLeftZ);
    fieldVal_z_right = LERP(t_right, minZ, topRightZ);

    double fieldVal_shade_left = LERP(t_left, minShade, topLeftShade);
    double fieldVal_shade_right = LERP(t_right, minShade, topRightShade);

    fieldVal_red_left = LERP(t_left, min_red, topLeft_red);
    fieldVal_green_left = LERP(t_left, min_green, topLeft_green);
    fieldVal_blue_left = LERP(t_left, min_blue, topLeft_blue);

    fieldVal_red_right = LERP(t_right, min_red, topRight_red);
    fieldVal_green_right = LERP(t_right, min_green, topRight_green);
    fieldVal_blue_right = LERP(t_right, min_blue, topRight_blue);

    for(double column = ceil_fun(xLeft); column <= floor_fun(xRight); column++)
    {
      double fieldVal_z, fieldVal_red, fieldVal_green, fieldVal_blue, fieldVal_shade;
      unsigned char color_red, color_green, color_blue;
      t = t_value(column, xLeft, xRight);
      fieldVal_z = LERP(t, fieldVal_z_left, fieldVal_z_right);
      fieldVal_shade = LERP(t, fieldVal_shade_left, fieldVal_shade_right);
      fieldVal_red = LERP(t, fieldVal_red_left, fieldVal_red_right);
      fieldVal_green = LERP(t, fieldVal_green_left, fieldVal_green_right);
      fieldVal_blue = LERP(t, fieldVal_blue_left, fieldVal_blue_right);

      double red_temp = (fieldVal_red * fieldVal_shade);
      double green_temp = (fieldVal_green * fieldVal_shade);
      double blue_temp = (fieldVal_blue * fieldVal_shade);

      color_red = ceil_fun(red_temp * 255);
      color_green = ceil_fun(green_temp * 255);
      color_blue = ceil_fun(blue_temp * 255);

      if(red_temp > 1)
      {
        color_red = 255;
      }
      if(green_temp > 1)
      {
        color_green = 255;
      }
      if(blue_temp > 1)
      {
        color_blue = 255;
      }

      int index = ((minRow * screen.width) + column) * 3;
      int Zindex = ((minRow * screen.width) + column);
      screen.Set_Z_Pixel(minRow, column, fieldVal_z, color_red, color_green, color_blue, index, Zindex, k);
    }
  }
}

void RasterizeGoingUpTriangle(Triangle T, Screen screen, int k)
{
  double maxX, bottomLeftX, bottomRightX;
  double maxY, bottomLeftY, bottomRightY;
  double maxZ, bottomLeftZ, bottomRightZ;
  double max_red, bottomLeft_red, bottomRight_red;
  double max_blue, bottomLeft_blue, bottomRight_blue;
  double max_green, bottomLeft_green, bottomRight_green;
  double t, t_left, t_right, fieldVal_z_left, fieldVal_red_left, fieldVal_green_left, fieldVal_blue_left;
  double fieldVal_z_right, fieldVal_red_right, fieldVal_green_right, fieldVal_blue_right;
  double maxShade, bottomLeftShade, bottomRightShade;

  maxY = std::max(T.Y[0], T.Y[1]);
  maxY = std::max(maxY, T.Y[2]);

  if (maxY == T.Y[0])
  {
    maxX = T.X[0];
    maxZ = T.Z[0];
    max_red = T.colors[0][0];
    max_green = T.colors[0][1];
    max_blue = T.colors[0][2];
    maxShade = T.shade[0];
    if(T.X[1] < T.X[2])
    {
      bottomLeftX = T.X[1];
      bottomRightX = T.X[2];
      bottomLeftY = T.Y[1];
      bottomRightY = T.Y[2];
      bottomLeftZ = T.Z[1];
      bottomRightZ = T.Z[2];
      bottomLeft_red = T.colors[1][0];
      bottomLeft_green = T.colors[1][1];
      bottomLeft_blue = T.colors[1][2];
      bottomRight_red = T.colors[2][0];
      bottomRight_green = T.colors[2][1];
      bottomRight_blue = T.colors[2][2];
      bottomLeftShade = T.shade[1];
      bottomRightShade = T.shade[2];

    }
    else
    {
      bottomLeftX = T.X[2];
      bottomRightX = T.X[1];
      bottomLeftY = T.Y[2];
      bottomRightY = T.Y[1];
      bottomLeftZ = T.Z[2];
      bottomRightZ = T.Z[1];
      bottomLeft_red = T.colors[2][0];
      bottomLeft_green = T.colors[2][1];
      bottomLeft_blue = T.colors[2][2];
      bottomRight_red = T.colors[1][0];
      bottomRight_green = T.colors[1][1];
      bottomRight_blue = T.colors[1][2];
      bottomLeftShade = T.shade[2];
      bottomRightShade = T.shade[1];

    }
  }
  if (maxY == T.Y[1])
  {
    maxX = T.X[1];
    maxZ = T.Z[1];
    max_red = T.colors[1][0];
    max_green = T.colors[1][1];
    max_blue = T.colors[1][2];
    maxShade = T.shade[1];
    if(T.X[0] < T.X[2])
    {
      bottomLeftX = T.X[0];
      bottomRightX = T.X[2];
      bottomLeftY = T.Y[0];
      bottomRightY = T.Y[2];
      bottomLeftZ = T.Z[0];
      bottomRightZ = T.Z[2];
      bottomLeft_red = T.colors[0][0];
      bottomLeft_green = T.colors[0][1];
      bottomLeft_blue = T.colors[0][2];
      bottomRight_red = T.colors[2][0];
      bottomRight_green = T.colors[2][1];
      bottomRight_blue = T.colors[2][2];
      bottomLeftShade = T.shade[0];
      bottomRightShade = T.shade[2];
    }
    else
    {
      bottomLeftX = T.X[2];
      bottomRightX = T.X[0];
      bottomLeftY = T.Y[2];
      bottomRightY = T.Y[0];
      bottomLeftZ = T.Z[2];
      bottomRightZ = T.Z[0];
      bottomLeft_red = T.colors[2][0];
      bottomLeft_green = T.colors[2][1];
      bottomLeft_blue = T.colors[2][2];
      bottomRight_red = T.colors[0][0];
      bottomRight_green = T.colors[0][1];
      bottomRight_blue = T.colors[0][2];
      bottomLeftShade = T.shade[2];
      bottomRightShade = T.shade[0];
    }
  }
  if (maxY == T.Y[2])
  {
    maxX = T.X[2];
    maxZ = T.Z[2];
    max_red = T.colors[2][0];
    max_green = T.colors[2][1];
    max_blue = T.colors[2][2];
    maxShade = T.shade[2];
    if(T.X[1] < T.X[0])
    {
      bottomLeftX = T.X[1];
      bottomRightX = T.X[0];
      bottomLeftY = T.Y[1];
      bottomRightY = T.Y[0];
      bottomLeftZ = T.Z[1];
      bottomRightZ = T.Z[0];
      bottomLeft_red = T.colors[1][0];
      bottomLeft_green = T.colors[1][1];
      bottomLeft_blue = T.colors[1][2];
      bottomRight_red = T.colors[0][0];
      bottomRight_green = T.colors[0][1];
      bottomRight_blue = T.colors[0][2];
      bottomLeftShade = T.shade[1];
      bottomRightShade = T.shade[0];
    }
    else
    {
      bottomLeftX = T.X[0];
      bottomRightX = T.X[1];
      bottomLeftY = T.Y[0];
      bottomRightY = T.Y[1];
      bottomLeftZ = T.Z[0];
      bottomRightZ = T.Z[1];
      bottomLeft_red = T.colors[0][0];
      bottomLeft_green = T.colors[0][1];
      bottomLeft_blue = T.colors[0][2];
      bottomRight_red = T.colors[1][0];
      bottomRight_green = T.colors[1][1];
      bottomRight_blue = T.colors[1][2];
      bottomLeftShade = T.shade[0];
      bottomRightShade = T.shade[1];
    }
  }

  double mLeft, mRight; //y = mx + b
  double bLeft, bRight;

  mLeft = ((maxY - bottomLeftY)/(maxX - bottomLeftX));
  mRight = ((maxY - bottomRightY)/(maxX - bottomRightX));
  bLeft = maxY - (mLeft * maxX);
  bRight = maxY - (mRight * maxX);

  double row = ceil_fun(bottomLeftY);
  double maxRow = floor_fun(maxY);

  for(double minRow = row; minRow <= maxRow; minRow++)
  {
    double xLeft;
    double xRight;

    xLeft = ((minRow - bLeft)/mLeft);
    xRight = ((minRow - bRight)/mRight);

    if (bottomLeftX == maxX)
    {
      xLeft = bottomLeftX;
    }
    if (bottomRightX == maxX)
    {
      xRight = bottomRightX;
    }

    t_left = t_value(xLeft, maxX, bottomLeftX);
    t_right = t_value(xRight, maxX, bottomRightX);

    if(xRight == maxX)
    {
      t_right = t_value(minRow, maxY, bottomRightY);
    }
    if(xLeft == maxX)
    {
      t_left = t_value(minRow, maxY, bottomLeftY);
    }

    fieldVal_z_left = LERP(t_left, maxZ, bottomLeftZ);
    fieldVal_z_right = LERP(t_right, maxZ, bottomRightZ);

    double fieldVal_shade_left = LERP(t_left, maxShade, bottomLeftShade);
    double fieldVal_shade_right = LERP(t_right, maxShade, bottomRightShade);

    fieldVal_red_left = LERP(t_left, max_red, bottomLeft_red);
    fieldVal_green_left = LERP(t_left, max_green, bottomLeft_green);
    fieldVal_blue_left = LERP(t_left, max_blue, bottomLeft_blue);

    fieldVal_red_right = LERP(t_right, max_red, bottomRight_red);
    fieldVal_green_right = LERP(t_right, max_green, bottomRight_green);
    fieldVal_blue_right = LERP(t_right, max_blue, bottomRight_blue);

    for(double column = ceil_fun(xLeft); column <= floor_fun(xRight); column++)
    {
      double fieldVal_z, fieldVal_red, fieldVal_green, fieldVal_blue, fieldVal_shade;
      unsigned char color_red, color_green, color_blue;
      t = t_value(column, xLeft, xRight);
      fieldVal_z = LERP(t, fieldVal_z_left, fieldVal_z_right);
      fieldVal_shade = LERP(t, fieldVal_shade_left, fieldVal_shade_right);
      fieldVal_red = LERP(t, fieldVal_red_left, fieldVal_red_right);
      fieldVal_green = LERP(t, fieldVal_green_left, fieldVal_green_right);
      fieldVal_blue = LERP(t, fieldVal_blue_left, fieldVal_blue_right);

      double red_temp = (fieldVal_red * fieldVal_shade);
      double green_temp = (fieldVal_green * fieldVal_shade);
      double blue_temp = (fieldVal_blue * fieldVal_shade);

      color_red = ceil_fun(red_temp * 255);
      color_green = ceil_fun(green_temp * 255);
      color_blue = ceil_fun(blue_temp * 255);

      if(red_temp > 1)
      {
        color_red = 255;
      }
      if(green_temp > 1)
      {
        color_green = 255;
      }
      if(blue_temp > 1)
      {
        color_blue = 255;
      }

      int index = ((minRow * screen.width) + column) * 3;
      int Zindex = (minRow *screen.width) + column;
      screen.Set_Z_Pixel(minRow, column, fieldVal_z, color_red, color_green, color_blue, index, Zindex, k);
    }
  }
}

void RasterizeArbitraryTriangle(Triangle triangle, Screen screen, int k)
{
  double mSlope, bVal, numberX;
  double minY, maxY, middleY, min_red, max_red, middle_red;
  double minX, maxX, middleX, min_green, max_green, middle_green;
  double minZ, maxZ, middleZ, min_blue, max_blue, middle_blue;
  double t, fieldVal_z;
  double fieldVal_red, fieldVal_green, fieldVal_blue;
  double minShade, middleShade, maxShade;

  minY = std::min(triangle.Y[0], triangle.Y[1]);
  minY = std::min(minY, triangle.Y[2]);

  maxY = std::max(triangle.Y[0], triangle.Y[1]);
  maxY = std::max(maxY, triangle.Y[2]);

  if(minY != triangle.Y[0] && maxY != triangle.Y[0])
  {
     middleY = triangle.Y[0];
  }
  if(minY != triangle.Y[1] && maxY != triangle.Y[1])
  {
     middleY = triangle.Y[1];
  }
  if(minY != triangle.Y[2] && maxY != triangle.Y[2])
  {
     middleY = triangle.Y[2];
  }

  if (minY == triangle.Y[0])
  {
    minX = triangle.X[0];
    minZ = triangle.Z[0];
    min_red = triangle.colors[0][0];
    min_green = triangle.colors[0][1];
    min_blue = triangle.colors[0][2];
    minShade = triangle.shade[0];
  }
  if (minY == triangle.Y[1])
  {
    minX = triangle.X[1];
    minZ = triangle.Z[1];
    min_red = triangle.colors[1][0];
    min_green = triangle.colors[1][1];
    min_blue = triangle.colors[1][2];
    minShade = triangle.shade[1];
  }
  if (minY == triangle.Y[2])
  {
    minX = triangle.X[2];
    minZ = triangle.Z[2];
    min_red = triangle.colors[2][0];
    min_green = triangle.colors[2][1];
    min_blue = triangle.colors[2][2];
    minShade = triangle.shade[2];
  }

  if (maxY == triangle.Y[0])
  {
    maxX = triangle.X[0];
    maxZ = triangle.Z[0];
    max_red = triangle.colors[0][0];
    max_green = triangle.colors[0][1];
    max_blue = triangle.colors[0][2];
    maxShade = triangle.shade[0];
  }
  if (maxY == triangle.Y[1])
  {
    maxX = triangle.X[1];
    maxZ = triangle.Z[1];
    max_red = triangle.colors[1][0];
    max_green = triangle.colors[1][1];
    max_blue = triangle.colors[1][2];
    maxShade = triangle.shade[1];
  }
  if (maxY == triangle.Y[2])
  {
    maxX = triangle.X[2];
    maxZ = triangle.Z[2];
    max_red = triangle.colors[2][0];
    max_green = triangle.colors[2][1];
    max_blue = triangle.colors[2][2];
    maxShade = triangle.shade[2];
  }

  if (middleY == triangle.Y[0])
  {
    middleX = triangle.X[0];
    middleZ = triangle.Z[0];
    middle_red = triangle.colors[0][0];
    middle_green = triangle.colors[0][1];
    middle_blue = triangle.colors[0][2];
    middleShade = triangle.shade[0];
  }
  if (middleY == triangle.Y[1])
  {
    middleX = triangle.X[1];
    middleZ = triangle.Z[1];
    middle_red = triangle.colors[1][0];
    middle_green = triangle.colors[1][1];
    middle_blue = triangle.colors[1][2];
    middleShade = triangle.shade[1];
  }
  if (middleY == triangle.Y[2])
  {
    middleX = triangle.X[2];
    middleZ = triangle.Z[2];
    middle_red = triangle.colors[2][0];
    middle_green = triangle.colors[2][1];
    middle_blue = triangle.colors[2][2];
    middleShade = triangle.shade[2];
  }

  mSlope = (maxY-minY)/(maxX-minX);
  bVal = minY - (mSlope * minX);
  numberX = (middleY - bVal)/mSlope; // (numberX, middleY)

  if(maxX == minX)
  {
    numberX = maxX;
  }

  Triangle UP;
  Triangle DOWN;

  t = t_value(middleY, maxY, minY);
  fieldVal_z = LERP(t, maxZ, minZ);

  UP.X[0] = maxX;
  UP.Y[0] = maxY;
  UP.Z[0] = maxZ;
  UP.shade[0] = maxShade;

  UP.X[2] = middleX;
  UP.Y[2] = middleY;
  UP.Z[2] = middleZ;
  UP.shade[2] = middleShade;

  UP.X[1] = numberX;
  UP.Y[1] = middleY;
  UP.Z[1] = fieldVal_z;
  UP.shade[1] = LERP(t, maxShade, minShade);

  UP.colors[0][0] = max_red;
  UP.colors[0][1] = max_green;
  UP.colors[0][2] = max_blue;

  UP.colors[2][0] = middle_red;
  UP.colors[2][1] = middle_green;
  UP.colors[2][2] = middle_blue;

  fieldVal_red = LERP(t, max_red, min_red);
  UP.colors[1][0] = fieldVal_red;

  fieldVal_green = LERP(t, max_green, min_green);
  UP.colors[1][1] = fieldVal_green;

  fieldVal_blue = LERP(t, max_blue, min_blue);
  UP.colors[1][2] = fieldVal_blue;

  DOWN.X[0] = minX;
  DOWN.Y[0] = minY;
  DOWN.Z[0] = minZ;
  DOWN.shade[0] = minShade;

  DOWN.X[1] = numberX;
  DOWN.Y[1] = middleY;
  DOWN.Z[1] = fieldVal_z;
  DOWN.shade[1] = LERP(t, maxShade, minShade);

  DOWN.X[2] = middleX;
  DOWN.Y[2] = middleY;
  DOWN.Z[2] = middleZ;
  DOWN.shade[2] = middleShade;

  DOWN.colors[0][0] = min_red;
  DOWN.colors[0][1] = min_green;
  DOWN.colors[0][2] = min_blue;

  DOWN.colors[1][0] = fieldVal_red;
  DOWN.colors[1][1] = fieldVal_green;
  DOWN.colors[1][2] = fieldVal_blue;

  DOWN.colors[2][0] = middle_red;
  DOWN.colors[2][1] = middle_green;
  DOWN.colors[2][2] = middle_blue;


  RasterizeGoingUpTriangle(UP, screen, k);
  RasterizeGoingDownTriangle(DOWN, screen, k);
}


void MainRasterize(Triangle t, Screen screen, int k)
{
  if(double_equals(t.Y[0],t.Y[1]))
  {
    if (t.Y[2] < t.Y[0] || t.Y[2] < t.Y[1])
    {
      RasterizeGoingDownTriangle(t, screen, k);
    }
    else
    {
      RasterizeGoingUpTriangle(t, screen, k);
    }
  }

  else if(double_equals(t.Y[0],t.Y[2]))
  {
    if (t.Y[1] < t.Y[0] || t.Y[1] < t.Y[2])
    {
      RasterizeGoingDownTriangle(t, screen, k);
    }
    else
    {
      RasterizeGoingUpTriangle(t, screen, k);
    }
  }

  else if(double_equals(t.Y[1],t.Y[2]))
  {
    if (t.Y[0] < t.Y[1] || t.Y[0] < t.Y[2])
    {
      RasterizeGoingDownTriangle(t, screen, k);
    }
    else
    {
      RasterizeGoingUpTriangle(t, screen, k);
    }
  }
  else
  {
    RasterizeArbitraryTriangle(t, screen, k); 
  }
}

Triangle ShadeTriangle(Triangle T, Camera C)
{
  double tempDot1, diffuse, specular, cosVal;
  double R[3], Rnorm, viewDir[3], viewDirNorm, shading[3];
  
  for(int i = 0; i < 3; i++)
  {
    tempDot1 = (lp.lightDir[0] * T.normals[i][0]) + (lp.lightDir[1] * T.normals[i][1]) + (lp.lightDir[2] * T.normals[i][2]);
    
    R[0] = ((2 * tempDot1 * T.normals[i][0]) - lp.lightDir[0]);
    R[1] = ((2 * tempDot1 * T.normals[i][1]) - lp.lightDir[1]);
    R[2] = ((2 * tempDot1 * T.normals[i][2]) - lp.lightDir[2]);

    Rnorm = sqrt((R[0]*R[0]) + (R[1]*R[1]) + (R[2]*R[2]));
    R[0] = R[0]/Rnorm;
    R[1] = R[1]/Rnorm;
    R[2] = R[2]/Rnorm;

    viewDir[0] = C.position[0] - T.X[i];
    viewDir[1] = C.position[1] - T.Y[i];
    viewDir[2] = C.position[2] - T.Z[i];

    viewDirNorm = sqrt((viewDir[0]*viewDir[0]) + (viewDir[1]*viewDir[1]) + (viewDir[2]*viewDir[2]));

    viewDir[0] = viewDir[0]/viewDirNorm;
    viewDir[1] = viewDir[1]/viewDirNorm;
    viewDir[2] = viewDir[2]/viewDirNorm;

    cosVal = ((viewDir[0] * R[0]) + (viewDir[1] * R[1]) + (viewDir[2] * R[2]));
    diffuse = lp.Kd * abs(tempDot1);
    specular = lp.Ks * (pow(cosVal, lp.alpha));

    if(cosVal < 0)
    {
      specular = 0;
    }

    shading[i] = lp.Ka + diffuse + specular;
    
    T.shade[i] = shading[i];
  }

  return T;
}


Matrix CameraTransform(Camera c)
{
  double o[3], w[3], v[3], u[3];

  for(int i = 0; i < 3; i++)
  {
    v[i] = c.up[i];
    o[i] = c.position[i];
    w[i] = o[i] - c.focus[i];
  }
  u[0] = (c.up[1]*w[2] - c.up[2]*w[1]);
  u[1] = (c.up[2]*w[0] - c.up[0]*w[2]);
  u[2] = (c.up[0]*w[1] - c.up[1]*w[0]);

  v[0] = (w[1]*u[2] - w[2]*u[1]);
  v[1] = (w[2]*u[0] - w[0]*u[2]);
  v[2] = (w[0]*u[1] - w[1]*u[0]);

  double vSquare = sqrt((v[0]*v[0]) + (v[1]*v[1]) + (v[2]*v[2]));
  double uSquare = sqrt((u[0]*u[0]) + (u[1]*u[1]) + (u[2]*u[2]));
  double wSquare = sqrt((w[0]*w[0]) + (w[1]*w[1]) + (w[2]*w[2]));

  for(int j = 0; j < 3; j++)
  {
    v[j] = v[j]/vSquare;
    w[j] = w[j]/wSquare;
    u[j] = u[j]/uSquare;
  }

  Matrix m;

  double t[3];
  t[0] = 0 - o[0];
  t[1] = 0 - o[1];
  t[2] = 0 - o[2];

  m.A[0][0] = u[0];
  m.A[0][1] = v[0];
  m.A[0][2] = w[0];
  m.A[0][3] = 0;

  m.A[1][0] = u[1];
  m.A[1][1] = v[1];
  m.A[1][2] = w[1];
  m.A[1][3] = 0;

  m.A[2][0] = u[2];
  m.A[2][1] = v[2];
  m.A[2][2] = w[2];
  m.A[2][3] = 0;

  m.A[3][0] = ((u[0]*t[0]) + (u[1]*t[1]) + (u[2]*t[2]));
  m.A[3][1] = ((v[0]*t[0]) + (v[1]*t[1]) + (v[2]*t[2]));
  m.A[3][2] = ((w[0]*t[0]) + (w[1]*t[1]) + (w[2]*t[2]));
  m.A[3][3] = 1;

    return m;
}

Matrix ViewTransform(Camera c)
{
  Matrix m;

  m.A[0][0] = 1/(tan(c.angle/2));
  m.A[0][1] = 0;
  m.A[0][2] = 0;
  m.A[0][3] = 0;

  m.A[1][0] = 0;
  m.A[1][1] = 1/(tan(c.angle/2));
  m.A[1][2] = 0;
  m.A[1][3] = 0;

  m.A[2][0] = 0;
  m.A[2][1] = 0;
  m.A[2][2] = ((c.far + c.near)/(c.far - c.near));
  m.A[2][3] = -1;

  m.A[3][0] = 0;
  m.A[3][1] = 0;
  m.A[3][2] = (2*c.far*c.near)/(c.far - c.near);
  m.A[3][3] = 0;

  return m;
}

Matrix DeviceTransform(Camera c, Screen screen)
{
  Matrix m;

  m.A[0][0] = (screen.width)/2;
  m.A[0][1] = 0;
  m.A[0][2] = 0;
  m.A[0][3] = 0;

  m.A[1][0] = 0;
  m.A[1][1] = (screen.height)/2;
  m.A[1][2] = 0;
  m.A[1][3] = 0;

  m.A[2][0] = 0;
  m.A[2][1] = 0;
  m.A[2][2] = 1;
  m.A[2][3] = 0;

  m.A[3][0] = (screen.width)/2;
  m.A[3][1] = (screen.height)/2;
  m.A[3][2] = 0;
  m.A[3][3] = 1;

  return m;
}

Triangle TransformTriangles(Matrix M, Triangle T, Camera C)
{
  double vertex1[4], vertex2[4], vertex3[4];

  vertex1[0] = T.X[0];
  vertex1[1] = T.Y[0];
  vertex1[2] = T.Z[0];
  vertex1[3] = 1;

  vertex2[0] = T.X[1];
  vertex2[1] = T.Y[1];
  vertex2[2] = T.Z[1];
  vertex2[3] = 1;

  vertex3[0] = T.X[2];
  vertex3[1] = T.Y[2];
  vertex3[2] = T.Z[2];
  vertex3[3] = 1;

  vertex1[0] = (T.X[0]*M.A[0][0] + T.Y[0]*M.A[1][0] + T.Z[0]*M.A[2][0] + 1*M.A[3][0]);
  vertex1[1] = (T.X[0]*M.A[0][1] + T.Y[0]*M.A[1][1] + T.Z[0]*M.A[2][1] + 1*M.A[3][1]);
  vertex1[2] = (T.X[0]*M.A[0][2] + T.Y[0]*M.A[1][2] + T.Z[0]*M.A[2][2] + 1*M.A[3][2]);
  vertex1[3] = (T.X[0]*M.A[0][3] + T.Y[0]*M.A[1][3] + T.Z[0]*M.A[2][3] + 1*M.A[3][3]); //needs to = 1

  if(vertex1[3] != 1)
  {
    vertex1[0] = vertex1[0]/vertex1[3];
    vertex1[1] = vertex1[1]/vertex1[3];
    vertex1[2] = vertex1[2]/vertex1[3];
  }

  vertex2[0] = (T.X[1]*M.A[0][0] + T.Y[1]*M.A[1][0] + T.Z[1]*M.A[2][0] + 1*M.A[3][0]);
  vertex2[1] = (T.X[1]*M.A[0][1] + T.Y[1]*M.A[1][1] + T.Z[1]*M.A[2][1] + 1*M.A[3][1]);
  vertex2[2] = (T.X[1]*M.A[0][2] + T.Y[1]*M.A[1][2] + T.Z[1]*M.A[2][2] + 1*M.A[3][2]);
  vertex2[3] = (T.X[1]*M.A[0][3] + T.Y[1]*M.A[1][3] + T.Z[1]*M.A[2][3] + 1*M.A[3][3]); //needs to = 1

  if(vertex2[3] != 1)
  {
    vertex2[0] = vertex2[0]/vertex2[3];
    vertex2[1] = vertex2[1]/vertex2[3];
    vertex2[2] = vertex2[2]/vertex2[3];
  }

  vertex3[0] = (T.X[2]*M.A[0][0] + T.Y[2]*M.A[1][0] + T.Z[2]*M.A[2][0] + 1*M.A[3][0]);
  vertex3[1] = (T.X[2]*M.A[0][1] + T.Y[2]*M.A[1][1] + T.Z[2]*M.A[2][1] + 1*M.A[3][1]);
  vertex3[2] = (T.X[2]*M.A[0][2] + T.Y[2]*M.A[1][2] + T.Z[2]*M.A[2][2] + 1*M.A[3][2]);
  vertex3[3] = (T.X[2]*M.A[0][3] + T.Y[2]*M.A[1][3] + T.Z[2]*M.A[2][3] + 1*M.A[3][3]); //needs to = 1

  if(vertex3[3] != 1)
  {
    vertex3[0] = vertex3[0]/vertex3[3];
    vertex3[1] = vertex3[1]/vertex3[3];
    vertex3[2] = vertex3[2]/vertex3[3];
  }

Triangle temp;
temp.X[0] = vertex1[0];
temp.Y[0] = vertex1[1];
temp.Z[0] = vertex1[2];

temp.X[1] = vertex2[0];
temp.Y[1] = vertex2[1];
temp.Z[1] = vertex2[2];

temp.X[2] = vertex3[0];
temp.Y[2] = vertex3[1];
temp.Z[2] = vertex3[2];

temp.colors[0][0] = T.colors[0][0];
temp.colors[0][1] = T.colors[0][1];
temp.colors[0][2] = T.colors[0][2];
temp.colors[1][0] = T.colors[1][0];
temp.colors[1][1] = T.colors[1][1];
temp.colors[1][2] = T.colors[1][2];
temp.colors[2][0] = T.colors[2][0];
temp.colors[2][1] = T.colors[2][1];
temp.colors[2][2] = T.colors[2][2];

T = ShadeTriangle(T, C);
temp.shade[0] = T.shade[0];
temp.shade[1] = T.shade[1];
temp.shade[2] = T.shade[2];

return temp;
}

int main()
{
   vtkImageData *image = NewImage(1000, 1000);
   unsigned char *buffer =
     (unsigned char *) image->GetScalarPointer(0,0,0);


   int npixels = 1000*1000;
   for (int i = 0 ; i < npixels*3 ; i++)
       buffer[i] = 0;

   std::vector<Triangle> triangles = GetTriangles();

   Screen screen;
   screen.buffer = buffer;
   screen.width = 1000;
   screen.height = 1000;
   screen.zee_buffer = new double[screen.width * screen.height];
   for (int i = 0; i < (screen.width * screen.height); i++)
   {
     screen.zee_buffer[i] = -1;
   }

   screen.triangleID = new int[1000*1000];
   for (int i = 0 ; i < npixels ; i++)
       screen.triangleID[i] = -1;


   // YOUR CODE GOES HERE TO DEPOSIT THE COLORS FROM TRIANGLES
   // INTO PIXELS USING THE SCANLINE ALGORITHM

   for (int i = 0; i < 250; i+=250)
   {
     int npixels = 1000*1000;
     for (int i = 0 ; i < npixels*3 ; i++)
         buffer[i] = 0;
     screen.buffer = buffer;

     for (int i = 0; i < npixels; i++)
     {
       screen.zee_buffer[i] = -1;
     }

     Camera cam = GetCamera(i,1000);

     Matrix CamT = CameraTransform(cam);
     Matrix ViewT = ViewTransform(cam);
     Matrix DevT = DeviceTransform(cam, screen);

     Matrix tempM = CamT.ComposeMatrices(CamT, ViewT);
     Matrix TotalM = tempM.ComposeMatrices(tempM, DevT);

     std::vector<Triangle> tempTriangles = GetTriangles();

     for(int tri = 0; tri < triangles.size(); tri++)
     {
       //printf("triangle %d\n", tri);
       tempTriangles[tri] = TransformTriangles(TotalM, triangles[tri], cam);
       //ShadeTriangle(triangles[tri], cam);
       MainRasterize(tempTriangles[tri], screen, tri);
     }

      WriteImage(image, "frame000");
   }
}