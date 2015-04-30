#ifndef SHADER_H
#define SHADER_H

#ifdef __APPLE__
#include <GLUT/glut.h>
#include <OpenGL/gl.h>
#include <OpenGL/glext.h>
#else
#include <GL/gl.h>
#include <GL/glext.h>
#include <GL/glut.h>
#endif

#include <string>

class Shader
{
public:
   Shader() 
   : vertexShader(0)
   , fragmentShader(0)
   , geometryShader(0)
   , program(0)
   , linked(false)
   { }


   ~Shader()
   {
      if (program) glDeleteProgram(program);
      if ( vertexShader) glDeleteShader(  vertexShader);
      if (fragmentShader) glDeleteShader(fragmentShader);
      if (geometryShader) glDeleteShader(geometryShader);
   }

   void loadVertex(const char* filename);

   void loadFragment(const char* filename);

   void loadGeometry(const char* filename);

   void enable();

   void disable() const;

   operator GLuint() const;

protected:
   void load(GLenum shaderType, const char* filename, GLuint& shader);
   bool readSource(const char* filename, std::string& source);

   GLuint vertexShader;
   GLuint fragmentShader;
   GLuint geometryShader;
   GLuint program;
   bool linked;
};

#endif // SHADER_H
