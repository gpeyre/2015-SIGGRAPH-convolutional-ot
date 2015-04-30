#include "Shader.h"

#include <fstream>
#include <iostream>

void Shader::loadVertex(const char* filename)
{
   load(GL_VERTEX_SHADER, filename, vertexShader);
}

void Shader::loadFragment(const char* filename)
{
   load(GL_FRAGMENT_SHADER, filename, fragmentShader);
}

void Shader::loadGeometry(const char* filename)
{
#ifdef GL_GEOMETRY_SHADER_EXT
   load(GL_GEOMETRY_SHADER_EXT, filename, geometryShader);
#else
   std::cerr << "Error: geometry shaders not supported!" << std::endl;
#endif
}

void Shader::enable()
{
   if (!linked)
   {
      glLinkProgram(program);
      linked = true;
   }
   glUseProgram(program);
}

void Shader::disable() const
{
   glUseProgram(0);
}

Shader::operator GLuint() const
{
   return program;
}

void Shader::load(GLenum shaderType, const char* filename, GLuint& shader)
{
   std::string source;
   if (!readSource(filename, source)) return;

   if (program == 0)
   {
      program = glCreateProgram();
   }

   if (shader != 0)
   {
      glDetachShader(program, shader);
   }

   shader = glCreateShader(shaderType);
   const char* source_c_str = source.c_str();
   glShaderSource(shader, 1, &(source_c_str), NULL);

   glCompileShader(shader);
   GLint compileStatus;
   glGetShaderiv(shader, GL_COMPILE_STATUS, &compileStatus);

   if (compileStatus == GL_TRUE)
   {
      glAttachShader(program, shader);
      linked = false;
   }
   else
   {
      GLsizei maxLength = 0;
      glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &maxLength);

      if (maxLength > 0)
      {
         GLchar* infoLog = new char[maxLength];
         GLsizei length;

         glGetShaderInfoLog(shader, maxLength, &length, infoLog);
         delete[] infoLog;
      }
   }
}

bool Shader::readSource(const char* filename, std::string& source)
{
   source = "";

   std::ifstream in(filename);
   if (!in.is_open())
   {
      std::cerr << "Error: could not open shader file ";
      std::cerr << filename;
      std::cerr << " for input!" << std::endl;
      return false;
   }

   std::string line;
   while (getline(in, line))
   {
      source += line;
   }
   return true;
}
