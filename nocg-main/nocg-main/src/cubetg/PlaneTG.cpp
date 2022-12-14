#include <glad/glad.h>
#define GLFW_INCLUDE_NONE
#include <GLFW/glfw3.h>

#include <glm/common.hpp>
#include <glm/matrix.hpp>
#include <glm/gtx/transform.hpp>

#include <math.h>
#include <iostream>

#include "Utils.h"

#include "PlaneTG.h"

using namespace nocg;
using glm::vec3;
using glm::mat4;
using std::cout;
using std::endl;

PlaneTG::PlaneTG(const vec2& dims)
    :_dims(dims)
{
    // load program 
    vector<string> shaderFiles = { "planeTG.vert", "planeTG.frag", "planeTG.geom"};
    _program = loadShaders(shaderFiles);

    // create geometry 
    _createGeometry();

    // load texture
    _textureID = loadTexture("text.png", true);
}

PlaneTG::~PlaneTG()
{

}

// hit test 
void PlaneTG::hitTest(double x, double y, 
    const glm::mat4& vMat, const glm::mat4& pMat, const glm::vec4& viewport)
{
    // unproject mouse (x, y) and find starting and ending points 
    // or ray passing through the 3D model
    mat4 mvMat = mvMat * _modelMat;
    vec3 p1 = glm::unProject(vec3(x, y, 0.0f), mvMat, pMat, viewport);
    vec3 p2 = glm::unProject(vec3(x, y, 1.0f), mvMat, pMat, viewport);

    GLuint ssbo;
    glGenBuffers(1, &ssbo);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, ssbo);
    float color[] = {1.f, 1.f, 0.f};
    glBufferData(GL_SHADER_STORAGE_BUFFER, 3*sizeof(float), color, GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, ssbo);
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);

    // render
    render(vMat, pMat);

    //GLuint buffer; // handle to buffer
    //std::vector<float> storage(n); // n is the size  
    //glGetNamedBufferSubData(buffer, 0, n * sizeof(float), storage.data());
    glGetNamedBufferSubData(ssbo, 0, 3 * sizeof(float), color);

    cout << color[0] << "," << color[1] << "," << color[2] << endl;
}

void PlaneTG::render(const glm::mat4& vMat, const glm::mat4& pMat)
{
    glUseProgram(_program);

    // set model matrix
    GLint mMatLoc = glGetUniformLocation(_program, "mMat");
    glUniformMatrix4fv(mMatLoc, 1, GL_FALSE, &_modelMat[0][0]);


    // set view matrix
    GLint vMatLoc = glGetUniformLocation(_program, "vMat");
    glUniformMatrix4fv(vMatLoc, 1, GL_FALSE, &vMat[0][0]);

    // set projection matrix 
    GLint pMatLoc = glGetUniformLocation(_program, "pMat");
    glUniformMatrix4fv(pMatLoc, 1, GL_FALSE, &pMat[0][0]);

    // texture settings
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _textureID);
    GLint samplerLoc = glGetUniformLocation(_program, "sampler");
    glUniform1i(samplerLoc, 0);

    glBindVertexArray(_vao);

    // draw 
    glDrawArrays(GL_LINE_STRIP_ADJACENCY, 0, _vertexCount);

    glBindVertexArray(0);
    glUseProgram(0);
    glBindTexture(GL_TEXTURE_2D, 0);
}

void PlaneTG::_createPlaneTG()
{
    // PlaneTG is centered on XY planeTG at Z=0
    _vertices = {   
        -_dims.x, _dims.y, 0.0f,
        _dims.x, _dims.y, 0.0f,
        -_dims.x, -_dims.y, 0.0f,
        _dims.x, -_dims.y, 0.0f,
    };
    _normals = { 
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f,
        0.0f, 0.0f, 1.0f
    };
    _texCoords = {
        0.0f, 1.0f,
        1.0f, 1.0f,
        0.0f, 0.0f,
        1.0f, 0.0f
    };
}

void PlaneTG::_createGeometry()
{
    glCreateVertexArrays(1, &_vao);
    glBindVertexArray(_vao);

    _createPlaneTG();
    _vertexCount = _vertices.size() / 3;

    GLuint buffer[4];

    // create buffers
    glCreateBuffers(4, buffer);

    // vertices 

    // Initialize the first buffer
    glNamedBufferStorage(buffer[0], sizeof(float) * _vertices.size(), _vertices.data(), 0);
    // Bind it to the vertex array - offset zero, stride = sizeof(vec3)
    glVertexArrayVertexBuffer(_vao, 0, buffer[0], 0, 3 * sizeof(float));
    // Tell OpenGL what the format of the attribute is
    glVertexArrayAttribFormat(_vao, 0, 3, GL_FLOAT, GL_FALSE, 0);
    // Tell OpenGL which vertex buffer binding to use for this attribute
    glVertexArrayAttribBinding(_vao, 0, 0);

    // Enable the attribute
    glEnableVertexArrayAttrib(_vao, 0);

    // normals 

    // Initialize the first buffer
    glNamedBufferStorage(buffer[1], sizeof(float) * _normals.size(), _normals.data(), 0);
    // Bind it to the vertex array - offset zero, stride = sizeof(vec3)
    glVertexArrayVertexBuffer(_vao, 1, buffer[1], 0, 3 * sizeof(float));
    // Tell OpenGL what the format of the attribute is
    glVertexArrayAttribFormat(_vao, 1, 3, GL_FLOAT, GL_FALSE, 0);
    // Tell OpenGL which vertex buffer binding to use for this attribute
    glVertexArrayAttribBinding(_vao, 1, 1);

    // Enable the attribute
    glEnableVertexArrayAttrib(_vao, 1);

    // tex coords 

    // Initialize the first buffer
    glNamedBufferStorage(buffer[2], sizeof(float) * _texCoords.size(), _texCoords.data(), 0);
    // Bind it to the vertex array - offset zero, stride = sizeof(vec3)
    glVertexArrayVertexBuffer(_vao, 2, buffer[2], 0, 2 * sizeof(float));
    // Tell OpenGL what the format of the attribute is
    glVertexArrayAttribFormat(_vao, 2, 2, GL_FLOAT, GL_FALSE, 0);
    // Tell OpenGL which vertex buffer binding to use for this attribute
    glVertexArrayAttribBinding(_vao, 2, 2);

    // Enable the attribute
    glEnableVertexArrayAttrib(_vao, 2);

    // tangents 
#if 0
    // Initialize the first buffer
    glNamedBufferStorage(buffer[3], sizeof(float) * _tangents.size(), _tangents.data(), 0);
    // Bind it to the vertex array - offset zero, stride = sizeof(vec3)
    glVertexArrayVertexBuffer(_vao, 3, buffer[3], 0, 3 * sizeof(float));
    // Tell OpenGL what the format of the attribute is
    glVertexArrayAttribFormat(_vao, 3, 3, GL_FLOAT, GL_FALSE, 0);
    // Tell OpenGL which vertex buffer binding to use for this attribute
    glVertexArrayAttribBinding(_vao, 3, 3);

    // Enable the attribute
    glEnableVertexArrayAttrib(_vao, 3);
#endif

    glBindVertexArray(0);

    // clear memory
    clear();
}

