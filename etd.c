/* Eternal Darkness by Team210
 * Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

// #define DEBUG

#include <SDL.h>
#undef main
#include <SDL_opengl.h>

#ifdef DEBUG
#include <stdio.h>
#include <stdlib.h>
#endif

SDL_Window *window;
SDL_GLContext gl_context;
int w,h;

typedef void (*glGetShaderiv_t)(GLuint,  GLenum,  GLint *);
glGetShaderiv_t glGetShaderiv;
typedef void (*glGetShaderInfoLog_t)(GLuint,  GLsizei, GLsizei,  GLchar *);
glGetShaderInfoLog_t glGetShaderInfoLog;
typedef GLuint (CALLBACK *glCreateShader_t)(GLenum);
glCreateShader_t glCreateShader;
typedef GLuint (*glCreateProgram_t)();
glCreateProgram_t glCreateProgram;
typedef void (*glShaderSource_t)(GLuint, GLsizei, GLchar **, GLint *);
glShaderSource_t glShaderSource;
typedef void (*glCompileShader_t)(GLuint);
glCompileShader_t glCompileShader;
typedef void (*glAttachShader_t)(GLuint, GLuint);
glAttachShader_t glAttachShader;
typedef void (*glLinkProgram_t)(GLuint);
glLinkProgram_t glLinkProgram;
typedef void (*glUseProgram_t)(GLuint);
glUseProgram_t glUseProgram;
typedef GLint (*glGetUniformLocation_t)(GLuint, const GLchar *);
glGetUniformLocation_t glGetUniformLocation;
typedef void (*glUniform2f_t)(GLint, GLfloat, GLfloat);
glUniform2f_t glUniform2f;
typedef void (*glUniform1f_t)(GLint, GLfloat);
glUniform1f_t glUniform1f;
typedef void (*glGenFramebuffers_t)(GLsizei, GLuint*);
glGenFramebuffers_t glGenFramebuffers;
typedef void (*glBindFramebuffer_t)(GLenum, GLuint);
glBindFramebuffer_t glBindFramebuffer;
typedef void (*glFramebufferTexture2D_t)(GLenum, GLenum, GLenum, GLuint, GLint);
glFramebufferTexture2D_t glFramebufferTexture2D;
typedef void (*glBlitFramebuffer_t)(GLint, GLint, GLint, GLint, GLint, GLint, GLint, GLint, GLbitfield, GLenum);
glBlitFramebuffer_t glBlitFramebuffer;
typedef void (*glNamedRenderbufferStorage_t) (GLuint, GLenum, GLsizei, GLsizei);
glNamedRenderbufferStorage_t glNamedRenderbufferStorage;

void *glfunc(const char *name)
{
    void *p = (void *)wglGetProcAddress(name);
    if(p == 0 ||
        (p == (void*)0x1) || (p == (void*)0x2) || (p == (void*)0x3) ||
        (p == (void*)-1) )
    {
        HMODULE module = LoadLibrary("opengl32.dll");
#ifdef DEBUG
        printf("OpenGL library address: %p\n", module);
#endif
        p = (void *)GetProcAddress(module, name);
    }
#ifdef DEBUG
    printf("Offset of symbol %s: %p\n", name, p);
#endif
    return p;
}

#ifdef DEBUG
void debug(int shader_handle)
{
    int compile_status = 1;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    printf("Compile status: %s\n", (compile_status == GL_TRUE)?"success":"failed");
    
    if(compile_status != GL_TRUE)
    {
        printf("FAILED.\n");
        int len;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("log length: %d\n", len);
        GLchar *CompileLog = (GLchar*)malloc(len*sizeof(GLchar));
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("error: %s\n", CompileLog);
        free(CompileLog);
    }
}
#endif

void render()
{
    SDL_GL_MakeCurrent(window, gl_context);
        
        glClearColor( 0.2f, 0.4f, 0.1f, 1.0f );
        glClear( GL_COLOR_BUFFER_BIT );
        
        glViewport(0,0,w,h);
        glBegin(GL_QUADS);
        glVertex3f(-1,-1,0);
        glVertex3f(-1,1,0);
        glVertex3f(1,1,0);
        glVertex3f(1,-1,0);
        glEnd();
        
        SDL_GL_SwapWindow(window);
}

void demo()
{
    // init SDL2
    SDL_Init(SDL_INIT_VIDEO | SDL_INIT_EVENTS);
    window = SDL_CreateWindow("title", SDL_WINDOWPOS_UNDEFINED, SDL_WINDOWPOS_UNDEFINED, 1920, 1080, SDL_WINDOW_OPENGL | SDL_WINDOW_FULLSCREEN);
    gl_context = SDL_GL_CreateContext(window);
    
    // gl extensions
    glGetShaderiv = (glGetShaderiv_t)glfunc("glGetShaderiv");
    glGetShaderInfoLog = (glGetShaderInfoLog_t) glfunc("glGetShaderInfoLog");
    glCreateShader = (glCreateShader_t) glfunc("glCreateShader");
    glCreateProgram = (glCreateProgram_t) glfunc("glCreateProgram");
    glShaderSource = (glShaderSource_t) glfunc("glShaderSource");
    glCompileShader = (glCompileShader_t) glfunc("glCompileShader");
    glAttachShader = (glAttachShader_t) glfunc("glAttachShader");
    glLinkProgram = (glLinkProgram_t) glfunc("glLinkProgram");
    glUseProgram = (glUseProgram_t) glfunc("glUseProgram");
    glGetUniformLocation = (glGetUniformLocation_t) glfunc("glGetUniformLocation");
    glUniform2f = (glUniform2f_t) glfunc("glUniform2f");
    glUniform1f = (glUniform1f_t) glfunc("glUniform1f");
    glGenFramebuffers = (glGenFramebuffers_t) glfunc("glGenFramebuffers");
    glBindFramebuffer = (glBindFramebuffer_t) glfunc("glBindFramebuffer");
    glFramebufferTexture2D = (glFramebufferTexture2D_t) glfunc("glFramebufferTexture2D");
    glBlitFramebuffer = (glBlitFramebuffer_t) glfunc("glBlitFramebuffer");
    glNamedRenderbufferStorage = (glNamedRenderbufferStorage_t) glfunc("glNamedRenderbufferStorage");

#ifdef DEBUG //TODO: finish all debug outputs
    printf("Addresses of opengl extensions:\n");
    printf("glGetShaderiv :: %p\n", glGetShaderiv);
#endif
    
    //parameters //TODO: add command line parser here? Or maybe just a config file?
    w = 1920, h = 1080;
    
#ifdef DEBUG
    printf("Team210 proudly presents \"Eternal Darkness\"\nWe are:\nQM :: Code/SFX\nNR4 :: Code/GFX\n");
#endif
    
    // load loading bar shader
#include "load.h"
    int lb_size = strlen(load_frag),
        lb_handle = glCreateShader(GL_FRAGMENT_SHADER),
        lb_program = glCreateProgram();
    glShaderSource(lb_handle, 1, (GLchar **)&load_frag, &lb_size);
    glCompileShader(lb_handle);
#ifdef DEBUG
    debug(lb_handle);
#endif
    glAttachShader(lb_program, lb_handle);
    glLinkProgram(lb_program);
    glUseProgram(lb_program);
    
    int lb_progress_location = glGetUniformLocation(lb_program, VAR_IPROGRESS),
        lb_time_location = glGetUniformLocation(lb_program, VAR_ITIME),
        lb_resolution_location = glGetUniformLocation(lb_program, VAR_IRESOLUTION);
    
    glUniform1f(lb_time_location, 1.2);
    glUniform1f(lb_progress_location, 0.);
    glUniform2f(lb_resolution_location, w, h);
    render();
    
    // load gfx shader
#undef VAR_IRESOLUTION
#undef VAR_ITIME
#include "gfx.h"
    int gfx_size = strlen(gfx_frag),
        gfx_handle = glCreateShader(GL_FRAGMENT_SHADER),
        gfx_program = glCreateProgram();
    glShaderSource(gfx_handle, 1, (GLchar **)&load_frag, &gfx_size);
    glCompileShader(gfx_handle);
#ifdef DEBUG
    debug(lb_handle);
#endif
    glAttachShader(gfx_program, gfx_handle);
    glLinkProgram(gfx_program);
    glUniform1f(lb_progress_location, .25);
    render();
    
    // load sfx shader
    
    
    int running = 1;
    SDL_Event event;
    while(running)
    {
        while(SDL_PollEvent(&event))
        {
            switch( event.type )
            {
                case SDL_KEYDOWN:
                    running = 0;
                    break;
            }
        }
        
        render();
    }
    
    ExitProcess(0);
    
    return 0;
}
