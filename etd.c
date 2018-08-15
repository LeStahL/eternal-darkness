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

int _fltused = 0;

#define WIN32_LEAN_AND_MEAN
#define VC_EXTRALEAN
#include <windows.h>

#include <GL/gl.h>
#include "glext.h"

//TODO: remove
#include <stdio.h>

// Standard library and CRT rewrite for saving executable size
void *memset(void *ptr, int value, size_t num)
{
    for(int i=num-1; i>=0; i--)
        ((unsigned char *)ptr)[i] = value;
    return ptr;
}

size_t strlen(const char *str)
{
    int len = 0;
    while(str[len] != '\0') ++len;
    return len;
}

// Windows globals
// HGLRC glrc;
// HDC hdc;

// OpenGL extensions
typedef void (*glGetShaderiv_t)(GLuint,  GLenum,  GLint *);
glGetShaderiv_t glGetShaderiv;
typedef void (*glGetShaderInfoLog_t)(GLuint,  GLsizei, GLsizei,  GLchar *);
glGetShaderInfoLog_t glGetShaderInfoLog;
typedef GLuint (*glCreateShader_t)(GLenum);
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

// Shader globals
int lb_program, lb_progress_location, lb_resolution_location, lb_time_location;

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    printf("wndproc, msg=%d\n", uMsg);
    
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

//TODO: remove
void debug(int shader_handle)
{
    printf("Debugging information for shader handle %d:\n", shader_handle);
    int compile_status = 0;
    glGetShaderiv(shader_handle, GL_COMPILE_STATUS, &compile_status);
    if(compile_status != GL_TRUE)
    {
        printf("FAILED.\n");
        int len = 12;
        glGetShaderiv(shader_handle, GL_INFO_LOG_LENGTH, &len);
        printf("log length: %d\n", len);
        GLchar *CompileLog = (GLchar *)malloc(len*sizeof(GLchar));
        glGetShaderInfoLog(shader_handle, len, NULL, CompileLog);
        printf("error: %s\n", CompileLog);
        free(CompileLog);
    } 
}

int WINAPI demo(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
{
    // TODO: remove 
    AllocConsole();
    freopen("CONIN$", "r", stdin);
    freopen("CONOUT$", "w", stdout);
    freopen("CONOUT$", "w", stderr);
    
    CHAR WindowClass[]  = "Team210 Demo Window";
    
    WNDCLASSEX wc = { 0 };
    wc.cbSize = sizeof(wc);
    wc.style = CS_OWNDC | CS_VREDRAW | CS_HREDRAW;
    wc.lpfnWndProc = &WindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
    wc.hIcon = LoadIcon(NULL, IDI_WINLOGO); 
    wc.hCursor = LoadCursor(NULL, IDC_ARROW);
    wc.hbrBackground = NULL;
    wc.lpszMenuName = NULL;
    wc.lpszClassName = WindowClass;
    wc.hIconSm = NULL;
    
    printf("register class: %d\n", RegisterClassEx(&wc));
    
    // Get full screen information
    HMONITOR hmon = MonitorFromWindow(0, MONITOR_DEFAULTTONEAREST);
    MONITORINFO mi = { sizeof(mi) };
    GetMonitorInfo(hmon, &mi);
    
    // Create the window.
    HWND hwnd = CreateWindowEx(
        0,                                                          // Optional window styles.
        WindowClass,                                                // Window class
        ":: NR4^QM/Team210 :: GO - MAKE A DEMO ::",                                 // Window text
        WS_POPUP | WS_VISIBLE,                                      // Window style
        mi.rcMonitor.left,
        mi.rcMonitor.top,
        mi.rcMonitor.right - mi.rcMonitor.left,
        mi.rcMonitor.bottom - mi.rcMonitor.top,                     // Size and position
        
        NULL,                                                       // Parent window    
        NULL,                                                       // Menu
        hInstance,                                                  // Instance handle
        0                                                           // Additional application data
    );
    
    // Show it
    ShowWindow(hwnd, TRUE);
    UpdateWindow(hwnd);
    
    // Create OpenGL context
    PIXELFORMATDESCRIPTOR pfd =
    {
        sizeof(PIXELFORMATDESCRIPTOR),
        1,
        PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,    //Flags
        PFD_TYPE_RGBA,        // The kind of framebuffer. RGBA or palette.
        32,                   // Colordepth of the framebuffer.
        0, 0, 0, 0, 0, 0,
        0,
        0,
        0,
        0, 0, 0, 0,
        24,                   // Number of bits for the depthbuffer
        8,                    // Number of bits for the stencilbuffer
        0,                    // Number of Aux buffers in the framebuffer.
        PFD_MAIN_PLANE,
        0,
        0, 0, 0
    };
    
    HDC hdc = GetDC(hwnd);
    
    int  pf = ChoosePixelFormat(hdc, &pfd); 
    SetPixelFormat(hdc, pf, &pfd);
    
    HGLRC glrc = wglCreateContext(hdc);
    wglMakeCurrent (hdc, glrc);
    
    // OpenGL extensions
    glGetShaderiv = (glGetShaderiv_t) wglGetProcAddress("glGetShaderiv");
    glGetShaderInfoLog = (glGetShaderInfoLog_t) wglGetProcAddress("glGetShaderInfoLog");
    glCreateShader = (glCreateShader_t) wglGetProcAddress("glCreateShader");
    glCreateProgram = (glCreateProgram_t) wglGetProcAddress("glCreateProgram");
    glShaderSource = (glShaderSource_t) wglGetProcAddress("glShaderSource");
    glCompileShader = (glCompileShader_t) wglGetProcAddress("glCompileShader");
    glAttachShader = (glAttachShader_t) wglGetProcAddress("glAttachShader");
    glLinkProgram = (glLinkProgram_t) wglGetProcAddress("glLinkProgram");
    glUseProgram = (glUseProgram_t) wglGetProcAddress("glUseProgram");
    glGetUniformLocation = (glGetUniformLocation_t) wglGetProcAddress("glGetUniformLocation");
    glUniform2f = (glUniform2f_t) wglGetProcAddress("glUniform2f");
    glUniform1f = (glUniform1f_t) wglGetProcAddress("glUniform1f");
    glGenFramebuffers = (glGenFramebuffers_t) wglGetProcAddress("glGenFramebuffers");
    glBindFramebuffer = (glBindFramebuffer_t) wglGetProcAddress("glBindFramebuffer");
    glFramebufferTexture2D = (glFramebufferTexture2D_t) wglGetProcAddress("glFramebufferTexture2D");
    glBlitFramebuffer = (glBlitFramebuffer_t) wglGetProcAddress("glBlitFramebuffer");
    glNamedRenderbufferStorage = (glNamedRenderbufferStorage_t) wglGetProcAddress("glNamedRenderbufferStorage");
    
    // TODO: add a way of configuring this in the future.
    int w = 1920, h = 1080;
    // Init loading bar.
    #include "load.h"
    int lb_size = strlen(load_frag);
    printf("lb size is: %d\n", lb_size);
    int lb_handle = glCreateShader(GL_FRAGMENT_SHADER);
    printf("lb handle is: %d\n", lb_handle);
    lb_program = glCreateProgram();
    printf("lb program handle is: %d\n", lb_program);
    glShaderSource(lb_handle, 1, (GLchar **)&load_frag, &lb_size);
    glCompileShader(lb_handle);
    //     debug(lb_handle);
    glAttachShader(lb_program, lb_handle);
    printf("shader attached to program.\n");
    glLinkProgram(lb_program);
    printf("linked program.\n");
    glUseProgram(lb_program);
    printf("used program.\n");
    lb_progress_location = glGetUniformLocation(lb_program, VAR_IPROGRESS);
    lb_time_location = glGetUniformLocation(lb_program, VAR_ITIME);
    lb_resolution_location = glGetUniformLocation(lb_program, VAR_IRESOLUTION);
    printf("have uniform locations:\n -> progress: %d\n -> time: %d\n --> resolution: %d\n", lb_progress_location, lb_time_location, lb_resolution_location);
    glUniform1f(lb_time_location, 1.2);
    glUniform1f(lb_progress_location, 0.);
    glUniform2f(lb_resolution_location, w, h);
    printf("set uniform values\n");
    printf("FINISHED WM_CREATE\n");
    UINT_PTR t = SetTimer(hwnd, 1, 1000, NULL);
  
    // Main loop
    MSG msg;
    while(GetMessage(&msg, hwnd, 0, 0) > 0 )
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
        
        UINT m = msg.message;
        printf("msg: %d\n",m);
        printf("keydown: %d\n", WM_KEYDOWN);
        switch(m)
        {
            case WM_KEYDOWN:
                ExitProcess(0);
                break;
                
            case WM_TIMER:
                //HDC hdc = GetDC(hwnd);
                printf("WM_TIMER\n");
                printf("rendering.\n");
                
                glUseProgram(lb_program);
                printf("used program.\n");
                //             
                glClearColor(0.,0.,0.,1.);
                glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
                printf("cleared.\n");
                glLoadIdentity();
                
                glColor3f( 1.0, 1., 0.0);
                glPolygonMode(GL_FRONT, GL_FILL);
                printf("starting to draw.\n");
                glBegin(GL_TRIANGLES);
                printf("began GL\n");
                glVertex3f(-1.,-1.,0.);
                glVertex3f(-1.,1.,0.);
                glVertex3f(1.,1.,0.);
                printf("first triangle\n");
                glVertex3f(1.,1.,0.);
                glVertex3f(1.,-1.,0.);
                glVertex3f(-1.,-1.,0.);
                glEnd();
                printf("now finsihed drawig vertices\n");
                glFlush();
                
                SwapBuffers(hdc);
                printf("making current.\n");
                //                         wglMakeCurrent (hdc, glrc);
                printf("rendering done.\n");
        }
        
    }
    
    return 0;
}

// void initdemo()
// {

// 
// #ifdef DEBUG //TODO: finish all debug outputs
//     printf("Addresses of opengl extensions:\n");
//     printf("glGetShaderiv :: %p\n", glGetShaderiv);
// #endif
//     
//     //parameters //TODO: add command line parser here? Or maybe just a config file?
//     w = 1920, h = 1080;
//     
// #ifdef DEBUG
//     printf("Team210 proudly presents \"Eternal Darkness\"\nWe are:\nQM :: Code/SFX\nNR4 :: Code/GFX\n");
// #endif
//     
//     // load loading bar shader
// #include "load.h"
//     int lb_size = strlen(load_frag),
//         lb_handle = glCreateShader(GL_FRAGMENT_SHADER),
//         lb_program = glCreateProgram();
//     glShaderSource(lb_handle, 1, (GLchar **)&load_frag, &lb_size);
//     glCompileShader(lb_handle);
// #ifdef DEBUG
//     printf("Compile Debug info for load.frag:");
//     debug(lb_handle);
// #endif
//     glAttachShader(lb_program, lb_handle);
//     glLinkProgram(lb_program);
// #ifdef DEBUG
//     printf("Link Debug info for load.frag:");
//     debug_link(lb_program);
// #endif
//     glUseProgram(lb_program);
//     int lb_progress_location = glGetUniformLocation(lb_program, VAR_IPROGRESS),
//         lb_time_location = glGetUniformLocation(lb_program, VAR_ITIME),
//         lb_resolution_location = glGetUniformLocation(lb_program, VAR_IRESOLUTION);
//     glUniform1f(lb_time_location, 1.2);
//     glUniform1f(lb_progress_location, 0.);
//     glUniform2f(lb_resolution_location, w, h);
//     render();
//     
//     // load gfx shader
// #undef VAR_IRESOLUTION
// #undef VAR_ITIME
// #include "gfx.h"
//     int gfx_size = strlen(gfx_frag),
//         gfx_handle = glCreateShader(GL_FRAGMENT_SHADER),
//         gfx_program = glCreateProgram();
//     glShaderSource(gfx_handle, 1, (GLchar **)&gfx_frag, &gfx_size);
//     glCompileShader(gfx_handle);
// #ifdef DEBUG
//     printf("Debug info for gfx.frag:");
//     debug(gfx_handle);
// #endif
//     glAttachShader(gfx_program, gfx_handle);
//     glLinkProgram(gfx_program);
// #ifdef DEBUG
//     printf("Link Debug info for gfx.frag:");
//     debug_link(gfx_program);
// #endif
//     int gfx_time_location =  glGetUniformLocation(gfx_program, VAR_ITIME),
//         gfx_resolution_location = glGetUniformLocation(gfx_program, VAR_IRESOLUTION);
//     glUniform1f(lb_progress_location, .25);
//     
//     render();
//     
//     // load sfx shader
// #include "sfx.h"
//     int sfx_size = strlen(sfx_frag),
//         sfx_handle = glCreateShader(GL_FRAGMENT_SHADER),
//         sfx_program = glCreateProgram();
//     glShaderSource(sfx_handle, 1, (GLchar **)&sfx_frag, &sfx_size);
//     glCompileShader(sfx_handle);
// #ifdef DEBUG
//     printf("Debug info for sfx.frag:");
//     debug(sfx_handle);
// #endif
//     glAttachShader(sfx_program, sfx_handle);
// //     glLinkProgram(sfx_program); //TODO: Sound shader has array indexing in for loop; this crashes the linker. FIXME
// #ifdef DEBUG
//     printf("Link Debug info for sfx.frag:");
//     debug_link(sfx_program);
// #endif
//     
//     glUniform1f(lb_progress_location, .5);
//     render();
//     
//     glUseProgram(gfx_program);
//     glUniform1f(gfx_time_location, 10.);
//     glUniform2f(gfx_resolution_location, w, h);
//     render();
// }
