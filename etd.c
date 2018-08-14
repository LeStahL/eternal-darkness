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

// Standard library und CRT rewrite for saving executable size
void *memset(void *ptr, int value, size_t num)
{
    for(int i=num-1; i>=0; i--)
        ((unsigned char *)ptr)[i] = value;
    return ptr;
}

HGLRC glrc;
HDC hdc;

LRESULT CALLBACK WindowProc(HWND hwnd, UINT uMsg, WPARAM wParam, LPARAM lParam)
{
    switch (uMsg)
    {
        case WM_KEYDOWN: 
        case VK_ESCAPE:
            ExitProcess(0);
            return 0;
            break;
            
        case WM_TIMER:
            glClearColor(0.,0.,0.,1.);
            glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
            
            glLoadIdentity();
            
            glColor3f( 1.0, 1., 0.0);
            glPolygonMode(GL_FRONT, GL_FILL);
            
            glBegin(GL_TRIANGLES);
            glVertex3f(-1.,-1.,0.);
            glVertex3f(-1.,1.,0.);
            glVertex3f(1.,1.,0.);
            
            glVertex3f(1.,1.,0.);
            glVertex3f(1.,-1.,0.);
            glVertex3f(-1.,-1.,0.);
            glEnd();
            
            glFlush();
            
            SwapBuffers(hdc);
            
            wglMakeCurrent (hdc, glrc);
            
            break;
    }
    return DefWindowProc(hwnd, uMsg, wParam, lParam);
}

int WINAPI demo(HINSTANCE hInstance, HINSTANCE hPrevInstance, PWSTR pCmdLine, int nCmdShow)
{
    CHAR WindowClass[]  = "Team210 Demo Window";
    
    WNDCLASSEX wc = { 0 };
    wc.cbSize = sizeof(wc);
    wc.style = CS_HREDRAW | CS_VREDRAW | CS_OWNDC;
    wc.lpfnWndProc = &WindowProc;
    wc.cbClsExtra = 0;
    wc.cbWndExtra = 0;
    wc.hInstance = hInstance;
    wc.hIcon = NULL; 
    wc.hCursor = NULL;
    wc.hbrBackground = NULL;
    wc.lpszMenuName = NULL;
    wc.lpszClassName = WindowClass;
    wc.hIconSm = NULL;
    
    RegisterClassEx(&wc);
    
    // Get full screen information
    HMONITOR hmon = MonitorFromWindow(0, MONITOR_DEFAULTTONEAREST);
    MONITORINFO mi = { sizeof(mi) };
    GetMonitorInfo(hmon, &mi);
    
    // Create the window.
    HWND hwnd = CreateWindowEx(
        0,                                                          // Optional window styles.
        WindowClass,                                                // Window class
        "Learn to Program Windows",                                 // Window text
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
    
    hdc = GetDC(hwnd);
    
    int  pf = ChoosePixelFormat(hdc, &pfd); 
    SetPixelFormat(hdc, pf, &pfd);
    
    glrc = wglCreateContext(hdc);
    wglMakeCurrent (hdc, glrc);
    
    // Show it
    ShowWindow(hwnd, TRUE);
    UpdateWindow(hwnd);
    
    // Add timer for regular rendering
    setTimer(hwnd, IDT_TIMER1, 1000./60., (TIMERPROC)NULL );
    
    MSG msg;
    while (GetMessage(&msg, NULL, 0, 0) != 0)
    {
        TranslateMessage(&msg);
        DispatchMessage(&msg);
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
