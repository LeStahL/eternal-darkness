/* File generated with Shader Minifier 1.1.5
 * http://www.ctrl-alt-test.fr
 */
#ifndef SFX_H_
# define SFX_H_
# define VAR_IBLOCKOFFSET "f"
# define VAR_ISAMPLERATE "x"
# define VAR_IVOLUME "v"

const char *sfx_frag =
 "#version 130\n"
 "uniform float f,x,v;"
 "vec2 r(float f)"
 "{"
   "float r=.1,x=1e-05;"
   "return vec2(tanh(sin(acos(-1.)*440.*f)));"
 "}"
 "void main()"
 "{"
   "float g=f+(gl_FragCoord.x-.5+(gl_FragCoord.y-.5)*512.)/x;"
   "vec2 y=v*r(g),a=floor((.5+.5*y)*65536.),d=mod(a,256.)/255.,o=floor(a/256.)/255.;"
   "gl_FragColor=vec4(d.x,o.x,d.y,o.y);"
 "}";

#endif // SFX_H_
