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
   "return vec2(sin(2764.56*f)*exp(-3.*f));"
 "}"
 "void main()"
 "{"
   "float g=f+(gl_FragCoord.x-.5+(gl_FragCoord.y-.5)*512.)/x;"
   "vec2 e=v*r(g),y=floor((.5+.5*e)*65536.),d=mod(y,256.)/255.,o=floor(y/256.)/255.;"
   "gl_FragColor=vec4(d.x,o.x,d.y,o.y);"
 "}";

#endif // SFX_H_
