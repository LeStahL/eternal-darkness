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
 
#version 130

uniform float iProgress;
uniform float iTime;
uniform vec2 iResolution;

const vec2 c = vec2(1.,0.);
const float pi = acos(-1.);

#define rand(a0) fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453)

float cr(vec2 x, float r, float w)
{
    return abs(length(x)-r)-w;
}

float cs(vec2 x, float r0, float w, float p0, float p1)
{
    float r = length(x), p = acos(x.x/r)*step(0.,x.y)-acos(x.x/r)*step(x.y,0.);
    p = clamp(p, p0, p1);
    vec2 y = r0*vec2(cos(p), sin(p));
    return length(x-y)-w;
}

float b(vec2 x, vec2 a, vec2 b, float w)
{
    vec2 d = b-a;
    return length(x-mix(a, b, clamp(dot(x-a, d)/dot(d,d), 0., 1.)))-w;
}

mat3 R(vec3 t)
{
    vec3 ct = cos(t), st = sin(t);
    return mat3(c.xyyy, ct.x, st.x, 0., -st.x, ct.x)
        *mat3(ct.y, 0., -st.y, c.yxy, st.y, 0., ct.y)
        *mat3(ct.z, st.z, 0., -st.z, ct.z, c.yyyx);
}

vec2 logo(vec2 x)
{
    vec2 sdf = c.xy;
    
    //210
    float sd = min(cr(x-.125*c.xy, .125, .04), cs(x+.125*c.xy, .125, .04, -pi/2., pi/2.));
    sd = min(sd, b(x, -.125*c.yx, .125*c.yx, .04));
    sdf = vec2(sd, 1.);
	vec2 sdb = vec2(abs(sd-.01)-.005, 2.);
    sdf = mix(sdb, sdf, step(sdf.x, sdb.x));
    
    //Team
    //T
    sd = min(
        b(x, vec2(-.1, -.25), vec2(-.1, -.2), .01),
        b(x, vec2(-.12,-.2), vec2(-.08, -.2), .01)
    );
    //E
    sd = min(sd, b(x, vec2(-.05, -.25), vec2(-.05, -.2), .01));
    sd = min(sd, b(x, vec2(-.05,-.2), vec2(-.02, -.2), .01));
    sd = min(sd, b(x, vec2(-.05,-.25), vec2(-.02, -.25), .01));
    sd = min(sd, b(x, vec2(-.05, -.225), vec2(-.03, -.225), .01));
    //A
    sd = min(sd, cs(x-vec2(.03, -.22), .02, .01, 0., pi));
    sd = min(sd, b(x, vec2(.01, -.22), vec2(.01, -.25), .01));
    sd = min(sd, b(x, vec2(.05, -.22), vec2(.05, -.25), .01));
    sd = min(sd, b(x, vec2(.01, -.23), vec2(.05, -.23), .01));
    //M
    sd = min(sd, b(x, vec2(.08, -.2), vec2(.08, -.25), .01));
    sd = min(sd, b(x, vec2(.08, -.2), vec2(.1, -.22), .01));
    sd = min(sd, b(x, vec2(.1, -.22), vec2(.12, -.2), .01));
    sd = min(sd, b(x, vec2(.12, -.2), vec2(.12, -.25), .01));
    sdb = vec2(abs(sd-.01)-.005, 2.);
    sdf = mix(sdb, sdf, step(sdf.x, sdb.x));
    
    sdb = vec2(sd, 1.);
    sdf = mix(sdb, sdf, step(sdf.x, sdb.x));
    
    //progress bar
    sd = b(x, vec2(-.4, -.35), vec2(mix(-.4,.4, iProgress), -.35), .02);
    sdb = vec2(sd, 3.);
    sdf = mix(sdb, sdf, step(sdf.x, sdb.x));
    sd = b(x, vec2(-.4, -.35), vec2(.4, -.35), .02);
	sdb = vec2(abs(sd-.01)-.005, 2.);
    sdf = mix(sdb, sdf, step(sdf.x, sdb.x));
    
    return sdf;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = fragCoord.xy/iResolution.yy-.5*c.xx-.33*c.xy;

    vec2 k = logo(uv);
    vec3 col = mix(c.xxy,c.xyx,k.x);
    fragColor = c.yyyx + .6*c.xxxy*smoothstep(.005, -.005, k.x)+.3*rand(uv-1.e-3*iTime*c.xx)*c.xxxx;
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
