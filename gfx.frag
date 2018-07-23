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

const float pi = acos(-1.);
const vec2 c = vec2(1., 0.);

float t;

float rand(vec2 a0)
{
    return -1.+2.*fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

float rand3(vec3 a0)
{
    return .33333*(rand(a0.xy)+rand(a0.yz)+rand(a0.zx));
}

float smoothstep_noise1d(float x)
{
    float d = floor(x);
    x = fract(x);
    float x0 = rand(d*c.xx), x1 = rand((d+1.)*c.xx);
    return mix(x0, x1, smoothstep(0., 1., x));
}

float smoothstep_noise2d(vec2 x)
{
    vec2 d = floor(x);
    x = fract(x);
    float x00 = rand(d),
        x01 = rand(d+c.yx),
        x10 = rand(d+c.xy), 
        x11 = rand(d+c.xx);
    return mix(mix(x00, x01, smoothstep(0.,1., x.y)), mix(x10, x11, smoothstep(0.,1., x.y)), smoothstep(0.,1., x.x));
}

//TODO: optimize compile time
float smoothstep_noise3d(vec3 x)
{
    vec3 d = floor(x);
    x = fract(x);
    float x000 = rand3(d),
        x100 = rand3(d+c.xyy),
        x010 = rand3(d+c.yxy),
        x001 = rand3(d+c.yyx),
        x110 = rand3(d+c.xxy),
        x011 = rand3(d+c.yxx),
        x101 = rand3(d+c.xyx),
        x111 = rand3(d+c.xxx);
    return mix(
                mix(
                    mix(x000, x001, smoothstep(0., 1., x.z)),
                    mix(x100, x101, smoothstep(0., 1., x.z)),
                    smoothstep(0., 1., x.x)),
       			mix(
                    mix(x010, x011, smoothstep(0., 1., x.z)),
                    mix(x110, x111, smoothstep(0., 1., x.z)),
                    smoothstep(0., 1., x.x)),
        		smoothstep(0., 1., x.y)
    		);
}

float mfsmoothstep_noise2d(vec2 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.2;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*smoothstep_noise2d(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
}

float mfsmoothstep_noise3d(vec3 x, float f0, float f1, float phi)
{
    float sum = 0.;
    float a = 1.;//1.2;
    
    for(float f = f0; f<f1; f = f*2.)
    {
        sum = a*smoothstep_noise3d(f*x) + sum;
        a = a*phi;
    }
    
    return sum;
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

void wave(in int i, out float st, out float am, out vec2 di, out float fr, out float sp)
{
    //setup wave params
	//st = abs(.035*rand(vec2(float(i))));//qi
	st = 55.*rand(float(-i)*c.xx);
    //am = .002*mod(float(i),4.5)+.005*mod(float(i),3.8)*rand(vec2(float(i+2)));//ai
    am = .002*sin(-3.*float(i)+12.*rand(float(i)*c.xx));
    //di = (4.*mod(float(5*i),2.2)+mod(float(i),2.6)*5.*vec2(1.7e0*rand(vec2(i,i+1)), 2.e0*rand(vec2(i+1,i))));//di
    di = 2.*vec2(rand(float(i)*c.xx), rand(float(i+1)*c.xx));
    fr = pow(.5, float(i))*(1.+.1*rand(float(i)*c.xx))*5.e1;
    sp = 12.e-1*mod(float(i*3),2.5)*rand(vec2(float(i+4)));//phi
}

float gerst(vec2 x, int nwaves)
{
    //vec3 val = vec3(x.xy, 0.);
    float val = 0.;
    
   	float st,fr,sp,am;
    vec2 di;
    
    for(int i=0; i<nwaves; ++i)
    {
   		wave(i, st, am, di, fr, sp);
        
        //gen values
        float d = dot(di, x);
        val += am*sin(-12.*fr*d+sp*t);
		//val += vec3(st*am*di*cos(fr*d+sp*t), am*sin(fr*d+sp*t));
    }

    return val;
}

//returns vec2(sdf, material)
vec2 z10presents(vec2 x)
{
    vec2 sda = vec2(b(x, vec2(-.5, -.4), vec2(iResolution.x/iResolution.y,-.4), .05), 1.),
        sdb = vec2(b(x, vec2(-.5, -.4), vec2(iResolution.x/iResolution.y,-.4), .02), 2.);
    return mix(sda, sdb, step(1., sdb.x));
}

//dusk scene
vec2 scene1(vec3 x)
{
    x += c.xyy*5.e-2*t;
    return vec2(x.z+.5-.5*mfsmoothstep_noise2d(x.xy, .9, 1., .1)+.05*mfsmoothstep_noise2d(x.xy, 120., 180., 1.7), 1.);
}

//lake scene
vec2 scene2(vec3 x)
{
    x += c.xxy * 14.+c.yxy*5.e-2*t;
    float k =- .5*mfsmoothstep_noise2d(x.xy, .9, 1., .1);
    if(x.z>-.46) k += .05*mfsmoothstep_noise2d(x.xy, 120., 180., 1.7);
    vec2 sda = vec2(x.z+k, 1.),
        sdb = vec2(x.z+.5-gerst(x.xy, 22), 2.);
    return mix(sda,sdb,step(sdb.x,sda.x));
}

vec2 scene(vec3 x)
{
    if(t < 2.) return scene1(x);//TODO: direction
    else if(t < 4000.) return scene2(x);
}

const float dx = 1.e-4;
vec3 normal(vec3 x)
{
    float s = scene(x).x;
    return normalize(vec3(
        scene(x + dx*c.xyy).x-s,
        scene(x + dx*c.yxy).x-s,
        scene(x + dx*c.yyx).x-s
    ));
}

vec3 bg(vec2 uv)
{
    uv += .01*iTime*c.yx;

    vec2 x = mod(uv, .05)-.025, y=uv-x;
    float sun = .17*exp(-2.e-1*t);
    return .15*c.xxx+mix(.4*c.xyy, .8*c.xxy, 1.-2.*(uv.y)-length(uv))*exp(-2.e-1*t)//sky
        +c.xxx*smoothstep(sun+.01, sun-.01,length(uv))//sun
        +2.*tanh(2.e-1*t)*smoothstep(.003, -.003, length(x-.02*vec2(rand(y), rand(.1*y+.1*c.xx)))-.002*rand(.1*y))
        *(.2*vec3(rand(.1*y),rand(.1*y+.1*c.xx),rand(y+.3*c.xx))+.8)//stars
    	+(.2+.1*mfsmoothstep_noise2d(uv+cos(.5*t)*sin(.5*uv.x)*sin(.6*uv.y), 2.5, 1.e2, .55))*exp(-5.e-2*t);//clouds
}

float B(float ton)
{
    return smoothstep(ton-1., ton, iTime)*(1.-smoothstep(ton+5., ton+6., iTime));
}

void fore(out vec4 fragColor, in vec2 uv, float time)
{
    t = time;
    
    vec2 s;
    vec3 x, o = .5+.5*mfsmoothstep_noise2d(c.yy, .9, 1., .1)-c.yxy+.5*c.yyx, t = .25*c.yyx, r = c.xyy, u = c.yyx, 
        rt = t + r * uv.x + u * uv.y, rd = normalize(rt-o);
    
    //raymarching
    float d = 0., vc = 0., cd = 0.;
    for(int i=0; i<200; ++i)
    {
        x = o + d * rd;
        s = scene(x);
        if(s.x<2.e-4)break;
        if((d>15.) || (i==199))
        {
            fragColor = vec4(mix(bg(uv),c.xxx,vc/cd), 1.);
            return;
        }
        d += s.x;
        
        //volumetric clouds
        cd += 5.e-1;
	    vc += smoothstep(0.,1.5, x.z)*mfsmoothstep_noise3d(o+cd*rd, 1., 100., .35) ;
    }
    
	//illumination
    vec3 n = normal(x), l = c.yxy, l2 = x+c.yyx, v = normalize(x-o), re = normalize(reflect(-l, n)), 
        re2 = normalize(reflect(-l2, n)), col;
    if(s.y == 1.)
    {
        col = .01*c.xyy*dot(l,n)+.2*c.xxx*pow(abs(dot(re,v)), 4.)+(.05*c.xyy+.05*c.yxy)*pow(abs(dot(re,v)), 2.);
    	col *= exp(-2.e-1*iTime);
    }
    else if(s.y == 2.)
    {
        col = .4*c.xxx*dot(l2,n)+.2*c.xxx*pow(abs(dot(re2,v)), 2.);
    }
    
    //fog
    col = mix(col, .15*c.xxx, cosh(-2.e-0*x.z)*tanh(1.09e-1*x.y));
    
    fragColor = vec4(col,1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //camera setup
    vec2 uv = fragCoord/iResolution.yy-.5, s;
    
    vec4 dt1, dt2, dt3;
    fore(dt1, uv, iTime);
    fore(dt2, uv, iTime+1.e-2);
    //fore(dt3, uv, iTime+2.e-2);
    fragColor = .333*(dt1+dt2);
    
    vec3 col = fragColor.xyz;
    
    //banner for text
    s = z10presents(uv);
    float sc = step(s.x, 0.) * B(5.); //objects
    if(s.y == 1.)
        col = mix(col, mix(col, c.xxx, .1), sc);
    
    //lens effect
    vec2 u = uv+(.01*iTime-.12)*c.yx;
    float ph = atan(abs(u.y/u.x)), ra = length(u);
    col += 1.4*(.5+.25*sin(2.3*pi*ph)+.25*sin(4.*pi*ph)+.25*sin(1.5*pi*ph))*exp(-3.*ra)*c.xxx*exp(-9.e-1*t)*(1.-smoothstep(1.4,1.6,iTime));
    
    //text
    float fx = b(uv, vec2(-.45, -.43), vec2(-.45, -.37), .005);//T
    fx = min(fx, b(uv, vec2(-.43,-.37), vec2(-.47, -.37), .005));//T
    fx = min(fx, b(uv, vec2(-.41,-.43), vec2(-.41,-.37), .005));//H
    fx = min(fx, b(uv, vec2(-.41,-.40), vec2(-.37,-.40), .005));//H
    fx = min(fx, b(uv, vec2(-.37, -.43), vec2(-.37, -.37), .005));//H
    fx = min(fx, b(uv, vec2(-.35, -.43), vec2(-.35, -.37), .005));//E
    fx = min(fx, b(uv, vec2(-.35, -.43), vec2(-.31, -.43), .005));//E
    fx = min(fx, b(uv, vec2(-.35, -.40), vec2(-.33, -.40), .005));//E
    fx = min(fx, b(uv, vec2(-.35, -.37), vec2(-.31, -.37), .005));//E
    fx = min(fx, b(uv, vec2(-.25, -.37), vec2(-.23, -.37), .005));//S
    fx = min(fx, cs(uv-vec2(-.25, -.385), .015, .005, .5*pi, 1.5*pi));//S
    fx = min(fx, cs(uv-vec2(-.25, -.385), .015, .005, -2.5*pi, -.5*pi));//S
    fx = min(fx, b(uv, vec2(-.25, -.40), vec2(-.24, -.40), .005));//S
    col = mix(col, .8*c.xxx, B(5.)*smoothstep(.002, -.002, fx ));
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
