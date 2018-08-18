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

uniform float iTime;
uniform vec2 iResolution;

const float pi = acos(-1.);
const vec3 c = vec3(1., 0., -1.);

vec2 vi;
float t;
#define T .5

float rand(vec2 a0)
{
    return -1.+2.*fract(sin(dot(a0.xy ,vec2(12.9898,78.233)))*43758.5453);
}

float rand3(vec3 a0)
{
    return .33333*(rand(a0.xy)+rand(a0.yz)+rand(a0.zx));
}

//distance to quadratic bezier spline with parameter t
float dist(vec2 p0,vec2 p1,vec2 p2,vec2 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum distance to quadratic bezier spline
float dsp(vec2 p0, vec2 p1, vec2 p2, vec2 x)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec2 E = x-p0, F = p2-2.*p1+p0, G = p1-p0;
    vec3 ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

	//discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        return dist(p0,p1,p2,x,ui.x+ui.y-tau);
    }
    
    //three distinct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    return min(
        dist(p0,p1,p2,x, t.x),
        min(
            dist(p0,p1,p2,x,t.y),
            dist(p0,p1,p2,x,t.z)
        )
    );
}

//distance to 3d spline with parameter t
float dist(vec3 p0,vec3 p1,vec3 p2,vec3 x,float t)
{
    t = clamp(t, 0., 1.);
    return length(x-pow(1.-t,2.)*p0-2.*(1.-t)*t*p1-t*t*p2);
}

//minimum distance to worm
float worm(vec3 p0, vec3 p1, vec3 p2, vec3 x)
{
    //coefficients for 0 = t^3 + a * t^2 + b * t + c
    vec3 E = x-p0, F = p2-2.*p1+p0, G = p1-p0,
    	ai = vec3(3.*dot(G,F), 2.*dot(G,G)-dot(E,F), -dot(E,G))/dot(F,F);

	//discriminant and helpers
    float tau = ai.x/3., p = ai.y-tau*ai.x, q = - tau*(tau*tau+p)+ai.z, dis = q*q/4.+p*p*p/27.;
    
    //triple real root
    if(dis > 0.) 
    {
        vec2 ki = -.5*q*c.xx+sqrt(dis)*c.xz, ui = sign(ki)*pow(abs(ki), c.xx/3.);
        return dist(p0,p1,p2,x,ui.x+ui.y-tau)-(.03+.01*mod(t,T)/T-.005*cos(1.7*pi*(ui.x+ui.y-tau))+.01*abs(sin(2.*pi*10.*(ui.x+ui.y-tau))));
    }
    
    //three distinct real roots
    float fac = sqrt(-4./3.*p), arg = acos(-.5*q*sqrt(-27./p/p/p))/3.;
    vec3 t = c.zxz*fac*cos(arg*c.xxx+c*pi/3.)-tau;
    return min(
        dist(p0,p1,p2,x, t.x)-(.03+.01*mod(iTime,T)/T-.005*cos(1.7*pi*t.x)+.01*abs(sin(2.*pi*10.*t.x))),
        min(
            dist(p0,p1,p2,x,t.y)-(.03+.01*mod(iTime,T)/T-.005*cos(1.7*pi*t.y)+.01*abs(sin(2.*pi*10.*t.y))),
            dist(p0,p1,p2,x,t.z)-(.03+.01*mod(iTime,T)/T-.005*cos(1.7*pi*t.z)+.01*abs(sin(2.*pi*10.*t.z)))
        )
    );
}

//minimum distance to linear bezier spline
float dsg(vec2 p0, vec2 p1, vec2 x)
{
    vec2 d = p1-p0;
    float t = clamp(dot(x-p0,d)/dot(d,d),0.,1.);
    return length(x-mix(p0,p1,t));
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
/*
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
}*/

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

/*
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
}*/

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
    am = 1.2e-4*pow(.9, float(i-1))*rand(float(i)*c.xx);
    //di = (4.*mod(float(5*i),2.2)+mod(float(i),2.6)*5.*vec2(1.7e0*rand(vec2(i,i+1)), 2.e0*rand(vec2(i+1,i))));//di
    //di = 2.*vec2(rand(float(i)*c.xx), rand(float(i+1)*c.xx));
    di = c.yx*rand(float(i)*c.xx)+.3*c.xy*rand(float(i)*c.xy);
    fr = pow(.5+.1*rand(float(i+1)*c.xx), float(i))*(1.+.2*rand(float(i)*c.xx))*7.e2;
    sp = 12.e0*mod(float(i*3),2.5)*rand(vec2(float(i+4)));//phi
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
        val += am*sin(fr*d+sp*t);
		//val += vec3(st*am*di*cos(fr*d+sp*t), am*sin(fr*d+sp*t));
    }

    return val;
}

/* compute voronoi distance and closest point.
 * x: coordinate
 * return value: vec3(distance, coordinate of control point)
 */
vec3 vor(vec2 x)
{
    x += 12.;// mod(x,6.*2.*pi);
    vec2 y = floor(x);
   	float ret = 1.;
    
    //find closest control point. ("In which cell am I?")
    vec2 pf=c.yy, p;
    float df=100., d;
    
    vec3 a = c.xyy*100.;
    
    for(float i=-2.; i<=2.; i+=1.)
        for(float j=-2.; j<=2.; j+=1.)
        {
            p = y + vec2(i, j);
            p += rand(p);
            
            d = length(x-p);
            
            a = mix(a, vec3(d,p), step(d, a.x));
        }
    
    //compute voronoi distance: minimum distance to any edge
    for(float i=-2.; i<=2.; i+=1.)
        for(float j=-2.; j<=2.; j+=1.)
        {
            p = y + vec2(i, j);
            p += rand(p);
            
            vec2 o = p - a.yz;
            d = dot(.5*o-x+a.yz, o)/length(o);
            ret = min(ret, d);
        }
    
    return vec3(ret, pf);
}

vec3 rot(vec3 x, vec3 theta)
{
    vec3 c = cos(theta), s = sin(theta);
    return mat3(c.x,0.,s.x,0.,1.,0.,-s.x,0.,c.x)
        *mat3(1., 0., 0., 0., c.y, -s.y, 0., s.y, c.y)
        *mat3(c.z,-s.z,0., s.z,c.z,0., 0.,0.,1.)*x;
}

// compute distance to regular star
float dstar(vec2 x, float N, vec2 R)
{
    float d = pi/N,
        p0 = acos(x.x/length(x)),
        p = mod(p0, d),
        i = mod(round((p-p0)/d),2.);
    x = length(x)*vec2(cos(p),sin(p));
    vec2 a = mix(R,R.yx,i),
    	p1 = a.x*c.xy,
        ff = a.y*vec2(cos(d),sin(d))-p1;
   	ff = ff.yx*c.zx;
    return dot(x-p1,ff)/length(ff);
}

// compute distance to regular polygon
float dpoly_min(vec2 x, float N, float R)
{
    float d = 2.*pi/N,
        t = mod(acos(x.x/length(x)), d)-.5*d;
    return R-length(x)*cos(t)/cos(.5*d);
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
    x += -13.*c.yxy + c.xxy * 25.+c.yxy*5.e-2*t;
    float k =- .5*mfsmoothstep_noise2d(x.xy, .91, 2., .1);
    if(x.z>-.19) k += .05*mfsmoothstep_noise2d(x.xy, 120., 180., 1.7);
    vec2 sda = vec2(x.z+k, 1.),
        sdb = vec2(x.z+.21-gerst(x.xy, 22), 2.);
    return mix(sda,sdb,step(sdb.x,sda.x));
}

//tree fog scene
vec2 scene3(vec3 x)
{
    x += c.yxy*5.e-1*t;
    vec3 y = vec3(mod(x.x, 2.)-1., mod(x.y, 2.)-1., x.z), z = x-y;
    vec2 sda = vec2(length(y.xy-.5*vec2(rand(z.xy),rand(z.xy*1.1)))-.1, 1.), 
        sdb = vec2(x.z-.2*mfsmoothstep_noise2d(x.xy, .9, 1., .1), 1.);
    float guard = -length(max(abs(y)-vec3(1.,2.,50.), 0.));
    guard = abs(guard)+2.*.2;
    sda.x = min(sda.x, guard);
    
    return mix(sdb, sda, step(sda.x, sdb.x));
}

//worm scene
vec2 scene4(vec3 x)
{
    vec3 x0 = x;
    x += 1.e0*t*c.yxy;
    
    vec3 y = x;
    x = mod(x,1.)-.5;
    
    vec3 i = y-x;
    float l = length(i);
    vec3
        dx = .25*vec3(rand(l*c.xx), rand(3.*l*c.xx), rand(6.*l*c.xx));
    
    vec2 sda = vec2(worm(c.yyy, vec3(.1*sin(5.*t+10.*rand(length(y-x)*c.xx)),0.,.1), .2*c.yyx, x-dx), 3.),
        sdb = vec2(length(x-dx)-(.036+.01*mod(t,T)/T),4.),
    	sdf = mix(sda,sdb, step(sdb.x,sda.x)),
        sdc = vec2(length(x-.2*c.yyx-dx)-(.036+.01*mod(t,T)/T),4.);
    sdf = mix(sdf, sdc, step(sdc.x, sdf.x));

   //vec2 sdf=c.xy,  sda;

    float phi =3.*acos(min(y.z,y.x)/length(y.xz));
    vec3 vp = .5*vor(vec2(y.y+phi,phi));//+.1*vor(2.*vec2(phi, y.y));//-.02*vor(5.*vec2(phi,y.y));
    //vec3 vp = .5*vor(y.xy+20.);
    float v = vp.x;
    vi = vp.yz;
    sda = vec2(abs(length(y.xz)-3.)-v, 5.);
    //sda = vec2(y.z+2.5-v, 5.);
    sdf = mix(sdf, sda, step(sda.x, sdf.x));
    

    float guard = -length(max(abs(x)-.5,0.));
    guard = abs(guard)+1.*.1;
    sdf.x = min(sdf.x, guard);

    return sdf;
}

//party scene
vec2 scene5(vec3 x)
{   
    vec3 x0 = x;
    x += 1.e0*t*c.yxy;
    
    vec3 y = x;
    x = vec3(mod(x.xy,1.)-.5, x.z);

    vec3 index = y-x;
    vi = index.xy;
    float ang = 2.*pi*rand(index.xy)*t;
    x.xy = mat2(cos(ang), sin(ang), -sin(ang), cos(ang))*x.xy;

    float ds = dstar(x.xy, 3.+floor(7.*abs(rand(index.xy+6.))), vec2(.2+.05*rand(vi),.3+.05*rand(vi+3.)));
    
    //vec2 sdf = vec2(length(vec2(min(ds,0.), min(abs(x.z+1.3), 0.05))) ,6.);
    vec2 d = abs(vec2(min(ds,0.),x.z+1.3)) -.4*c.yx- (.5+2.*mod(t,T))*.5*abs(rand(index.xy+t-mod(t,T)))*c.yx;
    vec2 sdf = vec2(min(max(d.x,d.y),0.) + length(max(d,0.0)),6.);
    

    //ds = dstar(x.xy, 3.+7.*rand(index), vec2(.2,.4));
    
    //vec2 sda = vec2(length(vec2(min(ds,0.), min((x.z-1.3+(.5+mod(t,T))*.5*rand(index.xy+t-mod(t,T))), 0.))) ,6.);
    d = abs(vec2(min(ds,0.),x.z-1.3)) -.4*c.yx- (.5+2.*mod(t,T))*.5*rand(index.xy+t-mod(t,T))*c.yx;
    vec2 sda = vec2(min(max(d.x,d.y),0.) + length(max(d,0.0)),6.);
    
    sdf = mix(sdf, sda, step(sda.x, sdf.x));
    
    
    float guard = -length(max(abs(x+50.*c.yyx)-vec3(.5,.5,100.),0.));
    guard = abs(guard)+1.*.1;
    sdf.x = min(sdf.x, guard);
    
    return sdf;
}

vec2 scene(vec3 x)
{
    if(t < 20.) return scene1(x);//TODO: direction
    else if(t < 30.) return scene2(x);
    else if(t < 40.) return scene3(x);
    else if(t < 50.) return scene4(x);
    else if(t < 7000.) return scene5(x);
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
    if(t<20.)
    {
    	uv += .01*t*c.yx;

    	vec2 x = mod(uv, .05)-.025, y=uv-x;
    	float sun = .17*exp(-2.e-1*t);
    	return .15*c.xxx+mix(.4*c.xyy, .8*c.xxy, 1.-2.*(uv.y)-length(uv))*exp(-2.e-1*t)//sky
        	+c.xxx*smoothstep(sun+.01, sun-.01,length(uv))//sun
        	+2.*tanh(2.e-1*t)*smoothstep(.003, -.003, length(x-.02*vec2(rand(y), rand(.1*y+.1*c.xx)))-.002*rand(.1*y))
        	*(.2*vec3(rand(.1*y),rand(.1*y+.1*c.xx),rand(y+.3*c.xx))+.8)//stars
    		+(.4+.2*mfsmoothstep_noise2d(uv+cos(.5*t)*sin(.5*uv.x)*sin(.6*uv.y), 2.5, 1.e2, .55))*exp(-5.e-2*t);//clouds
    }
    else 
    {
        return c.xxx;
    }
}

float B(float ton)
{
    return smoothstep(ton-2., ton-1., iTime)*(1.-smoothstep(ton+2., ton+3., iTime));
}

void fore(out vec4 fragColor, in vec2 uv, float time)
{
    t = time;
    
    vec2 s;
    vec3 x, o = .5+.5*mfsmoothstep_noise2d(c.yy, .9, 1., .1)-c.yxy+.5*c.yyx, ta = .25*c.yyx, r = c.xyy, u = c.yyx, 
        rt = ta + r * uv.x + u * uv.y, rd = normalize(rt-o);
    if(t>40.)
    {
        o = -c.yxy;
        ta=c.yyy;
    }
    //raymarching
    float d = 0., vc = 0., cd = 0., ci=15.;
    int ni = 100;
    if(time > 20.) ci = 55.;
    if(time > 30.) 
    {
        ci = 50.;
        ni = 200;
    }
    else if(time > 40.) ci = 200.;
    for(int i=0; i<ni; ++i)
    {
        x = o + d * rd;
        s = scene(x);
        if(s.x<5.e-3)break;
        if((d>ci) || (i==ni-1))
        {
            //fragColor = vec4(mix(bg(uv),c.xxx,vc/cd), 1.);
            fragColor = vec4(bg(uv), 1.);
            return;
        }
        d += s.x;
        
        //volumetric clouds
        //cd += 5.e-1;
	    //vc += smoothstep(0.,1.5, x.z)*mfsmoothstep_noise3d(o+cd*rd, 1., 100., .35) ;
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
    else if(s.y == 3.)//worm body
    {
        col = .6-(.8*c.yyx+.3*c.xyy)*dot(l,n)*mod(t,T)/T+.8*c.xxx*dot(l,n)+.2*c.xxx*pow(abs(dot(re,v)), 4.)+(.05*c.xyy+.05*c.yxy)*pow(abs(dot(re,v)), 2.);
    }
    else if(s.y == 4.)//worm ends
    {
        col = 1.+.8*c.xxx*dot(l,n)+.2*c.xxx*pow(abs(dot(re,v)), 4.)+(.05*c.yyx+.05*c.yyx)*pow(abs(dot(re,v)), 2.);
    }
    else if(s.y == 5.)//worm tunnel
    {
        l = normalize(x);
        re = normalize(reflect(-l,n));
        v = normalize(x-o);
        col = .6+.3*vec3(.75+.1*rand(vi),.75+.1*rand(2.*vi), .6+.4*rand(5.*vi))*dot(l,n)*mod(t,T)/T+.8*c.xxx*dot(l,n)+.05*c.xxx*pow(abs(dot(re,v)), 2.);
    }
    else if(s.y == 6.)//dance
    {
        l = normalize(x-c.yxy);
        re = normalize(reflect(-l,n));
        v = normalize(x-c.yxy-o);
        vec3 co = .5 + .5*cos(x+vec3(0.,2.,4.)+5.e-1*abs(rand(vi))+time), 
            co2 = .5 + .5*cos(x+vec3(0.,2.,4.)+6.e-1*abs(rand(vi+3.))+time),
            co3 = .5 + .5*cos(x+vec3(0.,2.,4.)+7.e-1*abs(rand(vi+7.))+time);
        col = co+co2*dot(l,n)+.1*co3*pow(abs(dot(re,v)), 2.);
    }
    
    else col = c.yyy;
    
    //fog
    if(time < 20.)
	    col = mix(col, .15*c.xxx, cosh(-2.e-0*x.z)*tanh(1.09e-1*x.y));
    else if(time < 40.)
        col = mix(col, c.xxx, tanh(1.e-1*x.y));
    else if(time < 7000.)
        col = mix(col, c.xxx, tanh(8.e-2*x.y));
    fragColor = vec4(col,1.0);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //camera setup
    vec2 uv = fragCoord/iResolution.yy-.5, s;
    
    vec4 dt1, dt2, dt3;
    fore(dt1, uv, iTime);
    if(iTime < 40.)
    {
        fore(dt2, uv, iTime+1.e-2);
        //fore(dt3, uv, iTime+2.e-2);
        fragColor = .333*(dt1+dt2);
    }
    else
    {
        fragColor = dt1;
    }

    vec3 col = fragColor.xyz;
    
    //banner for text
    if(iTime < 180.)
    {
        s = z10presents(uv);
        float sc = step(s.x, 0.) ; //objects
        if(s.y == 1.)
            col = mix(col, mix(col, c.xxx, .1), sc);
     }

    //lens effect
    vec2 u = uv+(.01*iTime-.12)*c.yx;
    float ph = atan(abs(u.y/u.x)), ra = length(u);
    col += 1.4*(.5+.25*sin(2.3*pi*ph)+.25*sin(4.*pi*ph)+.25*sin(1.5*pi*ph))*exp(-3.*ra)*c.xxx*exp(-9.e-1*t)*(1.-smoothstep(1.4,1.6,iTime));
    
    //text
    float d = 1.;
    if(iTime < 8.)
    {
        const vec2 lin[84] = vec2[84](vec2(-3.88e-01,-4.15e-01),vec2(-3.88e-01,-3.68e-01),vec2(-4.00e-01,-3.68e-01),vec2(-3.77e-01,-3.68e-01),vec2(-3.53e-01,-4.15e-01),vec2(-3.41e-01,-4.15e-01),vec2(-3.64e-01,-4.03e-01),vec2(-3.41e-01,-4.03e-01),vec2(-3.05e-01,-4.15e-01),vec2(-3.05e-01,-3.92e-01),vec2(-2.93e-01,-4.15e-01),vec2(-2.93e-01,-3.92e-01),vec2(-2.70e-01,-4.15e-01),vec2(-2.70e-01,-4.03e-01),vec2(-2.47e-01,-4.15e-01),vec2(-2.47e-01,-4.03e-01),vec2(-2.34e-01,-4.15e-01),vec2(-2.11e-01,-4.15e-01),vec2(-1.99e-01,-3.80e-01),vec2(-1.87e-01,-3.68e-01),vec2(-1.87e-01,-3.68e-01),vec2(-1.87e-01,-4.15e-01),vec2(-1.75e-01,-4.03e-01),vec2(-1.75e-01,-3.80e-01),vec2(-1.51e-01,-4.03e-01),vec2(-1.51e-01,-3.80e-01),vec2(-1.15e-01,-4.38e-01),vec2(-1.15e-01,-3.92e-01),vec2(-1.15e-01,-3.92e-01),vec2(-1.04e-01,-3.92e-01),vec2(-1.15e-01,-4.15e-01),vec2(-1.04e-01,-4.15e-01),vec2(-7.97e-02,-4.15e-01),vec2(-7.97e-02,-3.92e-01),vec2(-2.02e-02,-3.92e-01),vec2(-2.02e-02,-4.03e-01),vec2(3.08e-03,-3.92e-01),vec2(3.08e-03,-4.15e-01),vec2(2.70e-02,-3.92e-01),vec2(3.87e-02,-3.92e-01),vec2(2.70e-02,-4.15e-01),vec2(3.87e-02,-4.15e-01),vec2(3.87e-02,-4.15e-01),vec2(3.87e-02,-3.68e-01),vec2(5.09e-02,-4.03e-01),vec2(5.09e-02,-3.68e-01),vec2(7.48e-02,-3.92e-01),vec2(7.48e-02,-4.03e-01),vec2(9.82e-02,-3.92e-01),vec2(9.82e-02,-4.27e-01),vec2(8.65e-02,-4.38e-01),vec2(7.48e-02,-4.38e-01),vec2(1.34e-01,-4.38e-01),vec2(1.34e-01,-3.92e-01),vec2(1.34e-01,-3.92e-01),vec2(1.46e-01,-3.92e-01),vec2(1.34e-01,-4.15e-01),vec2(1.46e-01,-4.15e-01),vec2(1.70e-01,-4.15e-01),vec2(1.70e-01,-3.92e-01),vec2(2.06e-01,-4.15e-01),vec2(2.17e-01,-4.15e-01),vec2(1.94e-01,-4.03e-01),vec2(2.17e-01,-4.03e-01),vec2(2.29e-01,-4.15e-01),vec2(2.41e-01,-4.15e-01),vec2(2.41e-01,-3.92e-01),vec2(2.53e-01,-3.92e-01),vec2(2.77e-01,-4.15e-01),vec2(2.88e-01,-4.15e-01),vec2(2.65e-01,-4.03e-01),vec2(2.88e-01,-4.03e-01),vec2(3.01e-01,-4.15e-01),vec2(3.01e-01,-3.92e-01),vec2(3.24e-01,-4.03e-01),vec2(3.24e-01,-4.15e-01),vec2(3.48e-01,-4.03e-01),vec2(3.48e-01,-3.68e-01),vec2(3.36e-01,-3.92e-01),vec2(3.60e-01,-3.92e-01),vec2(3.72e-01,-4.15e-01),vec2(3.83e-01,-4.15e-01),vec2(3.83e-01,-3.92e-01),vec2(3.95e-01,-3.92e-01)),
        quad[141] = vec2[141](vec2(-3.41e-01,-4.03e-01),vec2(-3.41e-01,-3.92e-01),vec2(-3.53e-01,-3.92e-01),vec2(-3.53e-01,-3.92e-01),vec2(-3.64e-01,-3.92e-01),vec2(-3.64e-01,-4.03e-01),vec2(-3.64e-01,-4.03e-01),vec2(-3.64e-01,-4.15e-01),vec2(-3.53e-01,-4.15e-01),vec2(-3.17e-01,-3.92e-01),vec2(-3.29e-01,-3.92e-01),vec2(-3.29e-01,-4.03e-01),vec2(-3.29e-01,-4.03e-01),vec2(-3.29e-01,-4.15e-01),vec2(-3.17e-01,-4.15e-01),vec2(-3.17e-01,-4.15e-01),vec2(-3.05e-01,-4.15e-01),vec2(-3.05e-01,-4.03e-01),vec2(-3.17e-01,-3.92e-01),vec2(-3.05e-01,-3.92e-01),vec2(-3.05e-01,-4.03e-01),vec2(-2.93e-01,-4.03e-01),vec2(-2.93e-01,-3.92e-01),vec2(-2.82e-01,-3.92e-01),vec2(-2.82e-01,-3.92e-01),vec2(-2.70e-01,-3.92e-01),vec2(-2.70e-01,-4.03e-01),vec2(-2.70e-01,-4.03e-01),vec2(-2.70e-01,-3.92e-01),vec2(-2.58e-01,-3.92e-01),vec2(-2.58e-01,-3.92e-01),vec2(-2.47e-01,-3.92e-01),vec2(-2.47e-01,-4.03e-01),vec2(-2.34e-01,-3.68e-01),vec2(-1.99e-01,-3.57e-01),vec2(-2.34e-01,-4.15e-01),vec2(-1.75e-01,-4.03e-01),vec2(-1.75e-01,-4.15e-01),vec2(-1.63e-01,-4.15e-01),vec2(-1.63e-01,-4.15e-01),vec2(-1.51e-01,-4.15e-01),vec2(-1.51e-01,-4.03e-01),vec2(-1.51e-01,-3.80e-01),vec2(-1.51e-01,-3.68e-01),vec2(-1.63e-01,-3.68e-01),vec2(-1.63e-01,-3.68e-01),vec2(-1.75e-01,-3.68e-01),vec2(-1.75e-01,-3.80e-01),vec2(-1.04e-01,-4.15e-01),vec2(-9.20e-02,-4.15e-01),vec2(-9.20e-02,-4.03e-01),vec2(-9.20e-02,-4.03e-01),vec2(-9.20e-02,-3.92e-01),vec2(-1.04e-01,-3.92e-01),vec2(-7.97e-02,-4.03e-01),vec2(-7.97e-02,-3.92e-01),vec2(-6.81e-02,-3.92e-01),vec2(-4.42e-02,-4.15e-01),vec2(-5.58e-02,-4.15e-01),vec2(-5.58e-02,-4.03e-01),vec2(-5.58e-02,-4.03e-01),vec2(-5.58e-02,-3.92e-01),vec2(-4.42e-02,-3.92e-01),vec2(-4.42e-02,-3.92e-01),vec2(-3.25e-02,-3.92e-01),vec2(-3.25e-02,-4.03e-01),vec2(-3.25e-02,-4.03e-01),vec2(-3.25e-02,-4.15e-01),vec2(-4.42e-02,-4.15e-01),vec2(-2.02e-02,-4.03e-01),vec2(-2.02e-02,-4.15e-01),vec2(-8.58e-03,-4.15e-01),vec2(-8.58e-03,-4.15e-01),vec2(3.08e-03,-4.15e-01),vec2(3.08e-03,-4.03e-01),vec2(2.70e-02,-4.15e-01),vec2(1.53e-02,-4.15e-01),vec2(1.53e-02,-4.03e-01),vec2(1.53e-02,-4.03e-01),vec2(1.53e-02,-3.92e-01),vec2(2.70e-02,-3.92e-01),vec2(5.09e-02,-4.03e-01),vec2(5.09e-02,-4.15e-01),vec2(6.26e-02,-4.15e-01),vec2(7.48e-02,-4.03e-01),vec2(7.48e-02,-4.15e-01),vec2(8.65e-02,-4.15e-01),vec2(8.65e-02,-4.15e-01),vec2(9.82e-02,-4.15e-01),vec2(9.82e-02,-4.03e-01),vec2(9.82e-02,-4.27e-01),vec2(9.82e-02,-4.38e-01),vec2(8.65e-02,-4.38e-01),vec2(1.46e-01,-4.15e-01),vec2(1.58e-01,-4.15e-01),vec2(1.58e-01,-4.03e-01),vec2(1.58e-01,-4.03e-01),vec2(1.58e-01,-3.92e-01),vec2(1.46e-01,-3.92e-01),vec2(1.70e-01,-4.03e-01),vec2(1.70e-01,-3.92e-01),vec2(1.82e-01,-3.92e-01),vec2(2.17e-01,-4.03e-01),vec2(2.17e-01,-3.92e-01),vec2(2.06e-01,-3.92e-01),vec2(2.06e-01,-3.92e-01),vec2(1.94e-01,-3.92e-01),vec2(1.94e-01,-4.03e-01),vec2(1.94e-01,-4.03e-01),vec2(1.94e-01,-4.15e-01),vec2(2.06e-01,-4.15e-01),vec2(2.41e-01,-4.15e-01),vec2(2.64e-01,-4.09e-01),vec2(2.41e-01,-4.03e-01),vec2(2.41e-01,-4.03e-01),vec2(2.18e-01,-3.98e-01),vec2(2.41e-01,-3.92e-01),vec2(2.88e-01,-4.03e-01),vec2(2.88e-01,-3.92e-01),vec2(2.77e-01,-3.92e-01),vec2(2.77e-01,-3.92e-01),vec2(2.65e-01,-3.92e-01),vec2(2.65e-01,-4.03e-01),vec2(2.65e-01,-4.03e-01),vec2(2.65e-01,-4.15e-01),vec2(2.77e-01,-4.15e-01),vec2(3.12e-01,-3.92e-01),vec2(3.24e-01,-3.92e-01),vec2(3.24e-01,-4.03e-01),vec2(3.01e-01,-4.03e-01),vec2(3.01e-01,-3.92e-01),vec2(3.12e-01,-3.92e-01),vec2(3.48e-01,-4.03e-01),vec2(3.48e-01,-4.15e-01),vec2(3.60e-01,-4.15e-01),vec2(3.83e-01,-4.15e-01),vec2(4.07e-01,-4.09e-01),vec2(3.83e-01,-4.03e-01),vec2(3.83e-01,-4.03e-01),vec2(3.60e-01,-3.98e-01),vec2(3.83e-01,-3.92e-01));
        for(int i=0; i<42;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
        for(int i=0; i<47; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
        col = mix(col, .8*c.xxx, B(5.)*smoothstep(.005, .002, d ));
    }
	else if(iTime < 13.)
    {
        const vec2 lin[58] = vec2[58](vec2(-3.99e-01,-4.15e-01),vec2(-3.99e-01,-3.68e-01),vec2(-3.99e-01,-4.15e-01),vec2(-3.76e-01,-4.15e-01),vec2(-3.99e-01,-3.92e-01),vec2(-3.76e-01,-3.92e-01),vec2(-3.99e-01,-3.68e-01),vec2(-3.76e-01,-3.68e-01),vec2(-3.52e-01,-4.03e-01),vec2(-3.52e-01,-3.68e-01),vec2(-3.64e-01,-3.92e-01),vec2(-3.41e-01,-3.92e-01),vec2(-3.17e-01,-4.15e-01),vec2(-3.05e-01,-4.15e-01),vec2(-3.28e-01,-4.03e-01),vec2(-3.05e-01,-4.03e-01),vec2(-2.93e-01,-4.15e-01),vec2(-2.93e-01,-3.92e-01),vec2(-2.69e-01,-4.15e-01),vec2(-2.69e-01,-3.92e-01),vec2(-2.45e-01,-4.03e-01),vec2(-2.45e-01,-4.15e-01),vec2(-2.10e-01,-4.15e-01),vec2(-2.10e-01,-3.92e-01),vec2(-1.98e-01,-4.03e-01),vec2(-1.98e-01,-3.68e-01),vec2(-1.50e-01,-4.15e-01),vec2(-1.50e-01,-3.68e-01),vec2(-1.50e-01,-3.68e-01),vec2(-1.38e-01,-3.68e-01),vec2(-1.50e-01,-4.15e-01),vec2(-1.38e-01,-4.15e-01),vec2(-1.26e-01,-4.03e-01),vec2(-1.26e-01,-3.80e-01),vec2(-9.08e-02,-4.15e-01),vec2(-9.08e-02,-3.92e-01),vec2(-7.86e-02,-4.15e-01),vec2(-7.86e-02,-3.92e-01),vec2(-5.47e-02,-4.15e-01),vec2(-5.47e-02,-3.68e-01),vec2(-5.47e-02,-4.03e-01),vec2(-4.30e-02,-4.03e-01),vec2(-1.91e-02,-4.15e-01),vec2(-1.91e-02,-3.92e-01),vec2(4.25e-03,-4.03e-01),vec2(4.25e-03,-4.15e-01),vec2(2.82e-02,-4.15e-01),vec2(3.98e-02,-4.15e-01),vec2(1.65e-02,-4.03e-01),vec2(3.98e-02,-4.03e-01),vec2(5.21e-02,-4.15e-01),vec2(6.38e-02,-4.15e-01),vec2(6.38e-02,-3.92e-01),vec2(7.54e-02,-3.92e-01),vec2(8.77e-02,-4.15e-01),vec2(9.93e-02,-4.15e-01),vec2(9.93e-02,-3.92e-01),vec2(1.11e-01,-3.92e-01)),
        quad[90] = vec2[90](vec2(-3.52e-01,-4.03e-01),vec2(-3.52e-01,-4.15e-01),vec2(-3.41e-01,-4.15e-01),vec2(-3.05e-01,-4.03e-01),vec2(-3.05e-01,-3.92e-01),vec2(-3.17e-01,-3.92e-01),vec2(-3.17e-01,-3.92e-01),vec2(-3.28e-01,-3.92e-01),vec2(-3.28e-01,-4.03e-01),vec2(-3.28e-01,-4.03e-01),vec2(-3.28e-01,-4.15e-01),vec2(-3.17e-01,-4.15e-01),vec2(-2.93e-01,-4.03e-01),vec2(-2.93e-01,-3.92e-01),vec2(-2.81e-01,-3.92e-01),vec2(-2.57e-01,-3.92e-01),vec2(-2.45e-01,-3.92e-01),vec2(-2.45e-01,-4.03e-01),vec2(-2.69e-01,-4.03e-01),vec2(-2.69e-01,-3.92e-01),vec2(-2.57e-01,-3.92e-01),vec2(-2.22e-01,-3.92e-01),vec2(-2.33e-01,-3.92e-01),vec2(-2.33e-01,-4.03e-01),vec2(-2.33e-01,-4.03e-01),vec2(-2.33e-01,-4.15e-01),vec2(-2.22e-01,-4.15e-01),vec2(-2.22e-01,-4.15e-01),vec2(-2.10e-01,-4.15e-01),vec2(-2.10e-01,-4.03e-01),vec2(-2.22e-01,-3.92e-01),vec2(-2.10e-01,-3.92e-01),vec2(-2.10e-01,-4.03e-01),vec2(-1.98e-01,-4.03e-01),vec2(-1.98e-01,-4.15e-01),vec2(-1.86e-01,-4.15e-01),vec2(-1.38e-01,-3.68e-01),vec2(-1.26e-01,-3.68e-01),vec2(-1.26e-01,-3.80e-01),vec2(-1.38e-01,-4.15e-01),vec2(-1.26e-01,-4.15e-01),vec2(-1.26e-01,-4.03e-01),vec2(-1.02e-01,-3.92e-01),vec2(-1.14e-01,-3.92e-01),vec2(-1.14e-01,-4.03e-01),vec2(-1.14e-01,-4.03e-01),vec2(-1.14e-01,-4.15e-01),vec2(-1.02e-01,-4.15e-01),vec2(-1.02e-01,-4.15e-01),vec2(-9.08e-02,-4.15e-01),vec2(-9.08e-02,-4.03e-01),vec2(-1.02e-01,-3.92e-01),vec2(-9.08e-02,-3.92e-01),vec2(-9.08e-02,-4.03e-01),vec2(-7.86e-02,-4.03e-01),vec2(-7.86e-02,-3.92e-01),vec2(-6.69e-02,-3.92e-01),vec2(-4.30e-02,-4.03e-01),vec2(-3.13e-02,-4.03e-01),vec2(-3.13e-02,-3.92e-01),vec2(-4.30e-02,-4.03e-01),vec2(-3.13e-02,-4.03e-01),vec2(-3.13e-02,-4.15e-01),vec2(-7.42e-03,-3.92e-01),vec2(4.25e-03,-3.92e-01),vec2(4.25e-03,-4.03e-01),vec2(-1.91e-02,-4.03e-01),vec2(-1.91e-02,-3.92e-01),vec2(-7.42e-03,-3.92e-01),vec2(3.98e-02,-4.03e-01),vec2(3.98e-02,-3.92e-01),vec2(2.82e-02,-3.92e-01),vec2(2.82e-02,-3.92e-01),vec2(1.65e-02,-3.92e-01),vec2(1.65e-02,-4.03e-01),vec2(1.65e-02,-4.03e-01),vec2(1.65e-02,-4.15e-01),vec2(2.82e-02,-4.15e-01),vec2(6.38e-02,-4.15e-01),vec2(8.71e-02,-4.09e-01),vec2(6.38e-02,-4.03e-01),vec2(6.38e-02,-4.03e-01),vec2(4.04e-02,-3.98e-01),vec2(6.38e-02,-3.92e-01),vec2(9.93e-02,-4.15e-01),vec2(1.23e-01,-4.09e-01),vec2(9.93e-02,-4.03e-01),vec2(9.93e-02,-4.03e-01),vec2(7.60e-02,-3.98e-01),vec2(9.93e-02,-3.92e-01));
        for(int i=0; i<29;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
        for(int i=0; i<30; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
		col = mix(col, .8*c.xxx, B(10.)*smoothstep(.005, .002, d ));
    }
    else if(iTime < 18.)
    {
        const vec2 lin[18] = vec2[18](vec2(-4.00e-01,-3.68e-01),vec2(-4.00e-01,-4.03e-01),vec2(-3.77e-01,-3.68e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.53e-01,-3.68e-01),vec2(-3.53e-01,-4.03e-01),vec2(-3.29e-01,-4.15e-01),vec2(-3.18e-01,-4.15e-01),vec2(-3.41e-01,-4.03e-01),vec2(-3.18e-01,-4.03e-01),vec2(-2.58e-01,-4.15e-01),vec2(-2.58e-01,-3.92e-01),vec2(-2.46e-01,-4.15e-01),vec2(-2.46e-01,-3.92e-01),vec2(-2.10e-01,-4.15e-01),vec2(-1.99e-01,-4.15e-01),vec2(-2.22e-01,-4.03e-01),vec2(-1.99e-01,-4.03e-01)),
        quad[45] = vec2[45](vec2(-4.00e-01,-4.03e-01),vec2(-4.00e-01,-4.15e-01),vec2(-3.88e-01,-4.15e-01),vec2(-3.88e-01,-4.15e-01),vec2(-3.77e-01,-4.15e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.77e-01,-4.15e-01),vec2(-3.65e-01,-4.15e-01),vec2(-3.65e-01,-4.15e-01),vec2(-3.53e-01,-4.15e-01),vec2(-3.53e-01,-4.03e-01),vec2(-3.18e-01,-4.03e-01),vec2(-3.18e-01,-3.92e-01),vec2(-3.29e-01,-3.92e-01),vec2(-3.29e-01,-3.92e-01),vec2(-3.41e-01,-3.92e-01),vec2(-3.41e-01,-4.03e-01),vec2(-3.41e-01,-4.03e-01),vec2(-3.41e-01,-4.15e-01),vec2(-3.29e-01,-4.15e-01),vec2(-2.70e-01,-3.92e-01),vec2(-2.82e-01,-3.92e-01),vec2(-2.82e-01,-4.03e-01),vec2(-2.82e-01,-4.03e-01),vec2(-2.82e-01,-4.15e-01),vec2(-2.70e-01,-4.15e-01),vec2(-2.70e-01,-4.15e-01),vec2(-2.58e-01,-4.15e-01),vec2(-2.58e-01,-4.03e-01),vec2(-2.70e-01,-3.92e-01),vec2(-2.58e-01,-3.92e-01),vec2(-2.58e-01,-4.03e-01),vec2(-2.46e-01,-4.03e-01),vec2(-2.46e-01,-3.92e-01),vec2(-2.34e-01,-3.92e-01),vec2(-1.99e-01,-4.03e-01),vec2(-1.99e-01,-3.92e-01),vec2(-2.10e-01,-3.92e-01),vec2(-2.10e-01,-3.92e-01),vec2(-2.22e-01,-3.92e-01),vec2(-2.22e-01,-4.03e-01),vec2(-2.22e-01,-4.03e-01),vec2(-2.22e-01,-4.15e-01),vec2(-2.10e-01,-4.15e-01));
        for(int i=0; i<9;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
        for(int i=0; i<15; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
        col = mix(col, .8*c.xxx, B(15.)*smoothstep(.005, .002, d ));
    }
    else if(iTime < 23.)
    {
        const vec2 lin[128] = vec2[128](vec2(-4.00e-01,-4.03e-01),vec2(-4.00e-01,-3.80e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.77e-01,-3.80e-01),vec2(-3.64e-01,-4.15e-01),vec2(-3.64e-01,-3.68e-01),vec2(-3.18e-01,-4.15e-01),vec2(-3.18e-01,-3.68e-01),vec2(-3.64e-01,-3.68e-01),vec2(-3.41e-01,-3.92e-01),vec2(-3.41e-01,-3.92e-01),vec2(-3.18e-01,-3.68e-01),vec2(-2.82e-01,-4.03e-01),vec2(-2.82e-01,-4.03e-01),vec2(-2.82e-01,-3.80e-01),vec2(-2.82e-01,-3.80e-01),vec2(-2.69e-01,-4.03e-01),vec2(-2.69e-01,-4.03e-01),vec2(-2.69e-01,-3.80e-01),vec2(-2.69e-01,-3.80e-01),vec2(-2.10e-01,-4.15e-01),vec2(-2.21e-01,-4.15e-01),vec2(-2.10e-01,-3.68e-01),vec2(-2.21e-01,-3.68e-01),vec2(-2.33e-01,-4.03e-01),vec2(-2.33e-01,-3.80e-01),vec2(-1.50e-01,-3.92e-01),vec2(-1.39e-01,-3.92e-01),vec2(-1.50e-01,-4.15e-01),vec2(-1.39e-01,-4.15e-01),vec2(-1.39e-01,-4.15e-01),vec2(-1.39e-01,-3.68e-01),vec2(-1.15e-01,-4.15e-01),vec2(-1.03e-01,-4.15e-01),vec2(-1.26e-01,-4.03e-01),vec2(-1.03e-01,-4.03e-01),vec2(-9.08e-02,-4.15e-01),vec2(-6.75e-02,-3.68e-01),vec2(-5.52e-02,-4.15e-01),vec2(-4.36e-02,-4.15e-01),vec2(-3.19e-02,-3.68e-01),vec2(-4.36e-02,-3.68e-01),vec2(-1.97e-02,-4.15e-01),vec2(-1.97e-02,-3.68e-01),vec2(-1.97e-02,-3.68e-01),vec2(3.67e-03,-3.68e-01),vec2(-1.97e-02,-3.92e-01),vec2(3.67e-03,-3.92e-01),vec2(1.59e-02,-4.15e-01),vec2(3.93e-02,-3.68e-01),vec2(1.59e-02,-3.68e-01),vec2(3.93e-02,-4.15e-01),vec2(9.88e-02,-4.15e-01),vec2(9.88e-02,-3.92e-01),vec2(1.11e-01,-4.15e-01),vec2(1.11e-01,-3.92e-01),vec2(1.34e-01,-4.03e-01),vec2(1.34e-01,-4.15e-01),vec2(1.58e-01,-3.92e-01),vec2(1.70e-01,-3.92e-01),vec2(1.58e-01,-4.15e-01),vec2(1.70e-01,-4.15e-01),vec2(1.70e-01,-4.15e-01),vec2(1.70e-01,-3.68e-01),vec2(2.06e-01,-4.15e-01),vec2(2.06e-01,-3.68e-01),vec2(2.06e-01,-3.68e-01),vec2(2.29e-01,-4.15e-01),vec2(2.29e-01,-4.15e-01),vec2(2.29e-01,-3.68e-01),vec2(2.42e-01,-4.15e-01),vec2(2.42e-01,-3.68e-01),vec2(2.42e-01,-3.68e-01),vec2(2.53e-01,-3.68e-01),vec2(2.42e-01,-3.92e-01),vec2(2.53e-01,-3.92e-01),vec2(2.65e-01,-4.03e-01),vec2(2.65e-01,-4.15e-01),vec2(3.01e-01,-3.68e-01),vec2(2.77e-01,-4.03e-01),vec2(2.77e-01,-4.03e-01),vec2(3.01e-01,-4.03e-01),vec2(3.01e-01,-3.92e-01),vec2(3.01e-01,-4.15e-01),vec2(3.37e-01,-4.03e-01),vec2(3.37e-01,-4.03e-01),vec2(3.37e-01,-3.80e-01),vec2(3.37e-01,-3.80e-01),vec2(3.49e-01,-4.03e-01),vec2(3.49e-01,-4.03e-01),vec2(3.49e-01,-3.80e-01),vec2(3.49e-01,-3.80e-01),vec2(4.09e-01,-4.15e-01),vec2(3.97e-01,-4.15e-01),vec2(4.09e-01,-3.68e-01),vec2(3.97e-01,-3.68e-01),vec2(3.85e-01,-4.03e-01),vec2(3.85e-01,-3.80e-01),vec2(4.68e-01,-3.92e-01),vec2(4.80e-01,-3.92e-01),vec2(4.68e-01,-4.15e-01),vec2(4.80e-01,-4.15e-01),vec2(4.80e-01,-4.15e-01),vec2(4.80e-01,-3.68e-01),vec2(5.04e-01,-4.15e-01),vec2(5.15e-01,-4.15e-01),vec2(4.92e-01,-4.03e-01),vec2(5.15e-01,-4.03e-01),vec2(5.27e-01,-4.15e-01),vec2(5.51e-01,-3.68e-01),vec2(5.75e-01,-3.92e-01),vec2(5.86e-01,-3.92e-01),vec2(5.86e-01,-3.92e-01),vec2(5.86e-01,-4.03e-01),vec2(5.63e-01,-4.03e-01),vec2(5.63e-01,-3.80e-01),vec2(5.75e-01,-3.68e-01),vec2(5.86e-01,-3.68e-01),vec2(5.99e-01,-4.15e-01),vec2(5.99e-01,-3.68e-01),vec2(5.99e-01,-3.68e-01),vec2(6.22e-01,-3.68e-01),vec2(5.99e-01,-3.92e-01),vec2(6.22e-01,-3.92e-01),vec2(6.34e-01,-4.15e-01),vec2(6.58e-01,-3.68e-01),vec2(6.34e-01,-3.68e-01),vec2(6.58e-01,-4.15e-01)),
        quad[135] = vec2[135](vec2(-4.00e-01,-4.03e-01),vec2(-4.00e-01,-4.15e-01),vec2(-3.88e-01,-4.15e-01),vec2(-3.88e-01,-4.15e-01),vec2(-3.77e-01,-4.15e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.77e-01,-3.80e-01),vec2(-3.77e-01,-3.68e-01),vec2(-3.88e-01,-3.68e-01),vec2(-3.88e-01,-3.68e-01),vec2(-4.00e-01,-3.68e-01),vec2(-4.00e-01,-3.80e-01),vec2(-3.88e-01,-4.03e-01),vec2(-3.88e-01,-4.15e-01),vec2(-3.77e-01,-4.15e-01),vec2(-2.21e-01,-4.15e-01),vec2(-2.33e-01,-4.15e-01),vec2(-2.33e-01,-4.03e-01),vec2(-2.33e-01,-3.80e-01),vec2(-2.33e-01,-3.68e-01),vec2(-2.21e-01,-3.68e-01),vec2(-1.86e-01,-4.15e-01),vec2(-1.98e-01,-4.15e-01),vec2(-1.98e-01,-4.03e-01),vec2(-1.98e-01,-4.03e-01),vec2(-1.98e-01,-3.92e-01),vec2(-1.86e-01,-3.92e-01),vec2(-1.86e-01,-3.92e-01),vec2(-1.74e-01,-3.92e-01),vec2(-1.74e-01,-4.03e-01),vec2(-1.74e-01,-4.03e-01),vec2(-1.74e-01,-4.15e-01),vec2(-1.86e-01,-4.15e-01),vec2(-1.50e-01,-4.15e-01),vec2(-1.62e-01,-4.15e-01),vec2(-1.62e-01,-4.03e-01),vec2(-1.62e-01,-4.03e-01),vec2(-1.62e-01,-3.92e-01),vec2(-1.50e-01,-3.92e-01),vec2(-1.03e-01,-4.03e-01),vec2(-1.03e-01,-3.92e-01),vec2(-1.15e-01,-3.92e-01),vec2(-1.15e-01,-3.92e-01),vec2(-1.26e-01,-3.92e-01),vec2(-1.26e-01,-4.03e-01),vec2(-1.26e-01,-4.03e-01),vec2(-1.26e-01,-4.15e-01),vec2(-1.15e-01,-4.15e-01),vec2(-4.36e-02,-4.15e-01),vec2(-3.19e-02,-4.15e-01),vec2(-3.19e-02,-4.03e-01),vec2(-3.19e-02,-4.03e-01),vec2(-3.19e-02,-3.92e-01),vec2(-4.36e-02,-3.92e-01),vec2(-4.36e-02,-3.92e-01),vec2(-5.52e-02,-3.92e-01),vec2(-5.52e-02,-3.80e-01),vec2(-5.52e-02,-3.80e-01),vec2(-5.52e-02,-3.68e-01),vec2(-4.36e-02,-3.68e-01),vec2(8.71e-02,-3.92e-01),vec2(7.54e-02,-3.92e-01),vec2(7.54e-02,-4.03e-01),vec2(7.54e-02,-4.03e-01),vec2(7.54e-02,-4.15e-01),vec2(8.71e-02,-4.15e-01),vec2(8.71e-02,-4.15e-01),vec2(9.88e-02,-4.15e-01),vec2(9.88e-02,-4.03e-01),vec2(8.71e-02,-3.92e-01),vec2(9.88e-02,-3.92e-01),vec2(9.88e-02,-4.03e-01),vec2(1.23e-01,-3.92e-01),vec2(1.34e-01,-3.92e-01),vec2(1.34e-01,-4.03e-01),vec2(1.11e-01,-4.03e-01),vec2(1.11e-01,-3.92e-01),vec2(1.23e-01,-3.92e-01),vec2(1.58e-01,-4.15e-01),vec2(1.47e-01,-4.15e-01),vec2(1.47e-01,-4.03e-01),vec2(1.47e-01,-4.03e-01),vec2(1.47e-01,-3.92e-01),vec2(1.58e-01,-3.92e-01),vec2(2.53e-01,-3.68e-01),vec2(2.65e-01,-3.68e-01),vec2(2.65e-01,-3.80e-01),vec2(2.65e-01,-3.80e-01),vec2(2.65e-01,-3.92e-01),vec2(2.53e-01,-3.92e-01),vec2(2.53e-01,-3.92e-01),vec2(2.65e-01,-3.92e-01),vec2(2.65e-01,-4.03e-01),vec2(3.97e-01,-4.15e-01),vec2(3.85e-01,-4.15e-01),vec2(3.85e-01,-4.03e-01),vec2(3.85e-01,-3.80e-01),vec2(3.85e-01,-3.68e-01),vec2(3.97e-01,-3.68e-01),vec2(4.32e-01,-4.15e-01),vec2(4.21e-01,-4.15e-01),vec2(4.21e-01,-4.03e-01),vec2(4.21e-01,-4.03e-01),vec2(4.21e-01,-3.92e-01),vec2(4.32e-01,-3.92e-01),vec2(4.32e-01,-3.92e-01),vec2(4.44e-01,-3.92e-01),vec2(4.44e-01,-4.03e-01),vec2(4.44e-01,-4.03e-01),vec2(4.44e-01,-4.15e-01),vec2(4.32e-01,-4.15e-01),vec2(4.68e-01,-4.15e-01),vec2(4.56e-01,-4.15e-01),vec2(4.56e-01,-4.03e-01),vec2(4.56e-01,-4.03e-01),vec2(4.56e-01,-3.92e-01),vec2(4.68e-01,-3.92e-01),vec2(5.15e-01,-4.03e-01),vec2(5.15e-01,-3.92e-01),vec2(5.04e-01,-3.92e-01),vec2(5.04e-01,-3.92e-01),vec2(4.92e-01,-3.92e-01),vec2(4.92e-01,-4.03e-01),vec2(4.92e-01,-4.03e-01),vec2(4.92e-01,-4.15e-01),vec2(5.04e-01,-4.15e-01),vec2(5.63e-01,-3.80e-01),vec2(5.63e-01,-3.68e-01),vec2(5.75e-01,-3.68e-01),vec2(5.75e-01,-4.15e-01),vec2(5.63e-01,-4.15e-01),vec2(5.63e-01,-4.03e-01),vec2(5.75e-01,-4.15e-01),vec2(5.86e-01,-4.15e-01),vec2(5.86e-01,-4.03e-01));
        for(int i=0; i<64;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
        for(int i=0; i<45; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
		col = mix(col, .8*c.xxx, B(20.)*smoothstep(.005, .002, d ));
    }
    else if(iTime < 28.)
    {
        const vec2 lin[24] = vec2[24](vec2(-4.00e-01,-3.68e-01),vec2(-4.00e-01,-4.03e-01),vec2(-3.77e-01,-3.68e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.53e-01,-3.68e-01),vec2(-3.53e-01,-4.03e-01),vec2(-3.29e-01,-4.15e-01),vec2(-3.18e-01,-4.15e-01),vec2(-3.41e-01,-4.03e-01),vec2(-3.18e-01,-4.03e-01),vec2(-2.82e-01,-4.03e-01),vec2(-2.82e-01,-3.68e-01),vec2(-2.58e-01,-4.03e-01),vec2(-2.58e-01,-3.92e-01),vec2(-2.58e-01,-3.80e-01),vec2(-2.58e-01,-3.80e-01),vec2(-2.34e-01,-4.15e-01),vec2(-2.34e-01,-3.68e-01),vec2(-2.34e-01,-4.03e-01),vec2(-2.22e-01,-4.03e-01),vec2(-1.86e-01,-4.15e-01),vec2(-1.75e-01,-4.15e-01),vec2(-1.98e-01,-4.03e-01),vec2(-1.75e-01,-4.03e-01)),
        quad[42] = vec2[42](vec2(-4.00e-01,-4.03e-01),vec2(-4.00e-01,-4.15e-01),vec2(-3.88e-01,-4.15e-01),vec2(-3.88e-01,-4.15e-01),vec2(-3.77e-01,-4.15e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.77e-01,-4.03e-01),vec2(-3.77e-01,-4.15e-01),vec2(-3.65e-01,-4.15e-01),vec2(-3.65e-01,-4.15e-01),vec2(-3.53e-01,-4.15e-01),vec2(-3.53e-01,-4.03e-01),vec2(-3.18e-01,-4.03e-01),vec2(-3.18e-01,-3.92e-01),vec2(-3.29e-01,-3.92e-01),vec2(-3.29e-01,-3.92e-01),vec2(-3.41e-01,-3.92e-01),vec2(-3.41e-01,-4.03e-01),vec2(-3.41e-01,-4.03e-01),vec2(-3.41e-01,-4.15e-01),vec2(-3.29e-01,-4.15e-01),vec2(-2.82e-01,-4.03e-01),vec2(-2.82e-01,-4.15e-01),vec2(-2.70e-01,-4.15e-01),vec2(-2.58e-01,-4.03e-01),vec2(-2.58e-01,-4.15e-01),vec2(-2.46e-01,-4.15e-01),vec2(-2.22e-01,-4.03e-01),vec2(-2.10e-01,-4.03e-01),vec2(-2.10e-01,-3.92e-01),vec2(-2.22e-01,-4.03e-01),vec2(-2.10e-01,-4.03e-01),vec2(-2.10e-01,-4.15e-01),vec2(-1.75e-01,-4.03e-01),vec2(-1.75e-01,-3.92e-01),vec2(-1.86e-01,-3.92e-01),vec2(-1.86e-01,-3.92e-01),vec2(-1.98e-01,-3.92e-01),vec2(-1.98e-01,-4.03e-01),vec2(-1.98e-01,-4.03e-01),vec2(-1.98e-01,-4.15e-01),vec2(-1.86e-01,-4.15e-01));
        for(int i=0; i<12;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
        for(int i=0; i<14; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
		col = mix(col, .8*c.xxx, B(25.)*smoothstep(.005, .002, d ));
    }
    else if(iTime < 33.)
    {
        const vec2 lin[176] = vec2[176](vec2(-4.00e-01,-4.15e-01),vec2(-4.00e-01,-3.80e-01),vec2(-3.77e-01,-4.15e-01),vec2(-3.77e-01,-3.80e-01),vec2(-4.00e-01,-3.92e-01),vec2(-3.77e-01,-3.92e-01),vec2(-3.64e-01,-4.15e-01),vec2(-3.64e-01,-3.68e-01),vec2(-3.64e-01,-3.68e-01),vec2(-3.41e-01,-4.15e-01),vec2(-3.41e-01,-4.15e-01),vec2(-3.41e-01,-3.68e-01),vec2(-3.29e-01,-4.15e-01),vec2(-3.29e-01,-3.68e-01),vec2(-3.29e-01,-3.68e-01),vec2(-3.17e-01,-3.68e-01),vec2(-3.29e-01,-4.15e-01),vec2(-3.17e-01,-4.15e-01),vec2(-3.05e-01,-4.03e-01),vec2(-3.05e-01,-3.80e-01),vec2(-2.69e-01,-4.03e-01),vec2(-2.69e-01,-4.03e-01),vec2(-2.69e-01,-3.80e-01),vec2(-2.69e-01,-3.80e-01),vec2(-2.57e-01,-4.03e-01),vec2(-2.57e-01,-4.03e-01),vec2(-2.57e-01,-3.80e-01),vec2(-2.57e-01,-3.80e-01),vec2(-2.21e-01,-4.15e-01),vec2(-2.21e-01,-3.92e-01),vec2(-1.74e-01,-3.92e-01),vec2(-1.74e-01,-4.27e-01),vec2(-1.85e-01,-4.38e-01),vec2(-1.97e-01,-4.38e-01),vec2(-1.61e-01,-4.15e-01),vec2(-1.61e-01,-3.68e-01),vec2(-1.61e-01,-4.15e-01),vec2(-1.50e-01,-4.15e-01),vec2(-1.61e-01,-3.92e-01),vec2(-1.50e-01,-3.92e-01),vec2(-1.02e-01,-4.15e-01),vec2(-1.02e-01,-3.92e-01),vec2(-6.63e-02,-4.03e-01),vec2(-6.63e-02,-4.03e-01),vec2(-6.63e-02,-3.80e-01),vec2(-6.63e-02,-3.80e-01),vec2(-5.41e-02,-4.03e-01),vec2(-5.41e-02,-4.03e-01),vec2(-5.41e-02,-3.80e-01),vec2(-5.41e-02,-3.80e-01),vec2(-6.25e-03,-4.03e-01),vec2(-6.25e-03,-3.80e-01),vec2(-1.79e-02,-3.92e-01),vec2(5.42e-03,-3.92e-01),vec2(4.10e-02,-4.15e-01),vec2(4.10e-02,-3.92e-01),vec2(5.33e-02,-4.15e-01),vec2(5.33e-02,-3.92e-01),vec2(7.72e-02,-4.15e-01),vec2(7.72e-02,-3.68e-01),vec2(7.72e-02,-4.15e-01),vec2(8.88e-02,-4.15e-01),vec2(7.72e-02,-3.92e-01),vec2(8.88e-02,-3.92e-01),vec2(1.13e-01,-4.15e-01),vec2(1.13e-01,-3.92e-01),vec2(1.60e-01,-4.15e-01),vec2(1.60e-01,-3.92e-01),vec2(1.72e-01,-3.92e-01),vec2(1.72e-01,-4.03e-01),vec2(1.96e-01,-3.92e-01),vec2(1.96e-01,-4.15e-01),vec2(2.08e-01,-4.15e-01),vec2(2.20e-01,-4.15e-01),vec2(2.20e-01,-3.92e-01),vec2(2.31e-01,-3.92e-01),vec2(2.67e-01,-3.92e-01),vec2(2.55e-01,-3.92e-01),vec2(2.67e-01,-4.15e-01),vec2(2.55e-01,-4.15e-01),vec2(2.79e-01,-4.15e-01),vec2(2.79e-01,-3.68e-01),vec2(3.02e-01,-4.15e-01),vec2(3.02e-01,-4.03e-01),vec2(2.91e-01,-3.92e-01),vec2(2.79e-01,-3.92e-01),vec2(3.39e-01,-4.03e-01),vec2(3.39e-01,-4.03e-01),vec2(3.39e-01,-3.80e-01),vec2(3.39e-01,-3.80e-01),vec2(3.51e-01,-4.03e-01),vec2(3.51e-01,-4.03e-01),vec2(3.51e-01,-3.80e-01),vec2(3.51e-01,-3.80e-01),vec2(3.87e-01,-4.15e-01),vec2(3.87e-01,-3.92e-01),vec2(4.10e-01,-4.15e-01),vec2(4.10e-01,-4.03e-01),vec2(4.34e-01,-4.15e-01),vec2(4.34e-01,-4.03e-01),vec2(4.58e-01,-4.15e-01),vec2(4.69e-01,-4.15e-01),vec2(4.46e-01,-4.03e-01),vec2(4.69e-01,-4.03e-01),vec2(4.81e-01,-4.15e-01),vec2(4.81e-01,-3.92e-01),vec2(5.29e-01,-3.92e-01),vec2(5.17e-01,-3.92e-01),vec2(5.29e-01,-4.15e-01),vec2(5.17e-01,-4.15e-01),vec2(5.41e-01,-3.92e-01),vec2(5.41e-01,-4.03e-01),vec2(5.64e-01,-3.92e-01),vec2(5.64e-01,-4.15e-01),vec2(5.77e-01,-4.15e-01),vec2(5.77e-01,-3.92e-01),vec2(6.00e-01,-3.92e-01),vec2(6.00e-01,-4.03e-01),vec2(6.24e-01,-3.92e-01),vec2(6.24e-01,-4.27e-01),vec2(6.12e-01,-4.38e-01),vec2(6.00e-01,-4.38e-01),vec2(6.60e-01,-4.03e-01),vec2(6.60e-01,-4.03e-01),vec2(6.60e-01,-3.80e-01),vec2(6.60e-01,-3.80e-01),vec2(6.72e-01,-4.03e-01),vec2(6.72e-01,-4.03e-01),vec2(6.72e-01,-3.80e-01),vec2(6.72e-01,-3.80e-01),vec2(7.32e-01,-4.15e-01),vec2(7.20e-01,-4.15e-01),vec2(7.32e-01,-3.68e-01),vec2(7.20e-01,-3.68e-01),vec2(7.08e-01,-4.03e-01),vec2(7.08e-01,-3.80e-01),vec2(7.56e-01,-4.03e-01),vec2(7.56e-01,-3.68e-01),vec2(7.44e-01,-3.92e-01),vec2(7.67e-01,-3.92e-01),vec2(7.79e-01,-4.15e-01),vec2(7.79e-01,-3.92e-01),vec2(8.03e-01,-4.03e-01),vec2(8.03e-01,-3.68e-01),vec2(8.27e-01,-3.92e-01),vec2(8.51e-01,-3.92e-01),vec2(8.63e-01,-4.15e-01),vec2(8.63e-01,-3.80e-01),vec2(8.86e-01,-4.15e-01),vec2(8.86e-01,-3.80e-01),vec2(8.63e-01,-3.92e-01),vec2(8.86e-01,-3.92e-01),vec2(8.98e-01,-4.03e-01),vec2(8.98e-01,-3.68e-01),vec2(9.34e-01,-4.03e-01),vec2(9.34e-01,-3.68e-01),vec2(9.22e-01,-3.92e-01),vec2(9.46e-01,-3.92e-01),vec2(9.58e-01,-3.92e-01),vec2(9.81e-01,-3.92e-01),vec2(1.01e+00,-4.15e-01),vec2(1.01e+00,-3.68e-01),vec2(9.94e-01,-3.68e-01),vec2(1.02e+00,-3.68e-01),vec2(1.04e+00,-4.15e-01),vec2(1.05e+00,-4.15e-01),vec2(1.03e+00,-4.03e-01),vec2(1.05e+00,-4.03e-01),vec2(1.06e+00,-4.15e-01),vec2(1.08e+00,-4.15e-01),vec2(1.08e+00,-3.92e-01),vec2(1.09e+00,-3.92e-01),vec2(1.11e+00,-4.03e-01),vec2(1.11e+00,-3.68e-01),vec2(1.10e+00,-3.92e-01),vec2(1.12e+00,-3.92e-01)),
        quad[204] = vec2[204](vec2(-3.77e-01,-3.80e-01),vec2(-3.77e-01,-3.68e-01),vec2(-3.88e-01,-3.68e-01),vec2(-3.88e-01,-3.68e-01),vec2(-4.00e-01,-3.68e-01),vec2(-4.00e-01,-3.80e-01),vec2(-3.17e-01,-3.68e-01),vec2(-3.05e-01,-3.68e-01),vec2(-3.05e-01,-3.80e-01),vec2(-3.17e-01,-4.15e-01),vec2(-3.05e-01,-4.15e-01),vec2(-3.05e-01,-4.03e-01),vec2(-2.21e-01,-4.03e-01),vec2(-2.21e-01,-3.92e-01),vec2(-2.09e-01,-3.92e-01),vec2(-1.74e-01,-4.03e-01),vec2(-1.74e-01,-3.92e-01),vec2(-1.85e-01,-3.92e-01),vec2(-1.85e-01,-3.92e-01),vec2(-1.97e-01,-3.92e-01),vec2(-1.97e-01,-4.03e-01),vec2(-1.97e-01,-4.03e-01),vec2(-1.97e-01,-4.15e-01),vec2(-1.85e-01,-4.15e-01),vec2(-1.85e-01,-4.15e-01),vec2(-1.74e-01,-4.15e-01),vec2(-1.74e-01,-4.03e-01),vec2(-1.74e-01,-4.27e-01),vec2(-1.74e-01,-4.38e-01),vec2(-1.85e-01,-4.38e-01),vec2(-1.50e-01,-4.15e-01),vec2(-1.38e-01,-4.15e-01),vec2(-1.38e-01,-4.03e-01),vec2(-1.38e-01,-4.03e-01),vec2(-1.38e-01,-3.92e-01),vec2(-1.50e-01,-3.92e-01),vec2(-1.14e-01,-3.92e-01),vec2(-1.26e-01,-3.92e-01),vec2(-1.26e-01,-4.03e-01),vec2(-1.26e-01,-4.03e-01),vec2(-1.26e-01,-4.15e-01),vec2(-1.14e-01,-4.15e-01),vec2(-1.14e-01,-4.15e-01),vec2(-1.02e-01,-4.15e-01),vec2(-1.02e-01,-4.03e-01),vec2(-1.14e-01,-3.92e-01),vec2(-1.02e-01,-3.92e-01),vec2(-1.02e-01,-4.03e-01),vec2(-1.79e-02,-4.15e-01),vec2(-6.25e-03,-4.15e-01),vec2(-6.25e-03,-4.03e-01),vec2(-6.25e-03,-3.80e-01),vec2(-6.25e-03,-3.68e-01),vec2(5.42e-03,-3.68e-01),vec2(2.93e-02,-3.92e-01),vec2(1.77e-02,-3.92e-01),vec2(1.77e-02,-4.03e-01),vec2(1.77e-02,-4.03e-01),vec2(1.77e-02,-4.15e-01),vec2(2.93e-02,-4.15e-01),vec2(2.93e-02,-4.15e-01),vec2(4.10e-02,-4.15e-01),vec2(4.10e-02,-4.03e-01),vec2(2.93e-02,-3.92e-01),vec2(4.10e-02,-3.92e-01),vec2(4.10e-02,-4.03e-01),vec2(5.33e-02,-4.03e-01),vec2(5.33e-02,-3.92e-01),vec2(6.49e-02,-3.92e-01),vec2(8.88e-02,-4.15e-01),vec2(1.01e-01,-4.15e-01),vec2(1.01e-01,-4.03e-01),vec2(1.01e-01,-4.03e-01),vec2(1.01e-01,-3.92e-01),vec2(8.88e-02,-3.92e-01),vec2(1.13e-01,-4.03e-01),vec2(1.13e-01,-3.92e-01),vec2(1.24e-01,-3.92e-01),vec2(1.48e-01,-3.92e-01),vec2(1.37e-01,-3.92e-01),vec2(1.37e-01,-4.03e-01),vec2(1.37e-01,-4.03e-01),vec2(1.37e-01,-4.15e-01),vec2(1.48e-01,-4.15e-01),vec2(1.48e-01,-4.15e-01),vec2(1.60e-01,-4.15e-01),vec2(1.60e-01,-4.03e-01),vec2(1.48e-01,-3.92e-01),vec2(1.60e-01,-3.92e-01),vec2(1.60e-01,-4.03e-01),vec2(1.72e-01,-4.03e-01),vec2(1.72e-01,-4.15e-01),vec2(1.84e-01,-4.15e-01),vec2(1.84e-01,-4.15e-01),vec2(1.96e-01,-4.15e-01),vec2(1.96e-01,-4.03e-01),vec2(2.20e-01,-4.15e-01),vec2(2.43e-01,-4.09e-01),vec2(2.20e-01,-4.03e-01),vec2(2.20e-01,-4.03e-01),vec2(1.96e-01,-3.98e-01),vec2(2.20e-01,-3.92e-01),vec2(2.55e-01,-3.92e-01),vec2(2.43e-01,-3.92e-01),vec2(2.43e-01,-4.03e-01),vec2(2.43e-01,-4.03e-01),vec2(2.43e-01,-4.15e-01),vec2(2.55e-01,-4.15e-01),vec2(3.02e-01,-4.03e-01),vec2(3.02e-01,-3.92e-01),vec2(2.91e-01,-3.92e-01),vec2(3.87e-01,-4.03e-01),vec2(3.87e-01,-3.92e-01),vec2(3.99e-01,-3.92e-01),vec2(3.99e-01,-3.92e-01),vec2(4.10e-01,-3.92e-01),vec2(4.10e-01,-4.03e-01),vec2(4.10e-01,-4.03e-01),vec2(4.10e-01,-3.92e-01),vec2(4.22e-01,-3.92e-01),vec2(4.22e-01,-3.92e-01),vec2(4.34e-01,-3.92e-01),vec2(4.34e-01,-4.03e-01),vec2(4.69e-01,-4.03e-01),vec2(4.69e-01,-3.92e-01),vec2(4.58e-01,-3.92e-01),vec2(4.58e-01,-3.92e-01),vec2(4.46e-01,-3.92e-01),vec2(4.46e-01,-4.03e-01),vec2(4.46e-01,-4.03e-01),vec2(4.46e-01,-4.15e-01),vec2(4.58e-01,-4.15e-01),vec2(4.81e-01,-4.03e-01),vec2(4.81e-01,-3.92e-01),vec2(4.93e-01,-3.92e-01),vec2(5.17e-01,-3.92e-01),vec2(5.05e-01,-3.92e-01),vec2(5.05e-01,-4.03e-01),vec2(5.05e-01,-4.03e-01),vec2(5.05e-01,-4.15e-01),vec2(5.17e-01,-4.15e-01),vec2(5.41e-01,-4.03e-01),vec2(5.41e-01,-4.15e-01),vec2(5.53e-01,-4.15e-01),vec2(5.53e-01,-4.15e-01),vec2(5.64e-01,-4.15e-01),vec2(5.64e-01,-4.03e-01),vec2(5.77e-01,-4.03e-01),vec2(5.77e-01,-3.92e-01),vec2(5.88e-01,-3.92e-01),vec2(6.00e-01,-4.03e-01),vec2(6.00e-01,-4.15e-01),vec2(6.12e-01,-4.15e-01),vec2(6.12e-01,-4.15e-01),vec2(6.24e-01,-4.15e-01),vec2(6.24e-01,-4.03e-01),vec2(6.24e-01,-4.27e-01),vec2(6.24e-01,-4.38e-01),vec2(6.12e-01,-4.38e-01),vec2(7.20e-01,-4.15e-01),vec2(7.08e-01,-4.15e-01),vec2(7.08e-01,-4.03e-01),vec2(7.08e-01,-3.80e-01),vec2(7.08e-01,-3.68e-01),vec2(7.20e-01,-3.68e-01),vec2(7.56e-01,-4.03e-01),vec2(7.56e-01,-4.15e-01),vec2(7.67e-01,-4.15e-01),vec2(7.79e-01,-4.03e-01),vec2(7.79e-01,-3.92e-01),vec2(7.91e-01,-3.92e-01),vec2(8.03e-01,-4.03e-01),vec2(8.03e-01,-4.15e-01),vec2(8.15e-01,-4.15e-01),vec2(8.86e-01,-3.80e-01),vec2(8.86e-01,-3.68e-01),vec2(8.75e-01,-3.68e-01),vec2(8.75e-01,-3.68e-01),vec2(8.63e-01,-3.68e-01),vec2(8.63e-01,-3.80e-01),vec2(8.98e-01,-4.03e-01),vec2(8.98e-01,-4.15e-01),vec2(9.10e-01,-4.15e-01),vec2(9.34e-01,-4.03e-01),vec2(9.34e-01,-4.15e-01),vec2(9.46e-01,-4.15e-01),vec2(1.05e+00,-4.03e-01),vec2(1.05e+00,-3.92e-01),vec2(1.04e+00,-3.92e-01),vec2(1.04e+00,-3.92e-01),vec2(1.03e+00,-3.92e-01),vec2(1.03e+00,-4.03e-01),vec2(1.03e+00,-4.03e-01),vec2(1.03e+00,-4.15e-01),vec2(1.04e+00,-4.15e-01),vec2(1.08e+00,-4.15e-01),vec2(1.10e+00,-4.09e-01),vec2(1.08e+00,-4.03e-01),vec2(1.08e+00,-4.03e-01),vec2(1.05e+00,-3.98e-01),vec2(1.08e+00,-3.92e-01),vec2(1.11e+00,-4.03e-01),vec2(1.11e+00,-4.15e-01),vec2(1.12e+00,-4.15e-01));
        for(int i=0; i<88;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
        for(int i=0; i<68; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
		col = mix(col, .8*c.xxx, B(30.)*smoothstep(.005, .002, d ));
    }   
    
    if((iTime > 80.) && (iTime < 90.))
    {
        const vec2 lin[114] = vec2[114](vec2(1.67e-02,1.00e-01),vec2(8.33e-03,1.00e-01),vec2(1.67e-02,1.33e-01),vec2(8.33e-03,1.33e-01),vec2(0.00e+00,1.08e-01),vec2(0.00e+00,1.25e-01),vec2(2.54e-02,1.00e-01),vec2(2.54e-02,1.33e-01),vec2(4.21e-02,1.00e-01),vec2(4.21e-02,1.08e-01),vec2(3.37e-02,1.17e-01),vec2(2.54e-02,1.17e-01),vec2(5.92e-02,1.00e-01),vec2(6.75e-02,1.00e-01),vec2(5.08e-02,1.08e-01),vec2(6.75e-02,1.08e-01),vec2(8.46e-02,1.00e-01),vec2(9.29e-02,1.00e-01),vec2(7.62e-02,1.08e-01),vec2(9.29e-02,1.08e-01),vec2(1.02e-01,1.00e-01),vec2(1.02e-01,1.17e-01),vec2(1.19e-01,1.00e-01),vec2(1.27e-01,1.00e-01),vec2(1.27e-01,1.17e-01),vec2(1.35e-01,1.17e-01),vec2(1.70e-01,1.08e-01),vec2(1.70e-01,1.33e-01),vec2(1.61e-01,1.17e-01),vec2(1.78e-01,1.17e-01),vec2(2.29e-01,1.00e-01),vec2(2.29e-01,1.33e-01),vec2(2.63e-01,1.00e-01),vec2(2.63e-01,1.33e-01),vec2(2.29e-01,1.33e-01),vec2(2.46e-01,1.17e-01),vec2(2.46e-01,1.17e-01),vec2(2.63e-01,1.33e-01),vec2(2.80e-01,1.00e-01),vec2(2.80e-01,1.33e-01),vec2(2.71e-01,1.00e-01),vec2(2.88e-01,1.00e-01),vec2(2.71e-01,1.33e-01),vec2(2.88e-01,1.33e-01),vec2(3.13e-01,1.00e-01),vec2(3.05e-01,1.00e-01),vec2(3.13e-01,1.33e-01),vec2(3.05e-01,1.33e-01),vec2(2.97e-01,1.08e-01),vec2(2.97e-01,1.25e-01),vec2(3.56e-01,1.00e-01),vec2(3.56e-01,1.17e-01),vec2(3.65e-01,1.00e-01),vec2(3.65e-01,1.17e-01),vec2(3.81e-01,1.08e-01),vec2(3.81e-01,1.00e-01),vec2(3.98e-01,1.17e-01),vec2(4.07e-01,1.17e-01),vec2(3.98e-01,1.00e-01),vec2(4.07e-01,1.00e-01),vec2(4.07e-01,1.00e-01),vec2(4.07e-01,1.33e-01),vec2(4.32e-01,1.00e-01),vec2(4.41e-01,1.00e-01),vec2(4.49e-01,1.33e-01),vec2(4.41e-01,1.33e-01),vec2(4.75e-01,1.17e-01),vec2(4.66e-01,1.17e-01),vec2(4.75e-01,1.00e-01),vec2(4.66e-01,1.00e-01),vec2(4.83e-01,1.00e-01),vec2(4.83e-01,1.33e-01),vec2(5.00e-01,1.00e-01),vec2(5.00e-01,1.08e-01),vec2(4.92e-01,1.17e-01),vec2(4.83e-01,1.17e-01),vec2(5.09e-01,1.00e-01),vec2(5.09e-01,1.17e-01),vec2(5.25e-01,1.08e-01),vec2(5.25e-01,1.00e-01),vec2(5.51e-01,1.00e-01),vec2(5.51e-01,1.17e-01),vec2(5.60e-01,8.33e-02),vec2(5.60e-01,1.17e-01),vec2(5.60e-01,1.17e-01),vec2(5.68e-01,1.17e-01),vec2(5.60e-01,1.00e-01),vec2(5.68e-01,1.00e-01),vec2(5.85e-01,8.33e-02),vec2(5.85e-01,1.17e-01),vec2(5.85e-01,1.17e-01),vec2(5.93e-01,1.17e-01),vec2(5.85e-01,1.00e-01),vec2(5.93e-01,1.00e-01),vec2(6.10e-01,1.00e-01),vec2(6.19e-01,1.00e-01),vec2(6.19e-01,1.17e-01),vec2(6.27e-01,1.17e-01),vec2(6.52e-01,1.17e-01),vec2(6.52e-01,9.17e-02),vec2(6.44e-01,8.33e-02),vec2(6.36e-01,8.33e-02),vec2(6.61e-01,1.08e-01),vec2(6.61e-01,1.17e-01),vec2(6.61e-01,1.25e-01),vec2(6.61e-01,1.25e-01),vec2(6.78e-01,1.00e-01),vec2(6.78e-01,1.17e-01),vec2(6.95e-01,1.08e-01),vec2(6.95e-01,1.33e-01),vec2(7.12e-01,1.00e-01),vec2(7.21e-01,1.00e-01),vec2(7.21e-01,1.17e-01),vec2(7.29e-01,1.17e-01)),
		quad[168] = vec2[168](vec2(8.33e-03,1.00e-01),vec2(0.00e+00,1.00e-01),vec2(0.00e+00,1.08e-01),vec2(0.00e+00,1.25e-01),vec2(0.00e+00,1.33e-01),vec2(8.33e-03,1.33e-01),vec2(4.21e-02,1.08e-01),vec2(4.21e-02,1.17e-01),vec2(3.37e-02,1.17e-01),vec2(6.75e-02,1.08e-01),vec2(6.75e-02,1.17e-01),vec2(5.92e-02,1.17e-01),vec2(5.92e-02,1.17e-01),vec2(5.08e-02,1.17e-01),vec2(5.08e-02,1.08e-01),vec2(5.08e-02,1.08e-01),vec2(5.08e-02,1.00e-01),vec2(5.92e-02,1.00e-01),vec2(9.29e-02,1.08e-01),vec2(9.29e-02,1.17e-01),vec2(8.46e-02,1.17e-01),vec2(8.46e-02,1.17e-01),vec2(7.62e-02,1.17e-01),vec2(7.62e-02,1.08e-01),vec2(7.62e-02,1.08e-01),vec2(7.62e-02,1.00e-01),vec2(8.46e-02,1.00e-01),vec2(1.02e-01,1.08e-01),vec2(1.02e-01,1.17e-01),vec2(1.10e-01,1.17e-01),vec2(1.27e-01,1.00e-01),vec2(1.44e-01,1.04e-01),vec2(1.27e-01,1.08e-01),vec2(1.27e-01,1.08e-01),vec2(1.10e-01,1.13e-01),vec2(1.27e-01,1.17e-01),vec2(1.70e-01,1.08e-01),vec2(1.70e-01,1.00e-01),vec2(1.78e-01,1.00e-01),vec2(1.95e-01,1.00e-01),vec2(1.87e-01,1.00e-01),vec2(1.87e-01,1.08e-01),vec2(1.87e-01,1.08e-01),vec2(1.87e-01,1.17e-01),vec2(1.95e-01,1.17e-01),vec2(1.95e-01,1.17e-01),vec2(2.03e-01,1.17e-01),vec2(2.03e-01,1.08e-01),vec2(2.03e-01,1.08e-01),vec2(2.03e-01,1.00e-01),vec2(1.95e-01,1.00e-01),vec2(3.05e-01,1.00e-01),vec2(2.97e-01,1.00e-01),vec2(2.97e-01,1.08e-01),vec2(2.97e-01,1.25e-01),vec2(2.97e-01,1.33e-01),vec2(3.05e-01,1.33e-01),vec2(3.47e-01,1.17e-01),vec2(3.39e-01,1.17e-01),vec2(3.39e-01,1.08e-01),vec2(3.39e-01,1.08e-01),vec2(3.39e-01,1.00e-01),vec2(3.47e-01,1.00e-01),vec2(3.47e-01,1.00e-01),vec2(3.56e-01,1.00e-01),vec2(3.56e-01,1.08e-01),vec2(3.47e-01,1.17e-01),vec2(3.56e-01,1.17e-01),vec2(3.56e-01,1.08e-01),vec2(3.73e-01,1.17e-01),vec2(3.81e-01,1.17e-01),vec2(3.81e-01,1.08e-01),vec2(3.65e-01,1.08e-01),vec2(3.65e-01,1.17e-01),vec2(3.73e-01,1.17e-01),vec2(3.98e-01,1.00e-01),vec2(3.90e-01,1.00e-01),vec2(3.90e-01,1.08e-01),vec2(3.90e-01,1.08e-01),vec2(3.90e-01,1.17e-01),vec2(3.98e-01,1.17e-01),vec2(4.41e-01,1.00e-01),vec2(4.49e-01,1.00e-01),vec2(4.49e-01,1.08e-01),vec2(4.49e-01,1.08e-01),vec2(4.49e-01,1.17e-01),vec2(4.41e-01,1.17e-01),vec2(4.41e-01,1.17e-01),vec2(4.32e-01,1.17e-01),vec2(4.32e-01,1.25e-01),vec2(4.32e-01,1.25e-01),vec2(4.32e-01,1.33e-01),vec2(4.41e-01,1.33e-01),vec2(4.66e-01,1.17e-01),vec2(4.58e-01,1.17e-01),vec2(4.58e-01,1.08e-01),vec2(4.58e-01,1.08e-01),vec2(4.58e-01,1.00e-01),vec2(4.66e-01,1.00e-01),vec2(5.00e-01,1.08e-01),vec2(5.00e-01,1.17e-01),vec2(4.92e-01,1.17e-01),vec2(5.17e-01,1.17e-01),vec2(5.25e-01,1.17e-01),vec2(5.25e-01,1.08e-01),vec2(5.09e-01,1.08e-01),vec2(5.09e-01,1.17e-01),vec2(5.17e-01,1.17e-01),vec2(5.42e-01,1.17e-01),vec2(5.34e-01,1.17e-01),vec2(5.34e-01,1.08e-01),vec2(5.34e-01,1.08e-01),vec2(5.34e-01,1.00e-01),vec2(5.42e-01,1.00e-01),vec2(5.42e-01,1.00e-01),vec2(5.51e-01,1.00e-01),vec2(5.51e-01,1.08e-01),vec2(5.42e-01,1.17e-01),vec2(5.51e-01,1.17e-01),vec2(5.51e-01,1.08e-01),vec2(5.68e-01,1.00e-01),vec2(5.76e-01,1.00e-01),vec2(5.76e-01,1.08e-01),vec2(5.76e-01,1.08e-01),vec2(5.76e-01,1.17e-01),vec2(5.68e-01,1.17e-01),vec2(5.93e-01,1.00e-01),vec2(6.02e-01,1.00e-01),vec2(6.02e-01,1.08e-01),vec2(6.02e-01,1.08e-01),vec2(6.02e-01,1.17e-01),vec2(5.93e-01,1.17e-01),vec2(6.19e-01,1.00e-01),vec2(6.35e-01,1.04e-01),vec2(6.19e-01,1.08e-01),vec2(6.19e-01,1.08e-01),vec2(6.02e-01,1.13e-01),vec2(6.19e-01,1.17e-01),vec2(6.52e-01,1.08e-01),vec2(6.52e-01,1.17e-01),vec2(6.44e-01,1.17e-01),vec2(6.44e-01,1.17e-01),vec2(6.36e-01,1.17e-01),vec2(6.36e-01,1.08e-01),vec2(6.36e-01,1.08e-01),vec2(6.36e-01,1.00e-01),vec2(6.44e-01,1.00e-01),vec2(6.44e-01,1.00e-01),vec2(6.52e-01,1.00e-01),vec2(6.52e-01,1.08e-01),vec2(6.52e-01,9.17e-02),vec2(6.52e-01,8.33e-02),vec2(6.44e-01,8.33e-02),vec2(6.61e-01,1.08e-01),vec2(6.61e-01,1.00e-01),vec2(6.70e-01,1.00e-01),vec2(6.78e-01,1.08e-01),vec2(6.78e-01,1.17e-01),vec2(6.87e-01,1.17e-01),vec2(6.95e-01,1.08e-01),vec2(6.95e-01,1.00e-01),vec2(7.04e-01,1.00e-01),vec2(7.21e-01,1.00e-01),vec2(7.37e-01,1.04e-01),vec2(7.21e-01,1.08e-01),vec2(7.21e-01,1.08e-01),vec2(7.04e-01,1.13e-01),vec2(7.21e-01,1.17e-01));
		for(int i=0; i<57;++i) d=min(d,dsg(lin[2*i], lin[2*i+1], uv));
		for(int i=0; i<56; ++i) d=min(d,dsp(quad[3*i], quad[3*i+1], quad[3*i+2], uv));
        col = mix(col, c.yyy, B(80.)*smoothstep(.005, .002, d ));
    }
    
    fragColor = vec4(col,1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}
