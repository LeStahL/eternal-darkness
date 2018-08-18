/* File generated with Shader Minifier 1.1.5
 * http://www.ctrl-alt-test.fr
 */
#ifndef GFX_H_
# define GFX_H_

const char *gfx_frag =
 "#version 130\n"
 "uniform float iTime;"
 "uniform vec2 iResolution;"
 "const float pi=acos(-1.);"
 "const vec3 c=vec3(1.,0.,-1.);"
 "vec2 vi;"
 "float t;\n"
 "#define T.4286\n"
 "float rand(vec2 v)"
 "{"
   "return-1.+2.*fract(sin(dot(v.xy,vec2(12.9898,78.233)))*43758.5);"
 "}"
 "float rand3(vec3 v)"
 "{"
   "return.33333*(rand(v.xy)+rand(v.yz)+rand(v.zx));"
 "}"
 "float dist(vec2 x,vec2 y,vec2 s,vec2 v,float i)"
 "{"
   "return i=clamp(i,0.,1.),length(v-pow(1.-i,2.)*x-2.*(1.-i)*i*y-i*i*s);"
 "}"
 "float dsp(vec2 v,vec2 m,vec2 x,vec2 y)"
 "{"
   "vec2 i=y-v,d=x-2.*m+v,f=m-v;"
   "vec3 s=vec3(3.*dot(f,d),2.*dot(f,f)-dot(i,d),-dot(i,f))/dot(d,d);"
   "float r=s.x/3.,p=s.y-r*s.x,a=-r*(r*r+p)+s.z,z=a*a/4.+p*p*p/27.;"
   "if(z>0.)"
     "{"
       "vec2 l=-.5*a*c.xx+sqrt(z)*c.xz,n=sign(l)*pow(abs(l),c.xx/3.);"
       "return dist(v,m,x,y,n.x+n.y-r);"
     "}"
   "float e=sqrt(-4./3.*p),l=acos(-.5*a*sqrt(-27./p/p/p))/3.;"
   "vec3 n=c.zxz*e*cos(l*c.xxx+c*pi/3.)-r;"
   "return min(dist(v,m,x,y,n.x),min(dist(v,m,x,y,n.y),dist(v,m,x,y,n.z)));"
 "}"
 "float dist(vec3 x,vec3 y,vec3 v,vec3 i,float m)"
 "{"
   "return m=clamp(m,0.,1.),length(i-pow(1.-m,2.)*x-2.*(1.-m)*m*y-m*m*v);"
 "}"
 "float worm(vec3 v,vec3 m,vec3 x,vec3 y)"
 "{"
   "vec3 i=y-v,d=x-2.*m+v,f=m-v,s=vec3(3.*dot(f,d),2.*dot(f,f)-dot(i,d),-dot(i,f))/dot(d,d);"
   "float r=s.x/3.,p=s.y-r*s.x,a=-r*(r*r+p)+s.z,l=a*a/4.+p*p*p/27.;"
   "if(l>0.)"
     "{"
       "vec2 z=-.5*a*c.xx+sqrt(l)*c.xz,n=sign(z)*pow(abs(z),c.xx/3.);"
       "return dist(v,m,x,y,n.x+n.y-r)-(.03+.01*mod(t,T)/T-.005*cos(1.7*pi*(n.x+n.y-r))+.01*abs(sin(2.*pi*10.*(n.x+n.y-r))));"
     "}"
   "float e=sqrt(-4./3.*p),n=acos(-.5*a*sqrt(-27./p/p/p))/3.;"
   "vec3 B=c.zxz*e*cos(n*c.xxx+c*pi/3.)-r;"
   "return min(dist(v,m,x,y,B.x)-(.03+.01*mod(iTime,T)/T-.005*cos(1.7*pi*B.x)+.01*abs(sin(2.*pi*10.*B.x))),min(dist(v,m,x,y,B.y)-(.03+.01*mod(iTime,T)/T-.005*cos(1.7*pi*B.y)+.01*abs(sin(2.*pi*10.*B.y))),dist(v,m,x,y,B.z)-(.03+.01*mod(iTime,T)/T-.005*cos(1.7*pi*B.z)+.01*abs(sin(2.*pi*10.*B.z)))));"
 "}"
 "float dsg(vec2 v,vec2 m,vec2 i)"
 "{"
   "vec2 x=m-v;"
   "float y=clamp(dot(i-v,x)/dot(x,x),0.,1.);"
   "return length(i-mix(v,m,y));"
 "}"
 "float smoothstep_noise1d(float v)"
 "{"
   "float i=floor(v);"
   "v=fract(v);"
   "float x=rand(i*c.xx),y=rand((i+1.)*c.xx);"
   "return mix(x,y,smoothstep(0.,1.,v));"
 "}"
 "float smoothstep_noise2d(vec2 v)"
 "{"
   "vec2 i=floor(v);"
   "v=fract(v);"
   "float x=rand(i),y=rand(i+c.yx),f=rand(i+c.xy),m=rand(i+c.xx);"
   "return mix(mix(x,y,smoothstep(0.,1.,v.y)),mix(f,m,smoothstep(0.,1.,v.y)),smoothstep(0.,1.,v.x));"
 "}"
 "float mfsmoothstep_noise2d(vec2 x,float v,float y,float i)"
 "{"
   "float m=0.,s=1.2;"
   "for(float f=v;f<y;f=f*2.)"
     "m=s*smoothstep_noise2d(f*x)+m,s=s*i;"
   "return m;"
 "}"
 "float cs(vec2 v,float m,float x,float y,float i)"
 "{"
   "float d=length(v),f=acos(v.x/d)*step(0.,v.y)-acos(v.x/d)*step(v.y,0.);"
   "f=clamp(f,y,i);"
   "vec2 s=m*vec2(cos(f),sin(f));"
   "return length(v-s)-x;"
 "}"
 "float b(vec2 v,vec2 i,vec2 m,float x)"
 "{"
   "vec2 y=m-i;"
   "return length(v-mix(i,m,clamp(dot(v-i,y)/dot(y,y),0.,1.)))-x;"
 "}"
 "void wave(in int v,out float x,out float y,out vec2 m,out float i,out float f)"
 "{"
   "x=55.*rand(float(-v)*c.xx),y=.00012*pow(.9,float(v-1))*rand(float(v)*c.xx),m=c.yx*rand(float(v)*c.xx)+.3*c.xy*rand(float(v)*c.xy),i=pow(.5+.1*rand(float(v+1)*c.xx),float(v))*(1.+.2*rand(float(v)*c.xx))*700.,f=12.*mod(float(v*3),2.5)*rand(vec2(float(v+4)));"
 "}"
 "float gerst(vec2 v,int x)"
 "{"
   "float i=0.,y,m,f,s;"
   "vec2 d;"
   "for(int r=0;r<x;++r)"
     "{"
       "wave(r,y,s,d,m,f);"
       "float p=dot(d,v);"
       "i+=s*sin(m*p+f*t);"
     "}"
   "return i;"
 "}"
 "vec3 vor(vec2 v)"
 "{"
   "v+=12.;"
   "vec2 m=floor(v);"
   "float i=1.;"
   "vec2 x=c.yy,f;"
   "float s=100.,y;"
   "vec3 d=c.xyy*100.;"
   "for(float r=-2.;r<=2.;r+=1.)"
     "for(float e=-2.;e<=2.;e+=1.)"
       "f=m+vec2(r,e),f+=rand(f),y=length(v-f),d=mix(d,vec3(y,f),step(y,d.x));"
   "for(float r=-2.;r<=2.;r+=1.)"
     "for(float e=-2.;e<=2.;e+=1.)"
       "{"
         "f=m+vec2(r,e);"
         "f+=rand(f);"
         "vec2 p=f-d.yz;"
         "y=dot(.5*p-v+d.yz,p)/length(p);"
         "i=min(i,y);"
       "}"
   "return vec3(i,x);"
 "}"
 "vec3 rot(vec3 x,vec3 v)"
 "{"
   "vec3 i=cos(v),s=sin(v);"
   "return mat3(i.x,0.,s.x,0.,1.,0.,-s.x,0.,i.x)*mat3(1.,0.,0.,0.,i.y,-s.y,0.,s.y,i.y)*mat3(i.z,-s.z,0.,s.z,i.z,0.,0.,0.,1.)*x;"
 "}"
 "float dstar(vec2 v,float x,vec2 i)"
 "{"
   "float f=pi/x,y=acos(v.x/length(v)),d=mod(y,f),s=mod(round((d-y)/f),2.);"
   "v=length(v)*vec2(cos(d),sin(d));"
   "vec2 m=mix(i,i.yx,s),r=m.x*c.xy,n=m.y*vec2(cos(f),sin(f))-r;"
   "n=n.yx*c.zx;"
   "return dot(v-r,n)/length(n);"
 "}"
 "float dpoly_min(vec2 v,float x,float y)"
 "{"
   "float i=2.*pi/x,f=mod(acos(v.x/length(v)),i)-.5*i;"
   "return y-length(v)*cos(f)/cos(.5*i);"
 "}"
 "vec2 z10presents(vec2 v)"
 "{"
   "vec2 i=vec2(b(v,vec2(-.5,-.4),vec2(iResolution.x/iResolution.y,-.4),.05),1.),y=vec2(b(v,vec2(-.5,-.4),vec2(iResolution.x/iResolution.y,-.4),.02),2.);"
   "return mix(i,y,step(1.,y.x));"
 "}"
 "vec2 scene1(vec3 v)"
 "{"
   "return v+=c.xyy*.05*t,vec2(v.z+.5-.5*mfsmoothstep_noise2d(v.xy,.9,1.,.1)+.05*mfsmoothstep_noise2d(v.xy,120.,180.,1.7),1.);"
 "}"
 "vec2 scene2(vec3 v)"
 "{"
   "v+=-13.*c.yxy+c.xxy*25.+c.yxy*.05*t;"
   "float x=-.5*mfsmoothstep_noise2d(v.xy,.91,2.,.1);"
   "if(v.z>-.19)"
     "x+=.05*mfsmoothstep_noise2d(v.xy,120.,180.,1.7);"
   "vec2 i=vec2(v.z+x,1.),y=vec2(v.z+.21-gerst(v.xy,22),2.);"
   "return mix(i,y,step(y.x,i.x));"
 "}"
 "vec2 scene3(vec3 v)"
 "{"
   "v+=c.yxy*.5*t+.5*c.yyx;"
   "vec3 i=vec3(mod(v.x,2.)-1.,mod(v.y,2.)-1.,v.z),f=v-i;"
   "vec2 y=vec2(length(i.xy-.5*vec2(rand(f.xy),rand(f.xy*1.1)))-.1,1.),x=vec2(v.z-.2*mfsmoothstep_noise2d(v.xy,.9,1.,.1),1.);"
   "float d=-length(max(abs(i)-vec3(1.,2.,50.),0.));"
   "d=abs(d)+.4;"
   "y.x=min(y.x,d);"
   "return mix(x,y,step(y.x,x.x));"
 "}"
 "vec2 scene4(vec3 v)"
 "{"
   "vec3 i=v;"
   "v+=t*c.yxy;"
   "vec3 m=v;"
   "v=mod(v,1.)-.5;"
   "vec3 x=m-v;"
   "float s=length(x);"
   "vec3 y=.25*vec3(rand(s*c.xx),rand(3.*s*c.xx),rand(6.*s*c.xx));"
   "vec2 f=vec2(worm(c.yyy,vec3(.1*sin(5.*t+10.*rand(length(m-v)*c.xx)),0.,.1),.2*c.yyx,v-y),3.),d=vec2(length(v-y)-(.036+.01*mod(t,T)/T),4.),r=mix(f,d,step(d.x,f.x)),n=vec2(length(v-.2*c.yyx-y)-(.036+.01*mod(t,T)/T),4.);"
   "r=mix(r,n,step(n.x,r.x));"
   "float z=3.*acos(min(m.z,m.x)/length(m.xz));"
   "vec3 a=.5*vor(vec2(m.y+z,z));"
   "float l=a.x;"
   "vi=a.yz;"
   "f=vec2(abs(length(m.xz)-3.)-l,5.);"
   "r=mix(r,f,step(f.x,r.x));"
   "float e=-length(max(abs(v)-.5,0.));"
   "e=abs(e)+.1;"
   "r.x=min(r.x,e);"
   "return r;"
 "}"
 "vec2 scene5(vec3 v)"
 "{"
   "vec3 i=v;"
   "v+=t*c.yxy;"
   "vec3 m=v;"
   "v=vec3(mod(v.xy,1.)-.5,v.z);"
   "vec3 f=m-v;"
   "vi=f.xy;"
   "float x=2.*pi*rand(f.xy)*t;"
   "v.xy=mat2(cos(x),sin(x),-sin(x),cos(x))*v.xy;"
   "float y=dstar(v.xy,3.+floor(7.*abs(rand(f.xy+6.))),vec2(.2+.05*rand(vi),.3+.05*rand(vi+3.)));"
   "vec2 r=abs(vec2(min(y,0.),v.z+1.3))-.4*c.yx-(.5+2.*mod(t,T))*.5*abs(rand(f.xy+t-mod(t,T)))*c.yx,d=vec2(min(max(r.x,r.y),0.)+length(max(r,0.)),6.);"
   "r=abs(vec2(min(y,0.),v.z-1.3))-.4*c.yx-(.5+2.*mod(t,T))*.5*rand(f.xy+t-mod(t,T))*c.yx;"
   "vec2 s=vec2(min(max(r.x,r.y),0.)+length(max(r,0.)),6.);"
   "d=mix(d,s,step(s.x,d.x));"
   "float e=-length(max(abs(v+50.*c.yyx)-vec3(.5,.5,100.),0.));"
   "e=abs(e)+.1;"
   "d.x=min(d.x,e);"
   "return d;"
 "}"
 "vec2 scene(vec3 v)"
 "{"
   "if(t<20.5728)"
     "return scene1(v);"
   "else"
     " if(t<41.1456)"
       "return scene2(v);"
     "else"
       " if(t<61.7184)"
         "return scene3(v);"
       "else"
         " if(t<78.8624)"
           "return scene4(v);"
         "else"
           " if(t<7000.)"
             "return scene5(v);"
 "}"
 "const float dx=.0001;"
 "vec3 normal(vec3 v)"
 "{"
   "float x=scene(v).x;"
   "return normalize(vec3(scene(v+dx*c.xyy).x-x,scene(v+dx*c.yxy).x-x,scene(v+dx*c.yyx).x-x));"
 "}"
 "vec3 bg(vec2 v)"
 "{"
   "if(t<20.)"
     "{"
       "v+=.01*t*c.yx;"
       "vec2 i=mod(v,.05)-.025,f=v-i;"
       "float y=.17*exp(-.2*t);"
       "return.15*c.xxx+mix(.4*c.xyy,.8*c.xxy,1.-2.*v.y-length(v))*exp(-.2*t)+c.xxx*smoothstep(y+.01,y-.01,length(v))+2.*tanh(.2*t)*smoothstep(.003,-.003,length(i-.02*vec2(rand(f),rand(.1*f+.1*c.xx)))-.002*rand(.1*f))*(.2*vec3(rand(.1*f),rand(.1*f+.1*c.xx),rand(f+.3*c.xx))+.8)+(.4+.2*mfsmoothstep_noise2d(v+cos(.5*t)*sin(.5*v.x)*sin(.6*v.y),2.5,100.,.55))*exp(-.05*t);"
     "}"
   "else"
     " return c.xxx;"
 "}"
 "float B(float v)"
 "{"
   "return smoothstep(v-2.,v-1.,iTime)*(1.-smoothstep(v+2.,v+3.,iTime));"
 "}"
 "void fore(out vec4 v,in vec2 i,float m)"
 "{"
   "t=m;"
   "vec2 f;"
   "vec3 x,y=.5+.5*mfsmoothstep_noise2d(c.yy,.9,1.,.1)-c.yxy+.5*c.yyx,s=.25*c.yyx,d=c.xyy,r=c.yyx,n=s+d*i.x+r*i.y,p=normalize(n-y);"
   "if(t>61.92)"
     "y=-c.yxy,s=c.yyy;"
   "float e=0.,B=0.,l=0.,o=15.;"
   "int a=100;"
   "if(m>20.64)"
     "o=55.;"
   "if(m>41.28)"
     "o=50.,a=200;"
   "else"
     " if(m>61.92)"
       "o=200.;"
   "for(int z=0;z<a;++z)"
     "{"
       "x=y+e*p;"
       "f=scene(x);"
       "if(f.x<.005)"
         "break;"
       "if(e>o||z==a-1)"
         "{"
           "v=vec4(bg(i),1.);"
           "return;"
         "}"
       "e+=f.x;"
     "}"
   "vec3 z=normal(x),b=c.yxy,g=x+c.yyx,u=normalize(x-y),w=normalize(reflect(-b,z)),h=normalize(reflect(-g,z)),k;"
   "if(f.y==1.)"
     "k=.01*c.xyy*dot(b,z)+.2*c.xxx*pow(abs(dot(w,u)),4.)+(.05*c.xyy+.05*c.yxy)*pow(abs(dot(w,u)),2.),k*=exp(-.2*iTime);"
   "else"
     " if(f.y==2.)"
       "k=.4*c.xxx*dot(g,z)+.2*c.xxx*pow(abs(dot(h,u)),2.);"
     "else"
       " if(f.y==3.)"
         "k=.6-(.8*c.yyx+.3*c.xyy)*dot(b,z)*mod(t,T)/T+.8*c.xxx*dot(b,z)+.2*c.xxx*pow(abs(dot(w,u)),4.)+(.05*c.xyy+.05*c.yxy)*pow(abs(dot(w,u)),2.);"
       "else"
         " if(f.y==4.)"
           "k=1.+.8*c.xxx*dot(b,z)+.2*c.xxx*pow(abs(dot(w,u)),4.)+(.05*c.yyx+.05*c.yyx)*pow(abs(dot(w,u)),2.);"
         "else"
           " if(f.y==5.)"
             "b=normalize(x),w=normalize(reflect(-b,z)),u=normalize(x-y),k=.6+.3*vec3(.75+.1*rand(vi),.75+.1*rand(2.*vi),.6+.4*rand(5.*vi))*dot(b,z)*mod(t,T)/T+.8*c.xxx*dot(b,z)+.05*c.xxx*pow(abs(dot(w,u)),2.);"
           "else"
             " if(f.y==6.)"
               "{"
                 "b=normalize(x-c.yxy);"
                 "w=normalize(reflect(-b,z));"
                 "u=normalize(x-c.yxy-y);"
                 "vec3 R=.5+.5*cos(x+vec3(0.,2.,4.)+.5*abs(rand(vi))+m),q=.5+.5*cos(x+vec3(0.,2.,4.)+.6*abs(rand(vi+3.))+m),I=.5+.5*cos(x+vec3(0.,2.,4.)+.7*abs(rand(vi+7.))+m);"
                 "k=R+q*dot(b,z)+.1*I*pow(abs(dot(w,u)),2.);"
               "}"
             "else"
               " k=c.yyy;"
   "if(m<20.5728)"
     "k=mix(k,.15*c.xxx,cosh(-2.*x.z)*tanh(.109*x.y));"
   "else"
     " if(m<41.1456)"
       "k=mix(k,c.xxx,tanh(.1*x.y));"
     "else"
       " if(m<7000.)"
         "k=mix(k,c.xxx,tanh(.08*x.y));"
   "v=vec4(k,1.);"
 "}"
 "void mainImage(out vec4 v,in vec2 x)"
 "{"
   "vec2 y=x/iResolution.yy-.5,f;"
   "vec4 m,i,d;"
   "fore(m,y,iTime);"
   "if(iTime<40.)"
     "fore(i,y,iTime+.01),v=.333*(m+i);"
   "else"
     " v=m;"
   "vec3 r=v.xyz;"
   "if(iTime<180.)"
     "{"
       "f=z10presents(y);"
       "float s=step(f.x,0.);"
       "if(f.y==1.)"
         "r=mix(r,mix(r,c.xxx,.1),s);"
     "}"
   "vec2 s=y+(.01*iTime-.12)*c.yx;"
   "float p=atan(abs(s.y/s.x)),z=length(s);"
   "r+=1.4*(.5+.25*sin(2.3*pi*p)+.25*sin(4.*pi*p)+.25*sin(1.5*pi*p))*exp(-3.*z)*c.xxx*exp(-.9*t)*(1.-smoothstep(1.4,1.6,iTime));"
   "float e=1.;"
   "if(iTime<8.)"
     "{"
       "const vec2 n[84]=vec2[84](vec2(-.388,-.415),vec2(-.388,-.368),vec2(-.4,-.368),vec2(-.377,-.368),vec2(-.353,-.415),vec2(-.341,-.415),vec2(-.364,-.403),vec2(-.341,-.403),vec2(-.305,-.415),vec2(-.305,-.392),vec2(-.293,-.415),vec2(-.293,-.392),vec2(-.27,-.415),vec2(-.27,-.403),vec2(-.247,-.415),vec2(-.247,-.403),vec2(-.234,-.415),vec2(-.211,-.415),vec2(-.199,-.38),vec2(-.187,-.368),vec2(-.187,-.368),vec2(-.187,-.415),vec2(-.175,-.403),vec2(-.175,-.38),vec2(-.151,-.403),vec2(-.151,-.38),vec2(-.115,-.438),vec2(-.115,-.392),vec2(-.115,-.392),vec2(-.104,-.392),vec2(-.115,-.415),vec2(-.104,-.415),vec2(-.0797,-.415),vec2(-.0797,-.392),vec2(-.0202,-.392),vec2(-.0202,-.403),vec2(.00308,-.392),vec2(.00308,-.415),vec2(.027,-.392),vec2(.0387,-.392),vec2(.027,-.415),vec2(.0387,-.415),vec2(.0387,-.415),vec2(.0387,-.368),vec2(.0509,-.403),vec2(.0509,-.368),vec2(.0748,-.392),vec2(.0748,-.403),vec2(.0982,-.392),vec2(.0982,-.427),vec2(.0865,-.438),vec2(.0748,-.438),vec2(.134,-.438),vec2(.134,-.392),vec2(.134,-.392),vec2(.146,-.392),vec2(.134,-.415),vec2(.146,-.415),vec2(.17,-.415),vec2(.17,-.392),vec2(.206,-.415),vec2(.217,-.415),vec2(.194,-.403),vec2(.217,-.403),vec2(.229,-.415),vec2(.241,-.415),vec2(.241,-.392),vec2(.253,-.392),vec2(.277,-.415),vec2(.288,-.415),vec2(.265,-.403),vec2(.288,-.403),vec2(.301,-.415),vec2(.301,-.392),vec2(.324,-.403),vec2(.324,-.415),vec2(.348,-.403),vec2(.348,-.368),vec2(.336,-.392),vec2(.36,-.392),vec2(.372,-.415),vec2(.383,-.415),vec2(.383,-.392),vec2(.395,-.392)),a[141]=vec2[141](vec2(-.341,-.403),vec2(-.341,-.392),vec2(-.353,-.392),vec2(-.353,-.392),vec2(-.364,-.392),vec2(-.364,-.403),vec2(-.364,-.403),vec2(-.364,-.415),vec2(-.353,-.415),vec2(-.317,-.392),vec2(-.329,-.392),vec2(-.329,-.403),vec2(-.329,-.403),vec2(-.329,-.415),vec2(-.317,-.415),vec2(-.317,-.415),vec2(-.305,-.415),vec2(-.305,-.403),vec2(-.317,-.392),vec2(-.305,-.392),vec2(-.305,-.403),vec2(-.293,-.403),vec2(-.293,-.392),vec2(-.282,-.392),vec2(-.282,-.392),vec2(-.27,-.392),vec2(-.27,-.403),vec2(-.27,-.403),vec2(-.27,-.392),vec2(-.258,-.392),vec2(-.258,-.392),vec2(-.247,-.392),vec2(-.247,-.403),vec2(-.234,-.368),vec2(-.199,-.357),vec2(-.234,-.415),vec2(-.175,-.403),vec2(-.175,-.415),vec2(-.163,-.415),vec2(-.163,-.415),vec2(-.151,-.415),vec2(-.151,-.403),vec2(-.151,-.38),vec2(-.151,-.368),vec2(-.163,-.368),vec2(-.163,-.368),vec2(-.175,-.368),vec2(-.175,-.38),vec2(-.104,-.415),vec2(-.092,-.415),vec2(-.092,-.403),vec2(-.092,-.403),vec2(-.092,-.392),vec2(-.104,-.392),vec2(-.0797,-.403),vec2(-.0797,-.392),vec2(-.0681,-.392),vec2(-.0442,-.415),vec2(-.0558,-.415),vec2(-.0558,-.403),vec2(-.0558,-.403),vec2(-.0558,-.392),vec2(-.0442,-.392),vec2(-.0442,-.392),vec2(-.0325,-.392),vec2(-.0325,-.403),vec2(-.0325,-.403),vec2(-.0325,-.415),vec2(-.0442,-.415),vec2(-.0202,-.403),vec2(-.0202,-.415),vec2(-.00858,-.415),vec2(-.00858,-.415),vec2(.00308,-.415),vec2(.00308,-.403),vec2(.027,-.415),vec2(.0153,-.415),vec2(.0153,-.403),vec2(.0153,-.403),vec2(.0153,-.392),vec2(.027,-.392),vec2(.0509,-.403),vec2(.0509,-.415),vec2(.0626,-.415),vec2(.0748,-.403),vec2(.0748,-.415),vec2(.0865,-.415),vec2(.0865,-.415),vec2(.0982,-.415),vec2(.0982,-.403),vec2(.0982,-.427),vec2(.0982,-.438),vec2(.0865,-.438),vec2(.146,-.415),vec2(.158,-.415),vec2(.158,-.403),vec2(.158,-.403),vec2(.158,-.392),vec2(.146,-.392),vec2(.17,-.403),vec2(.17,-.392),vec2(.182,-.392),vec2(.217,-.403),vec2(.217,-.392),vec2(.206,-.392),vec2(.206,-.392),vec2(.194,-.392),vec2(.194,-.403),vec2(.194,-.403),vec2(.194,-.415),vec2(.206,-.415),vec2(.241,-.415),vec2(.264,-.409),vec2(.241,-.403),vec2(.241,-.403),vec2(.218,-.398),vec2(.241,-.392),vec2(.288,-.403),vec2(.288,-.392),vec2(.277,-.392),vec2(.277,-.392),vec2(.265,-.392),vec2(.265,-.403),vec2(.265,-.403),vec2(.265,-.415),vec2(.277,-.415),vec2(.312,-.392),vec2(.324,-.392),vec2(.324,-.403),vec2(.301,-.403),vec2(.301,-.392),vec2(.312,-.392),vec2(.348,-.403),vec2(.348,-.415),vec2(.36,-.415),vec2(.383,-.415),vec2(.407,-.409),vec2(.383,-.403),vec2(.383,-.403),vec2(.36,-.398),vec2(.383,-.392));"
       "for(int l=0;l<42;++l)"
         "e=min(e,dsg(n[2*l],n[2*l+1],y));"
       "for(int l=0;l<47;++l)"
         "e=min(e,dsp(a[3*l],a[3*l+1],a[3*l+2],y));"
       "r=mix(r,.8*c.xxx,B(5.)*smoothstep(.005,.002,e));"
     "}"
   "else"
     " if(iTime<13.)"
       "{"
         "const vec2 a[58]=vec2[58](vec2(-.399,-.415),vec2(-.399,-.368),vec2(-.399,-.415),vec2(-.376,-.415),vec2(-.399,-.392),vec2(-.376,-.392),vec2(-.399,-.368),vec2(-.376,-.368),vec2(-.352,-.403),vec2(-.352,-.368),vec2(-.364,-.392),vec2(-.341,-.392),vec2(-.317,-.415),vec2(-.305,-.415),vec2(-.328,-.403),vec2(-.305,-.403),vec2(-.293,-.415),vec2(-.293,-.392),vec2(-.269,-.415),vec2(-.269,-.392),vec2(-.245,-.403),vec2(-.245,-.415),vec2(-.21,-.415),vec2(-.21,-.392),vec2(-.198,-.403),vec2(-.198,-.368),vec2(-.15,-.415),vec2(-.15,-.368),vec2(-.15,-.368),vec2(-.138,-.368),vec2(-.15,-.415),vec2(-.138,-.415),vec2(-.126,-.403),vec2(-.126,-.38),vec2(-.0908,-.415),vec2(-.0908,-.392),vec2(-.0786,-.415),vec2(-.0786,-.392),vec2(-.0547,-.415),vec2(-.0547,-.368),vec2(-.0547,-.403),vec2(-.043,-.403),vec2(-.0191,-.415),vec2(-.0191,-.392),vec2(.00425,-.403),vec2(.00425,-.415),vec2(.0282,-.415),vec2(.0398,-.415),vec2(.0165,-.403),vec2(.0398,-.403),vec2(.0521,-.415),vec2(.0638,-.415),vec2(.0638,-.392),vec2(.0754,-.392),vec2(.0877,-.415),vec2(.0993,-.415),vec2(.0993,-.392),vec2(.111,-.392)),n[90]=vec2[90](vec2(-.352,-.403),vec2(-.352,-.415),vec2(-.341,-.415),vec2(-.305,-.403),vec2(-.305,-.392),vec2(-.317,-.392),vec2(-.317,-.392),vec2(-.328,-.392),vec2(-.328,-.403),vec2(-.328,-.403),vec2(-.328,-.415),vec2(-.317,-.415),vec2(-.293,-.403),vec2(-.293,-.392),vec2(-.281,-.392),vec2(-.257,-.392),vec2(-.245,-.392),vec2(-.245,-.403),vec2(-.269,-.403),vec2(-.269,-.392),vec2(-.257,-.392),vec2(-.222,-.392),vec2(-.233,-.392),vec2(-.233,-.403),vec2(-.233,-.403),vec2(-.233,-.415),vec2(-.222,-.415),vec2(-.222,-.415),vec2(-.21,-.415),vec2(-.21,-.403),vec2(-.222,-.392),vec2(-.21,-.392),vec2(-.21,-.403),vec2(-.198,-.403),vec2(-.198,-.415),vec2(-.186,-.415),vec2(-.138,-.368),vec2(-.126,-.368),vec2(-.126,-.38),vec2(-.138,-.415),vec2(-.126,-.415),vec2(-.126,-.403),vec2(-.102,-.392),vec2(-.114,-.392),vec2(-.114,-.403),vec2(-.114,-.403),vec2(-.114,-.415),vec2(-.102,-.415),vec2(-.102,-.415),vec2(-.0908,-.415),vec2(-.0908,-.403),vec2(-.102,-.392),vec2(-.0908,-.392),vec2(-.0908,-.403),vec2(-.0786,-.403),vec2(-.0786,-.392),vec2(-.0669,-.392),vec2(-.043,-.403),vec2(-.0313,-.403),vec2(-.0313,-.392),vec2(-.043,-.403),vec2(-.0313,-.403),vec2(-.0313,-.415),vec2(-.00742,-.392),vec2(.00425,-.392),vec2(.00425,-.403),vec2(-.0191,-.403),vec2(-.0191,-.392),vec2(-.00742,-.392),vec2(.0398,-.403),vec2(.0398,-.392),vec2(.0282,-.392),vec2(.0282,-.392),vec2(.0165,-.392),vec2(.0165,-.403),vec2(.0165,-.403),vec2(.0165,-.415),vec2(.0282,-.415),vec2(.0638,-.415),vec2(.0871,-.409),vec2(.0638,-.403),vec2(.0638,-.403),vec2(.0404,-.398),vec2(.0638,-.392),vec2(.0993,-.415),vec2(.123,-.409),vec2(.0993,-.403),vec2(.0993,-.403),vec2(.076,-.398),vec2(.0993,-.392));"
         "for(int l=0;l<29;++l)"
           "e=min(e,dsg(a[2*l],a[2*l+1],y));"
         "for(int l=0;l<30;++l)"
           "e=min(e,dsp(n[3*l],n[3*l+1],n[3*l+2],y));"
         "r=mix(r,.8*c.xxx,B(10.)*smoothstep(.005,.002,e));"
       "}"
     "else"
       " if(iTime<18.)"
         "{"
           "const vec2 n[18]=vec2[18](vec2(-.4,-.368),vec2(-.4,-.403),vec2(-.377,-.368),vec2(-.377,-.403),vec2(-.353,-.368),vec2(-.353,-.403),vec2(-.329,-.415),vec2(-.318,-.415),vec2(-.341,-.403),vec2(-.318,-.403),vec2(-.258,-.415),vec2(-.258,-.392),vec2(-.246,-.415),vec2(-.246,-.392),vec2(-.21,-.415),vec2(-.199,-.415),vec2(-.222,-.403),vec2(-.199,-.403)),a[45]=vec2[45](vec2(-.4,-.403),vec2(-.4,-.415),vec2(-.388,-.415),vec2(-.388,-.415),vec2(-.377,-.415),vec2(-.377,-.403),vec2(-.377,-.403),vec2(-.377,-.415),vec2(-.365,-.415),vec2(-.365,-.415),vec2(-.353,-.415),vec2(-.353,-.403),vec2(-.318,-.403),vec2(-.318,-.392),vec2(-.329,-.392),vec2(-.329,-.392),vec2(-.341,-.392),vec2(-.341,-.403),vec2(-.341,-.403),vec2(-.341,-.415),vec2(-.329,-.415),vec2(-.27,-.392),vec2(-.282,-.392),vec2(-.282,-.403),vec2(-.282,-.403),vec2(-.282,-.415),vec2(-.27,-.415),vec2(-.27,-.415),vec2(-.258,-.415),vec2(-.258,-.403),vec2(-.27,-.392),vec2(-.258,-.392),vec2(-.258,-.403),vec2(-.246,-.403),vec2(-.246,-.392),vec2(-.234,-.392),vec2(-.199,-.403),vec2(-.199,-.392),vec2(-.21,-.392),vec2(-.21,-.392),vec2(-.222,-.392),vec2(-.222,-.403),vec2(-.222,-.403),vec2(-.222,-.415),vec2(-.21,-.415));"
           "for(int l=0;l<9;++l)"
             "e=min(e,dsg(n[2*l],n[2*l+1],y));"
           "for(int l=0;l<15;++l)"
             "e=min(e,dsp(a[3*l],a[3*l+1],a[3*l+2],y));"
           "r=mix(r,.8*c.xxx,B(15.)*smoothstep(.005,.002,e));"
         "}"
       "else"
         " if(iTime<23.)"
           "{"
             "const vec2 a[128]=vec2[128](vec2(-.4,-.403),vec2(-.4,-.38),vec2(-.377,-.403),vec2(-.377,-.38),vec2(-.364,-.415),vec2(-.364,-.368),vec2(-.318,-.415),vec2(-.318,-.368),vec2(-.364,-.368),vec2(-.341,-.392),vec2(-.341,-.392),vec2(-.318,-.368),vec2(-.282,-.403),vec2(-.282,-.403),vec2(-.282,-.38),vec2(-.282,-.38),vec2(-.269,-.403),vec2(-.269,-.403),vec2(-.269,-.38),vec2(-.269,-.38),vec2(-.21,-.415),vec2(-.221,-.415),vec2(-.21,-.368),vec2(-.221,-.368),vec2(-.233,-.403),vec2(-.233,-.38),vec2(-.15,-.392),vec2(-.139,-.392),vec2(-.15,-.415),vec2(-.139,-.415),vec2(-.139,-.415),vec2(-.139,-.368),vec2(-.115,-.415),vec2(-.103,-.415),vec2(-.126,-.403),vec2(-.103,-.403),vec2(-.0908,-.415),vec2(-.0675,-.368),vec2(-.0552,-.415),vec2(-.0436,-.415),vec2(-.0319,-.368),vec2(-.0436,-.368),vec2(-.0197,-.415),vec2(-.0197,-.368),vec2(-.0197,-.368),vec2(.00367,-.368),vec2(-.0197,-.392),vec2(.00367,-.392),vec2(.0159,-.415),vec2(.0393,-.368),vec2(.0159,-.368),vec2(.0393,-.415),vec2(.0988,-.415),vec2(.0988,-.392),vec2(.111,-.415),vec2(.111,-.392),vec2(.134,-.403),vec2(.134,-.415),vec2(.158,-.392),vec2(.17,-.392),vec2(.158,-.415),vec2(.17,-.415),vec2(.17,-.415),vec2(.17,-.368),vec2(.206,-.415),vec2(.206,-.368),vec2(.206,-.368),vec2(.229,-.415),vec2(.229,-.415),vec2(.229,-.368),vec2(.242,-.415),vec2(.242,-.368),vec2(.242,-.368),vec2(.253,-.368),vec2(.242,-.392),vec2(.253,-.392),vec2(.265,-.403),vec2(.265,-.415),vec2(.301,-.368),vec2(.277,-.403),vec2(.277,-.403),vec2(.301,-.403),vec2(.301,-.392),vec2(.301,-.415),vec2(.337,-.403),vec2(.337,-.403),vec2(.337,-.38),vec2(.337,-.38),vec2(.349,-.403),vec2(.349,-.403),vec2(.349,-.38),vec2(.349,-.38),vec2(.409,-.415),vec2(.397,-.415),vec2(.409,-.368),vec2(.397,-.368),vec2(.385,-.403),vec2(.385,-.38),vec2(.468,-.392),vec2(.48,-.392),vec2(.468,-.415),vec2(.48,-.415),vec2(.48,-.415),vec2(.48,-.368),vec2(.504,-.415),vec2(.515,-.415),vec2(.492,-.403),vec2(.515,-.403),vec2(.527,-.415),vec2(.551,-.368),vec2(.575,-.392),vec2(.586,-.392),vec2(.586,-.392),vec2(.586,-.403),vec2(.563,-.403),vec2(.563,-.38),vec2(.575,-.368),vec2(.586,-.368),vec2(.599,-.415),vec2(.599,-.368),vec2(.599,-.368),vec2(.622,-.368),vec2(.599,-.392),vec2(.622,-.392),vec2(.634,-.415),vec2(.658,-.368),vec2(.634,-.368),vec2(.658,-.415)),n[135]=vec2[135](vec2(-.4,-.403),vec2(-.4,-.415),vec2(-.388,-.415),vec2(-.388,-.415),vec2(-.377,-.415),vec2(-.377,-.403),vec2(-.377,-.38),vec2(-.377,-.368),vec2(-.388,-.368),vec2(-.388,-.368),vec2(-.4,-.368),vec2(-.4,-.38),vec2(-.388,-.403),vec2(-.388,-.415),vec2(-.377,-.415),vec2(-.221,-.415),vec2(-.233,-.415),vec2(-.233,-.403),vec2(-.233,-.38),vec2(-.233,-.368),vec2(-.221,-.368),vec2(-.186,-.415),vec2(-.198,-.415),vec2(-.198,-.403),vec2(-.198,-.403),vec2(-.198,-.392),vec2(-.186,-.392),vec2(-.186,-.392),vec2(-.174,-.392),vec2(-.174,-.403),vec2(-.174,-.403),vec2(-.174,-.415),vec2(-.186,-.415),vec2(-.15,-.415),vec2(-.162,-.415),vec2(-.162,-.403),vec2(-.162,-.403),vec2(-.162,-.392),vec2(-.15,-.392),vec2(-.103,-.403),vec2(-.103,-.392),vec2(-.115,-.392),vec2(-.115,-.392),vec2(-.126,-.392),vec2(-.126,-.403),vec2(-.126,-.403),vec2(-.126,-.415),vec2(-.115,-.415),vec2(-.0436,-.415),vec2(-.0319,-.415),vec2(-.0319,-.403),vec2(-.0319,-.403),vec2(-.0319,-.392),vec2(-.0436,-.392),vec2(-.0436,-.392),vec2(-.0552,-.392),vec2(-.0552,-.38),vec2(-.0552,-.38),vec2(-.0552,-.368),vec2(-.0436,-.368),vec2(.0871,-.392),vec2(.0754,-.392),vec2(.0754,-.403),vec2(.0754,-.403),vec2(.0754,-.415),vec2(.0871,-.415),vec2(.0871,-.415),vec2(.0988,-.415),vec2(.0988,-.403),vec2(.0871,-.392),vec2(.0988,-.392),vec2(.0988,-.403),vec2(.123,-.392),vec2(.134,-.392),vec2(.134,-.403),vec2(.111,-.403),vec2(.111,-.392),vec2(.123,-.392),vec2(.158,-.415),vec2(.147,-.415),vec2(.147,-.403),vec2(.147,-.403),vec2(.147,-.392),vec2(.158,-.392),vec2(.253,-.368),vec2(.265,-.368),vec2(.265,-.38),vec2(.265,-.38),vec2(.265,-.392),vec2(.253,-.392),vec2(.253,-.392),vec2(.265,-.392),vec2(.265,-.403),vec2(.397,-.415),vec2(.385,-.415),vec2(.385,-.403),vec2(.385,-.38),vec2(.385,-.368),vec2(.397,-.368),vec2(.432,-.415),vec2(.421,-.415),vec2(.421,-.403),vec2(.421,-.403),vec2(.421,-.392),vec2(.432,-.392),vec2(.432,-.392),vec2(.444,-.392),vec2(.444,-.403),vec2(.444,-.403),vec2(.444,-.415),vec2(.432,-.415),vec2(.468,-.415),vec2(.456,-.415),vec2(.456,-.403),vec2(.456,-.403),vec2(.456,-.392),vec2(.468,-.392),vec2(.515,-.403),vec2(.515,-.392),vec2(.504,-.392),vec2(.504,-.392),vec2(.492,-.392),vec2(.492,-.403),vec2(.492,-.403),vec2(.492,-.415),vec2(.504,-.415),vec2(.563,-.38),vec2(.563,-.368),vec2(.575,-.368),vec2(.575,-.415),vec2(.563,-.415),vec2(.563,-.403),vec2(.575,-.415),vec2(.586,-.415),vec2(.586,-.403));"
             "for(int l=0;l<64;++l)"
               "e=min(e,dsg(a[2*l],a[2*l+1],y));"
             "for(int l=0;l<45;++l)"
               "e=min(e,dsp(n[3*l],n[3*l+1],n[3*l+2],y));"
             "r=mix(r,.8*c.xxx,B(20.)*smoothstep(.005,.002,e));"
           "}"
         "else"
           " if(iTime<28.)"
             "{"
               "const vec2 n[24]=vec2[24](vec2(-.4,-.368),vec2(-.4,-.403),vec2(-.377,-.368),vec2(-.377,-.403),vec2(-.353,-.368),vec2(-.353,-.403),vec2(-.329,-.415),vec2(-.318,-.415),vec2(-.341,-.403),vec2(-.318,-.403),vec2(-.282,-.403),vec2(-.282,-.368),vec2(-.258,-.403),vec2(-.258,-.392),vec2(-.258,-.38),vec2(-.258,-.38),vec2(-.234,-.415),vec2(-.234,-.368),vec2(-.234,-.403),vec2(-.222,-.403),vec2(-.186,-.415),vec2(-.175,-.415),vec2(-.198,-.403),vec2(-.175,-.403)),a[42]=vec2[42](vec2(-.4,-.403),vec2(-.4,-.415),vec2(-.388,-.415),vec2(-.388,-.415),vec2(-.377,-.415),vec2(-.377,-.403),vec2(-.377,-.403),vec2(-.377,-.415),vec2(-.365,-.415),vec2(-.365,-.415),vec2(-.353,-.415),vec2(-.353,-.403),vec2(-.318,-.403),vec2(-.318,-.392),vec2(-.329,-.392),vec2(-.329,-.392),vec2(-.341,-.392),vec2(-.341,-.403),vec2(-.341,-.403),vec2(-.341,-.415),vec2(-.329,-.415),vec2(-.282,-.403),vec2(-.282,-.415),vec2(-.27,-.415),vec2(-.258,-.403),vec2(-.258,-.415),vec2(-.246,-.415),vec2(-.222,-.403),vec2(-.21,-.403),vec2(-.21,-.392),vec2(-.222,-.403),vec2(-.21,-.403),vec2(-.21,-.415),vec2(-.175,-.403),vec2(-.175,-.392),vec2(-.186,-.392),vec2(-.186,-.392),vec2(-.198,-.392),vec2(-.198,-.403),vec2(-.198,-.403),vec2(-.198,-.415),vec2(-.186,-.415));"
               "for(int l=0;l<12;++l)"
                 "e=min(e,dsg(n[2*l],n[2*l+1],y));"
               "for(int l=0;l<14;++l)"
                 "e=min(e,dsp(a[3*l],a[3*l+1],a[3*l+2],y));"
               "r=mix(r,.8*c.xxx,B(25.)*smoothstep(.005,.002,e));"
             "}"
           "else"
             " if(iTime<33.)"
               "{"
                 "const vec2 a[176]=vec2[176](vec2(-.4,-.415),vec2(-.4,-.38),vec2(-.377,-.415),vec2(-.377,-.38),vec2(-.4,-.392),vec2(-.377,-.392),vec2(-.364,-.415),vec2(-.364,-.368),vec2(-.364,-.368),vec2(-.341,-.415),vec2(-.341,-.415),vec2(-.341,-.368),vec2(-.329,-.415),vec2(-.329,-.368),vec2(-.329,-.368),vec2(-.317,-.368),vec2(-.329,-.415),vec2(-.317,-.415),vec2(-.305,-.403),vec2(-.305,-.38),vec2(-.269,-.403),vec2(-.269,-.403),vec2(-.269,-.38),vec2(-.269,-.38),vec2(-.257,-.403),vec2(-.257,-.403),vec2(-.257,-.38),vec2(-.257,-.38),vec2(-.221,-.415),vec2(-.221,-.392),vec2(-.174,-.392),vec2(-.174,-.427),vec2(-.185,-.438),vec2(-.197,-.438),vec2(-.161,-.415),vec2(-.161,-.368),vec2(-.161,-.415),vec2(-.15,-.415),vec2(-.161,-.392),vec2(-.15,-.392),vec2(-.102,-.415),vec2(-.102,-.392),vec2(-.0663,-.403),vec2(-.0663,-.403),vec2(-.0663,-.38),vec2(-.0663,-.38),vec2(-.0541,-.403),vec2(-.0541,-.403),vec2(-.0541,-.38),vec2(-.0541,-.38),vec2(-.00625,-.403),vec2(-.00625,-.38),vec2(-.0179,-.392),vec2(.00542,-.392),vec2(.041,-.415),vec2(.041,-.392),vec2(.0533,-.415),vec2(.0533,-.392),vec2(.0772,-.415),vec2(.0772,-.368),vec2(.0772,-.415),vec2(.0888,-.415),vec2(.0772,-.392),vec2(.0888,-.392),vec2(.113,-.415),vec2(.113,-.392),vec2(.16,-.415),vec2(.16,-.392),vec2(.172,-.392),vec2(.172,-.403),vec2(.196,-.392),vec2(.196,-.415),vec2(.208,-.415),vec2(.22,-.415),vec2(.22,-.392),vec2(.231,-.392),vec2(.267,-.392),vec2(.255,-.392),vec2(.267,-.415),vec2(.255,-.415),vec2(.279,-.415),vec2(.279,-.368),vec2(.302,-.415),vec2(.302,-.403),vec2(.291,-.392),vec2(.279,-.392),vec2(.339,-.403),vec2(.339,-.403),vec2(.339,-.38),vec2(.339,-.38),vec2(.351,-.403),vec2(.351,-.403),vec2(.351,-.38),vec2(.351,-.38),vec2(.387,-.415),vec2(.387,-.392),vec2(.41,-.415),vec2(.41,-.403),vec2(.434,-.415),vec2(.434,-.403),vec2(.458,-.415),vec2(.469,-.415),vec2(.446,-.403),vec2(.469,-.403),vec2(.481,-.415),vec2(.481,-.392),vec2(.529,-.392),vec2(.517,-.392),vec2(.529,-.415),vec2(.517,-.415),vec2(.541,-.392),vec2(.541,-.403),vec2(.564,-.392),vec2(.564,-.415),vec2(.577,-.415),vec2(.577,-.392),vec2(.6,-.392),vec2(.6,-.403),vec2(.624,-.392),vec2(.624,-.427),vec2(.612,-.438),vec2(.6,-.438),vec2(.66,-.403),vec2(.66,-.403),vec2(.66,-.38),vec2(.66,-.38),vec2(.672,-.403),vec2(.672,-.403),vec2(.672,-.38),vec2(.672,-.38),vec2(.732,-.415),vec2(.72,-.415),vec2(.732,-.368),vec2(.72,-.368),vec2(.708,-.403),vec2(.708,-.38),vec2(.756,-.403),vec2(.756,-.368),vec2(.744,-.392),vec2(.767,-.392),vec2(.779,-.415),vec2(.779,-.392),vec2(.803,-.403),vec2(.803,-.368),vec2(.827,-.392),vec2(.851,-.392),vec2(.863,-.415),vec2(.863,-.38),vec2(.886,-.415),vec2(.886,-.38),vec2(.863,-.392),vec2(.886,-.392),vec2(.898,-.403),vec2(.898,-.368),vec2(.934,-.403),vec2(.934,-.368),vec2(.922,-.392),vec2(.946,-.392),vec2(.958,-.392),vec2(.981,-.392),vec2(1.01,-.415),vec2(1.01,-.368),vec2(.994,-.368),vec2(1.02,-.368),vec2(1.04,-.415),vec2(1.05,-.415),vec2(1.03,-.403),vec2(1.05,-.403),vec2(1.06,-.415),vec2(1.08,-.415),vec2(1.08,-.392),vec2(1.09,-.392),vec2(1.11,-.403),vec2(1.11,-.368),vec2(1.1,-.392),vec2(1.12,-.392)),n[204]=vec2[204](vec2(-.377,-.38),vec2(-.377,-.368),vec2(-.388,-.368),vec2(-.388,-.368),vec2(-.4,-.368),vec2(-.4,-.38),vec2(-.317,-.368),vec2(-.305,-.368),vec2(-.305,-.38),vec2(-.317,-.415),vec2(-.305,-.415),vec2(-.305,-.403),vec2(-.221,-.403),vec2(-.221,-.392),vec2(-.209,-.392),vec2(-.174,-.403),vec2(-.174,-.392),vec2(-.185,-.392),vec2(-.185,-.392),vec2(-.197,-.392),vec2(-.197,-.403),vec2(-.197,-.403),vec2(-.197,-.415),vec2(-.185,-.415),vec2(-.185,-.415),vec2(-.174,-.415),vec2(-.174,-.403),vec2(-.174,-.427),vec2(-.174,-.438),vec2(-.185,-.438),vec2(-.15,-.415),vec2(-.138,-.415),vec2(-.138,-.403),vec2(-.138,-.403),vec2(-.138,-.392),vec2(-.15,-.392),vec2(-.114,-.392),vec2(-.126,-.392),vec2(-.126,-.403),vec2(-.126,-.403),vec2(-.126,-.415),vec2(-.114,-.415),vec2(-.114,-.415),vec2(-.102,-.415),vec2(-.102,-.403),vec2(-.114,-.392),vec2(-.102,-.392),vec2(-.102,-.403),vec2(-.0179,-.415),vec2(-.00625,-.415),vec2(-.00625,-.403),vec2(-.00625,-.38),vec2(-.00625,-.368),vec2(.00542,-.368),vec2(.0293,-.392),vec2(.0177,-.392),vec2(.0177,-.403),vec2(.0177,-.403),vec2(.0177,-.415),vec2(.0293,-.415),vec2(.0293,-.415),vec2(.041,-.415),vec2(.041,-.403),vec2(.0293,-.392),vec2(.041,-.392),vec2(.041,-.403),vec2(.0533,-.403),vec2(.0533,-.392),vec2(.0649,-.392),vec2(.0888,-.415),vec2(.101,-.415),vec2(.101,-.403),vec2(.101,-.403),vec2(.101,-.392),vec2(.0888,-.392),vec2(.113,-.403),vec2(.113,-.392),vec2(.124,-.392),vec2(.148,-.392),vec2(.137,-.392),vec2(.137,-.403),vec2(.137,-.403),vec2(.137,-.415),vec2(.148,-.415),vec2(.148,-.415),vec2(.16,-.415),vec2(.16,-.403),vec2(.148,-.392),vec2(.16,-.392),vec2(.16,-.403),vec2(.172,-.403),vec2(.172,-.415),vec2(.184,-.415),vec2(.184,-.415),vec2(.196,-.415),vec2(.196,-.403),vec2(.22,-.415),vec2(.243,-.409),vec2(.22,-.403),vec2(.22,-.403),vec2(.196,-.398),vec2(.22,-.392),vec2(.255,-.392),vec2(.243,-.392),vec2(.243,-.403),vec2(.243,-.403),vec2(.243,-.415),vec2(.255,-.415),vec2(.302,-.403),vec2(.302,-.392),vec2(.291,-.392),vec2(.387,-.403),vec2(.387,-.392),vec2(.399,-.392),vec2(.399,-.392),vec2(.41,-.392),vec2(.41,-.403),vec2(.41,-.403),vec2(.41,-.392),vec2(.422,-.392),vec2(.422,-.392),vec2(.434,-.392),vec2(.434,-.403),vec2(.469,-.403),vec2(.469,-.392),vec2(.458,-.392),vec2(.458,-.392),vec2(.446,-.392),vec2(.446,-.403),vec2(.446,-.403),vec2(.446,-.415),vec2(.458,-.415),vec2(.481,-.403),vec2(.481,-.392),vec2(.493,-.392),vec2(.517,-.392),vec2(.505,-.392),vec2(.505,-.403),vec2(.505,-.403),vec2(.505,-.415),vec2(.517,-.415),vec2(.541,-.403),vec2(.541,-.415),vec2(.553,-.415),vec2(.553,-.415),vec2(.564,-.415),vec2(.564,-.403),vec2(.577,-.403),vec2(.577,-.392),vec2(.588,-.392),vec2(.6,-.403),vec2(.6,-.415),vec2(.612,-.415),vec2(.612,-.415),vec2(.624,-.415),vec2(.624,-.403),vec2(.624,-.427),vec2(.624,-.438),vec2(.612,-.438),vec2(.72,-.415),vec2(.708,-.415),vec2(.708,-.403),vec2(.708,-.38),vec2(.708,-.368),vec2(.72,-.368),vec2(.756,-.403),vec2(.756,-.415),vec2(.767,-.415),vec2(.779,-.403),vec2(.779,-.392),vec2(.791,-.392),vec2(.803,-.403),vec2(.803,-.415),vec2(.815,-.415),vec2(.886,-.38),vec2(.886,-.368),vec2(.875,-.368),vec2(.875,-.368),vec2(.863,-.368),vec2(.863,-.38),vec2(.898,-.403),vec2(.898,-.415),vec2(.91,-.415),vec2(.934,-.403),vec2(.934,-.415),vec2(.946,-.415),vec2(1.05,-.403),vec2(1.05,-.392),vec2(1.04,-.392),vec2(1.04,-.392),vec2(1.03,-.392),vec2(1.03,-.403),vec2(1.03,-.403),vec2(1.03,-.415),vec2(1.04,-.415),vec2(1.08,-.415),vec2(1.1,-.409),vec2(1.08,-.403),vec2(1.08,-.403),vec2(1.05,-.398),vec2(1.08,-.392),vec2(1.11,-.403),vec2(1.11,-.415),vec2(1.12,-.415));"
                 "for(int l=0;l<88;++l)"
                   "e=min(e,dsg(a[2*l],a[2*l+1],y));"
                 "for(int l=0;l<68;++l)"
                   "e=min(e,dsp(n[3*l],n[3*l+1],n[3*l+2],y));"
                 "r=mix(r,.8*c.xxx,B(30.)*smoothstep(.005,.002,e));"
               "}"
   "if(iTime>80.&&iTime<90.)"
     "{"
       "const vec2 n[114]=vec2[114](vec2(.0167,.1),vec2(.00833,.1),vec2(.0167,.133),vec2(.00833,.133),vec2(0.,.108),vec2(0.,.125),vec2(.0254,.1),vec2(.0254,.133),vec2(.0421,.1),vec2(.0421,.108),vec2(.0337,.117),vec2(.0254,.117),vec2(.0592,.1),vec2(.0675,.1),vec2(.0508,.108),vec2(.0675,.108),vec2(.0846,.1),vec2(.0929,.1),vec2(.0762,.108),vec2(.0929,.108),vec2(.102,.1),vec2(.102,.117),vec2(.119,.1),vec2(.127,.1),vec2(.127,.117),vec2(.135,.117),vec2(.17,.108),vec2(.17,.133),vec2(.161,.117),vec2(.178,.117),vec2(.229,.1),vec2(.229,.133),vec2(.263,.1),vec2(.263,.133),vec2(.229,.133),vec2(.246,.117),vec2(.246,.117),vec2(.263,.133),vec2(.28,.1),vec2(.28,.133),vec2(.271,.1),vec2(.288,.1),vec2(.271,.133),vec2(.288,.133),vec2(.313,.1),vec2(.305,.1),vec2(.313,.133),vec2(.305,.133),vec2(.297,.108),vec2(.297,.125),vec2(.356,.1),vec2(.356,.117),vec2(.365,.1),vec2(.365,.117),vec2(.381,.108),vec2(.381,.1),vec2(.398,.117),vec2(.407,.117),vec2(.398,.1),vec2(.407,.1),vec2(.407,.1),vec2(.407,.133),vec2(.432,.1),vec2(.441,.1),vec2(.449,.133),vec2(.441,.133),vec2(.475,.117),vec2(.466,.117),vec2(.475,.1),vec2(.466,.1),vec2(.483,.1),vec2(.483,.133),vec2(.5,.1),vec2(.5,.108),vec2(.492,.117),vec2(.483,.117),vec2(.509,.1),vec2(.509,.117),vec2(.525,.108),vec2(.525,.1),vec2(.551,.1),vec2(.551,.117),vec2(.56,.0833),vec2(.56,.117),vec2(.56,.117),vec2(.568,.117),vec2(.56,.1),vec2(.568,.1),vec2(.585,.0833),vec2(.585,.117),vec2(.585,.117),vec2(.593,.117),vec2(.585,.1),vec2(.593,.1),vec2(.61,.1),vec2(.619,.1),vec2(.619,.117),vec2(.627,.117),vec2(.652,.117),vec2(.652,.0917),vec2(.644,.0833),vec2(.636,.0833),vec2(.661,.108),vec2(.661,.117),vec2(.661,.125),vec2(.661,.125),vec2(.678,.1),vec2(.678,.117),vec2(.695,.108),vec2(.695,.133),vec2(.712,.1),vec2(.721,.1),vec2(.721,.117),vec2(.729,.117)),a[168]=vec2[168](vec2(.00833,.1),vec2(0.,.1),vec2(0.,.108),vec2(0.,.125),vec2(0.,.133),vec2(.00833,.133),vec2(.0421,.108),vec2(.0421,.117),vec2(.0337,.117),vec2(.0675,.108),vec2(.0675,.117),vec2(.0592,.117),vec2(.0592,.117),vec2(.0508,.117),vec2(.0508,.108),vec2(.0508,.108),vec2(.0508,.1),vec2(.0592,.1),vec2(.0929,.108),vec2(.0929,.117),vec2(.0846,.117),vec2(.0846,.117),vec2(.0762,.117),vec2(.0762,.108),vec2(.0762,.108),vec2(.0762,.1),vec2(.0846,.1),vec2(.102,.108),vec2(.102,.117),vec2(.11,.117),vec2(.127,.1),vec2(.144,.104),vec2(.127,.108),vec2(.127,.108),vec2(.11,.113),vec2(.127,.117),vec2(.17,.108),vec2(.17,.1),vec2(.178,.1),vec2(.195,.1),vec2(.187,.1),vec2(.187,.108),vec2(.187,.108),vec2(.187,.117),vec2(.195,.117),vec2(.195,.117),vec2(.203,.117),vec2(.203,.108),vec2(.203,.108),vec2(.203,.1),vec2(.195,.1),vec2(.305,.1),vec2(.297,.1),vec2(.297,.108),vec2(.297,.125),vec2(.297,.133),vec2(.305,.133),vec2(.347,.117),vec2(.339,.117),vec2(.339,.108),vec2(.339,.108),vec2(.339,.1),vec2(.347,.1),vec2(.347,.1),vec2(.356,.1),vec2(.356,.108),vec2(.347,.117),vec2(.356,.117),vec2(.356,.108),vec2(.373,.117),vec2(.381,.117),vec2(.381,.108),vec2(.365,.108),vec2(.365,.117),vec2(.373,.117),vec2(.398,.1),vec2(.39,.1),vec2(.39,.108),vec2(.39,.108),vec2(.39,.117),vec2(.398,.117),vec2(.441,.1),vec2(.449,.1),vec2(.449,.108),vec2(.449,.108),vec2(.449,.117),vec2(.441,.117),vec2(.441,.117),vec2(.432,.117),vec2(.432,.125),vec2(.432,.125),vec2(.432,.133),vec2(.441,.133),vec2(.466,.117),vec2(.458,.117),vec2(.458,.108),vec2(.458,.108),vec2(.458,.1),vec2(.466,.1),vec2(.5,.108),vec2(.5,.117),vec2(.492,.117),vec2(.517,.117),vec2(.525,.117),vec2(.525,.108),vec2(.509,.108),vec2(.509,.117),vec2(.517,.117),vec2(.542,.117),vec2(.534,.117),vec2(.534,.108),vec2(.534,.108),vec2(.534,.1),vec2(.542,.1),vec2(.542,.1),vec2(.551,.1),vec2(.551,.108),vec2(.542,.117),vec2(.551,.117),vec2(.551,.108),vec2(.568,.1),vec2(.576,.1),vec2(.576,.108),vec2(.576,.108),vec2(.576,.117),vec2(.568,.117),vec2(.593,.1),vec2(.602,.1),vec2(.602,.108),vec2(.602,.108),vec2(.602,.117),vec2(.593,.117),vec2(.619,.1),vec2(.635,.104),vec2(.619,.108),vec2(.619,.108),vec2(.602,.113),vec2(.619,.117),vec2(.652,.108),vec2(.652,.117),vec2(.644,.117),vec2(.644,.117),vec2(.636,.117),vec2(.636,.108),vec2(.636,.108),vec2(.636,.1),vec2(.644,.1),vec2(.644,.1),vec2(.652,.1),vec2(.652,.108),vec2(.652,.0917),vec2(.652,.0833),vec2(.644,.0833),vec2(.661,.108),vec2(.661,.1),vec2(.67,.1),vec2(.678,.108),vec2(.678,.117),vec2(.687,.117),vec2(.695,.108),vec2(.695,.1),vec2(.704,.1),vec2(.721,.1),vec2(.737,.104),vec2(.721,.108),vec2(.721,.108),vec2(.704,.113),vec2(.721,.117));"
       "for(int l=0;l<57;++l)"
         "e=min(e,dsg(n[2*l],n[2*l+1],y));"
       "for(int l=0;l<56;++l)"
         "e=min(e,dsp(a[3*l],a[3*l+1],a[3*l+2],y));"
       "r=mix(r,c.yyy,B(80.)*smoothstep(.005,.002,e));"
     "}"
   "v=vec4(r,1.);"
 "}"
 "void main()"
 "{"
   "mainImage(gl_FragColor,gl_FragCoord.xy);"
 "}";

#endif // GFX_H_
