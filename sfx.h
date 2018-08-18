/* File generated with Shader Minifier 1.1.5
 * http://www.ctrl-alt-test.fr
 */
#ifndef SFX_H_
# define SFX_H_

const char *sfx_frag =
 "#version 130\n"
 "uniform float iBlockOffset,iSampleRate,iVolume;"
 "const float s=acos(-1.);"
 "float f(float f)"
 "{"
   "return sin(2.*s*mod(f,1.));"
 "}"
 "float f(float s,float i)"
 "{"
   "return.5*f(s)+.5*f((1.+i)*s);"
 "}"
 "float p(float f)"
 "{"
   "return sign(2.*fract(f)-1.);"
 "}"
 "float p(float f,float i)"
 "{"
   "return sign(2.*fract(f)-1.+i);"
 "}"
 "float t(float f)"
 "{"
   "return 4.*abs(fract(f)-.5)-1.;"
 "}"
 "float n(float f)"
 "{"
   "return 2.*fract(f)-1.;"
 "}"
 "float f(float f,float i,float s)"
 "{"
   "return floor(i*f+.5)*s;"
 "}"
 "float n(float f,float i)"
 "{"
   "return floor(i*f+.5)/i;"
 "}"
 "float r(float f)"
 "{"
   "return clamp(f,-1.,1.);"
 "}"
 "float i(float f)"
 "{"
   "return smoothstep(0.,.01,f);"
 "}"
 "float w(float f)"
 "{"
   "return 32.7*pow(2.,f/12.);"
 "}"
 "float d(int f)"
 "{"
   "return 1.-2.*float(f%2);"
 "}"
 "const float m=17.5,g=m/60.,e=60./m;"
 "float d(float f,float i,float e,float s)"
 "{"
   "return smoothstep(-1e-05,i,f)-(1.-s)*smoothstep(0.,e,f-i);"
 "}"
 "float c(float f)"
 "{"
   "return 2./s*atan(f);"
 "}"
 "float a(float f)"
 "{"
   "return clamp(c(f)-.1*cos(.9*f*exp(f)),-1.,1.);"
 "}"
 "float a(float f,float i)"
 "{"
   "return abs(f)<i?f:floor(4.*f+.5)*.25;"
 "}"
 "float a(float s,float e,int i,float m,float c,float a,float x,float n)"
 "{"
   "float r=0.;"
   "int p=8;"
   "for(int g=0;g<=i;g+=p)"
     "{"
       "float v=2.*float(g)+1.,t=1./v,w=(g%2==1?-1.:1.)*t*t,l=t,y=pow(1.+pow(float(g)*c,2.*a),-.5)+x*exp(-pow(float(g)*c*n,2.));"
       "r+=(m*w+(1.-m)*l)*y*f(v*e*s);"
     "}"
   "return r;"
 "}"
 "float a(float f,float i,float s,float e,float m,float c)"
 "{"
   "float r=pow(f/s,8.),g=m+(1.-m)*exp(-(f-s)/e),t=f<i-c?1.:pow((i-f)/c,4.);"
   "return(f<s?r:g)*t;"
 "}"
 "float c(float e,float i,float m,float c,float g,float x,float v,float a)"
 "{"
   "float r=0.,t=1./m,n=1./a,y=i*e;"
   "for(int p=1;p<=200;p++)"
     "{"
       "float l=2./s*(1.-2.*float(p%2))/float(p),w=pow(1.+pow(float(p)*i*t,c),-.5)+v*exp(-pow((float(p)*i-m)*n,2.));"
       "if(abs(w*l)<1e-06)"
         "break;"
       "if(g>0.||x>0.)"
         "r+=.33*(f(float(p)*y)+f(float(p)*y*(1.+g))+f(float(p)*y*(1.+x)));"
       "else"
         " r+=w*l*f(float(p)*y);"
     "}"
   "return r;"
 "}"
 "float a(float f,float i,float s)"
 "{"
   "f=f-min(f,i);"
   "float p=6000.,g=800.,m=350.,e=.01,c=.01,w=.1,y=t(f*(m+(p-g)*smoothstep(-e,0.,-f)+(g-m)*smoothstep(-e-c,-e,-f)))*smoothstep(-w,-e-c,-f),l=2.*fract(sin(f*90.)*45000.)*d(f,.05,.3,.3);"
   "return s*clamp(1.7*(2.*y+l),-1.5,1.5)*d(f,0.,.25,.3);"
 "}"
 "float c(float f,float i,float s)"
 "{"
   "f=f-min(f,i);"
   "float e=fract(sin(f*90.)*45000.);"
   "e=1./(1.+e);"
   "return s*2.*e*d(f,0.,.12,0.);"
 "}"
 "float d(float f,float i,float s)"
 "{"
   "return f=f-min(f,i),s*.5*fract(sin(f*90.)*45000.)*d(f,.03,.15,.15);"
 "}"
 "float v(float f)"
 "{"
   "float i=floor(f*(1500.*exp(-f*.1)));"
   "vec2 s=fract(vec2(i*5.3983,i*5.4427));"
   "s+=dot(s.yx,s.xy+vec2(21.5351,14.3137));"
   "return fract(s.x*s.y*3.4337)*.5*smoothstep(-.3,0.,-f);"
 "}"
 "float i(float i,float s,float p)"
 "{"
   "i=i-min(i,s);"
   "float e=50.+150.*smoothstep(-.12,0.,-i),w=smoothstep(0.,.015,i)*smoothstep(-.08,0.,.16-i),g=w*a(i,e,3,1.,.8,8.,4.,1.),y=.4*step(i,.03)*f(i*1100.*n(i*800.)),t=(1.-exp(-1000.*i))*exp(-40.*i)*f((400.-200.*i)*i*f(e*i));"
   "return p*(g+t+.1*y);"
 "}"
 "float n(float i,float s,float p)"
 "{"
   "i=i-min(i,s);"
   "float e=60.+150.*smoothstep(-.3,0.,-i),w=smoothstep(0.,.01,i)*smoothstep(-.1,.2,.3-i),g=w*.1*a(i,e,100,1.,1.,.1,16.,10.);"
   "g+=.7*(smoothstep(0.,.01,i)*smoothstep(-.2,.2,.3-i))*f(i*e*.5);"
   "float r=1.5*step(i,.05)*f(i*5000.*n(i*1000.));"
   "r=c(40.*(1.-exp(-1000.*i))*exp(-80.*i)*f((1200.-1000.*sin(1000.*i*sin(30.*i)))*i));"
   "float m=a(10.*(1.-exp(-1000.*i))*exp(-30.*i)*f((300.-300.*i)*i));"
   "return p*2.*clamp(g+m+r,-1.5,1.5);"
 "}"
 "float c(float s,float p,float m,float g,float r,int l)"
 "{"
   "float y=p-m,t=y/(g-m),n=e*(p-m),x=w(r),v=i(p-m)*i(g-p),b=clamp(1.1*f(w(r)*s),-.999,.999);"
   "if(l==-1)"
     "return 0.;"
   "if(l==0)"
     "return v*b;"
   "else"
     " if(l==61)"
       "{"
         "v=smoothstep(0.,.0002,y)*smoothstep(.05,0.,p-g);"
         "float d=20.,o=2000.+1000.*a(16.*y,g-m,1.5,2.5,.2,10.);"
         "b=.9*c(s,.5*x,o,d,.01,.02,0.,0.)+.4*c(s,.499*x,o,d,.01,.02,0.,0.);"
       "}"
     "else"
       " if(l==84)"
         "{"
           "v=smoothstep(0.,1e-05,y)*smoothstep(.05,0.,p-g);"
           "float d=100.,o=1500.+1000.*smoothstep(0.,.25,t);"
           "b=.9*c(s,x,o,d,.01,.02,.3,3.);"
         "}"
   "if(l==76)"
     "b=a(s,w(r),160,.4+.3*sin(32.*p*(1.+sin(3.*p))),.1+.1*sin(24.*p),.5,.2,.1);"
   "return clamp(v,0.,1.)*clamp(b,-1.,1.);"
 "}"
 "float o(float f)"
 "{"
   "int s=6,m[7]=int[7](0,19,34,46,49,85,92),t[6]=int[6](76,84,61,61,6,6);"
   "float p[92]=float[92](0.,2.,6.,8.,10.,12.,14.,16.,18.,20.,21.,22.,23.,25.,27.,29.,31.,33.,35.,4.,6.,8.,10.,14.,16.,18.,20.,21.,22.,27.,29.,31.,33.,35.,6.,8.,10.,12.,16.,18.,23.,25.,27.,29.,31.,33.,6.,8.,10.,12.,14.,16.,18.,18.5,19.,19.5,20.,20.5,21.,21.5,22.,22.5,23.,23.5,24.,24.5,25.,25.5,26.,26.5,27.,27.5,28.,28.5,29.,29.5,30.,30.5,31.5,32.,32.5,33.,33.5,34.,34.5,24.5,26.5,28.5,30.5,32.5,34.5,36.),w[92]=float[92](2.,4.,8.,10.,12.,14.,16.,18.,20.,21.,22.,23.,25.,27.,29.,31.,33.,35.,37.,6.,8.,10.,12.,16.,18.,20.,21.,22.,23.,29.,31.,33.,35.,37.,8.,10.,12.,14.,18.,20.,25.,27.,29.,31.,33.,35.,8.,10.,12.,14.,16.,18.,18.5,19.,19.5,20.,20.5,21.,21.5,22.,22.5,23.,23.5,24.,24.5,25.,25.5,26.,26.5,27.,27.5,28.,28.5,29.,29.5,30.,30.5,31.,32.,32.5,33.,33.5,34.,34.5,35.,25.,27.,29.,31.,33.,35.,36.5);"
   "int r[92]=int[92](1,1,5,5,5,6,6,7,8,9,9,9,10,10,10,10,10,10,0,2,2,2,2,7,7,2,9,9,9,10,10,10,10,16,3,3,3,6,7,3,11,11,11,11,11,11,4,4,4,12,12,12,13,13,13,13,13,13,13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,15,15,15,15,15,15,15);"
   "float v[92]=float[92](24.,24.,0.,0.,0.,12.,12.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,-12.,0.,0.,-12.,-12.,0.,0.,0.,0.,0.,0.,24.,12.,0.,0.,0.,0.,0.,12.,12.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.),y=.166667,l=37.;"
   "int x=6;"
   "float o=11.;"
   "int b=17,u[18]=int[18](0,2,12,43,74,107,132,160,183,208,222,254,310,347,357,365,367,373);"
   "float h[373]=float[373](0.,.5,0.,.25,.375,.5,.625,1.,1.25,1.375,1.5,1.625,0.,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,0.,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,0.,0.,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.75,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,.0625,.125,.1875,.3125,.375,.4375,.5625,.625,.6875,.8125,.875,.9375,1.0625,1.125,1.1875,1.3125,1.375,1.4375,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,0.,0.,0.,0.,.375,.375,.375,.375,.5,.5,.5,.5,.875,.875,.875,.875,1.,1.,1.,1.,1.375,1.375,1.375,1.375,1.5,1.5,1.5,1.5,.0625,.125,.1875,.25,.3125,.375,.5625,.625,.6875,.75,.8125,.875,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.5,1.5625,1.625,1.8125,1.875,.0625,.125,.1875,.3125,.375,.4375,.5625,.625,.6875,.8125,.875,.9375,1.0625,1.125,1.1875,1.3125,1.375,1.4375,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,.0625,.125,.1875,.25,.3125,.375,.4375,.5625,.625,.6875,.75,.8125,.875,.9375,0.,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.75,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,.0625,.0625,.0625,.125,.1875,.1875,.1875,.3125,.3125,.3125,.375,.4375,.4375,.4375,.5625,.5625,.5625,.625,.6875,.6875,.6875,.8125,.8125,.8125,.875,.9375,.9375,.9375,1.0625,1.0625,1.0625,1.125,1.1875,1.1875,1.1875,1.3125,1.3125,1.3125,1.375,1.4375,1.4375,1.4375,1.5625,1.5625,1.5625,1.625,1.6875,1.6875,1.6875,1.8125,1.8125,1.8125,1.875,1.9375,1.9375,1.9375,0.,.03125,.09375,.09375,.125,.28125,.3125,.34375,.375,.5,.53125,.5625,.59375,.75,.78125,.84375,.875,1.,1.,1.03125,1.0625,1.09375,1.28125,1.3125,1.34375,1.375,1.5,1.5,1.53125,1.5625,1.59375,1.6875,1.75,1.78125,1.8125,1.84375,1.875,0.,0.,.0625,.125,.1875,.25,.25,.3125,.375,.4375,0.,.0625,.125,.1875,.25,.3125,.375,.4375,.3125,.4375,0.,0.,0.,.5,.5,.5),k[373]=float[373](.5,2.,.25,.5,.4375,1.,.6875,1.25,1.5,1.4375,2.,1.6875,.25,.125,.1875,.25,.5,.375,.4375,.5,1.,.625,.6875,.75,.875,.9375,1.,1.25,1.125,1.1875,1.25,1.5,1.375,1.4375,1.5,2.,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,.25,.125,.1875,.25,.5,.375,.4375,.5,1.,.625,.6875,.75,.875,.9375,1.,1.25,1.125,1.1875,1.25,1.5,1.375,1.4375,1.5,2.,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,.0625,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.75,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,.125,.1875,.25,.375,.4375,.5,.625,.6875,.75,.875,.9375,1.,1.125,1.1875,1.25,1.375,1.4375,1.5,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,.25,.25,.25,.25,.5,.5,.5,.5,.75,.75,.75,.75,1.,1.,1.,1.,1.25,1.25,1.25,1.25,1.5,1.5,1.5,1.5,1.75,1.75,1.75,1.75,.09375,.1875,.25,.28125,.34375,.5,.59375,.6875,.75,.78125,.84375,1.,1.09375,1.1875,1.25,1.28125,1.34375,1.5,1.5625,1.59375,1.75,1.84375,2.,.125,.1875,.25,.375,.4375,.5,.625,.6875,.75,.875,.9375,1.,1.125,1.1875,1.25,1.375,1.4375,1.5,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,.125,.1875,.25,.3125,.375,.4375,.5,.625,.6875,.75,.8125,.875,.9375,1.,.0625,.125,.1875,.25,.3125,.375,.4375,.5,.5625,.625,.6875,.75,.8125,.875,.9375,1.,1.0625,1.125,1.1875,1.25,1.3125,1.375,1.4375,1.5,1.5625,1.625,1.6875,1.75,1.8125,1.875,1.9375,2.,.125,.125,.125,.1875,.25,.25,.25,.375,.375,.375,.4375,.5,.5,.5,.625,.625,.625,.6875,.75,.75,.75,.875,.875,.875,.9375,1.,1.,1.,1.125,1.125,1.125,1.1875,1.25,1.25,1.25,1.375,1.375,1.375,1.4375,1.5,1.5,1.5,1.625,1.625,1.625,1.6875,1.75,1.75,1.75,1.875,1.875,1.875,1.9375,2.,2.,2.,.0625,.0625,.125,.15625,.1875,.3125,.375,.375,.4375,.5625,.5625,.625,.625,.8125,.8125,.875,.9375,1.0625,1.0625,1.0625,1.125,1.125,1.3125,1.375,1.375,1.4375,1.5625,1.5625,1.5625,1.625,1.625,1.75,1.8125,1.8125,1.875,1.875,1.9375,.0625,.0625,.125,.1875,.25,.3125,.3125,.375,.4375,.5,.125,.125,.25,.25,.375,.375,.5,.5,.375,.5,.5,.46875,.46875,1.875,2.,1.96875),F[373]=float[373](22.,17.,29.,36.,29.,37.,32.,29.,35.,32.,31.,24.,29.,17.,24.,17.,36.,24.,29.,24.,37.,25.,33.,25.,25.,34.,33.,30.,17.,24.,17.,35.,19.,26.,19.,31.,13.,24.,13.,10.,12.,13.,15.,41.,29.,36.,29.,48.,36.,41.,36.,49.,37.,45.,37.,37.,46.,45.,42.,29.,36.,29.,47.,31.,38.,31.,43.,25.,36.,25.,22.,24.,25.,29.,36.,36.,41.,48.,36.,37.,41.,44.,49.,37.,41.,44.,49.,34.,37.,41.,43.,44.,36.,37.,41.,43.,36.,44.,36.,43.,32.,31.,29.,36.,32.,31.,25.,17.,24.,17.,24.,29.,24.,25.,33.,25.,25.,34.,33.,17.,24.,17.,19.,26.,19.,13.,24.,13.,10.,12.,13.,15.,17.,32.,41.,48.,28.,12.,36.,43.,34.,6.,37.,49.,34.,10.,26.,53.,15.,56.,43.,27.,10.,26.,34.,53.,12.,53.,36.,31.,48.,46.,44.,43.,41.,43.,44.,43.,41.,40.,37.,40.,44.,46.,44.,43.,41.,44.,48.,49.,43.,41.,43.,29.,36.,29.,36.,41.,36.,37.,45.,37.,37.,46.,45.,29.,36.,29.,31.,38.,31.,25.,36.,25.,26.,28.,28.,28.,17.,29.,17.,20.,24.,28.,23.,17.,24.,17.,20.,25.,32.,27.,17.,17.,17.,17.,24.,24.,24.,24.,17.,17.,17.,17.,16.,16.,16.,16.,17.,17.,17.,17.,19.,19.,19.,19.,20.,20.,20.,20.,23.,23.,23.,23.,41.,17.,36.,29.,46.,17.,32.,36.,24.,41.,28.,40.,43.,24.,44.,41.,17.,24.,44.,36.,17.,43.,36.,16.,23.,43.,36.,16.,44.,36.,17.,24.,44.,36.,17.,43.,46.,19.,12.,46.,40.,19.,48.,44.,25.,22.,48.,44.,25.,52.,43.,29.,24.,52.,43.,29.,23.,26.,26.,23.,25.,26.,23.,26.,25.,25.,26.,23.,26.,23.,26.,26.,25.,25.,23.,26.,23.,26.,26.,23.,26.,25.,25.,23.,26.,23.,26.,23.,23.,26.,23.,26.,25.,23.,37.,35.,35.,34.,37.,34.,35.,35.,37.,23.,26.,23.,3.,23.,26.,23.,3.,25.,25.,22.,41.,53.,53.,17.,44.),C=1.,V=.3,S=.6,R[7]=float[7](1.,1.,.9,.9,1.,.9,1.7),O=0.,B=0.,Z=mod(g*f,l);"
   "if(Z>l)"
     "return O;"
   "float Y=.3,X=1.,W=.99,U=.6,T=0.,Q=0.;"
   "for(int P=0;P<s;P++)"
     "{"
       "int N=m[P+1]-m[P],M=N;"
       "for(int L=0;L<N;L++)"
         "if(Z<w[m[P]+L])"
           "{"
             "M=L;"
             "break;"
           "}"
       "if(M==N)"
         "continue;"
       "float L=Z-p[m[P]+M];"
       "int K=r[m[P]+M],J=u[K+1]-u[K],I=J-1;"
       "for(int H=0;H<J-1;H++)"
         "if(L<h[u[K]+H+1])"
           "{"
             "I=H;"
             "break;"
           "}"
       "int H=J-1;"
       "for(int G=0;G<J-1;G++)"
         "if(L<=k[u[K]+G]+V)"
           "{"
             "H=G;"
             "break;"
           "}"
       "for(int G=H;G<=I;G++)"
         "{"
           "T=h[u[K]+G];"
           "Q=k[u[K]+G];"
           "if(t[P]==x)"
             "{"
               "float E=mod(F[u[K]+G],o),D=1.,A=1.-exp(-1000.*(L-T)),z=0.;"
               "if(E<.01)"
                 "X=A-W*i(L-T)*smoothstep(-U,0.,T-L);"
               "else"
                 " if(E<1.01)"
                   "X=A-W*i(L-T)*smoothstep(-U,0.,T-L),X*=.5,z=i(L*e,T*e,D);"
                 "else"
                   " if(E<2.01)"
                     "X=A-W*i(L-T)*smoothstep(-U,0.,T-L),X*=.5,z=n(L*e,T*e,D);"
                   "else"
                     " if(E<3.01)"
                       "z=a(L*e,T*e,D);"
                     "else"
                       " if(E<4.01)"
                         "z=c(L*e,T*e,D);"
                       "else"
                         " if(E<5.01)"
                           "z=d(L*e,T*e,D);"
                         "else"
                           " if(E<6.01)"
                             ";"
               "B+=R[P]*z;"
             "}"
           "else"
             " O+=R[P]*c(f,L,T,Q,F[u[K]+G]+v[M],t[P]);"
         "}"
     "}"
   "B*=S;"
   "O*=S;"
   "X=1.;"
   "Y=.5;"
   "float L=c((1.-Y)*X*O+Y*B);"
   "return c(L);"
 "}"
 "vec2 l(float f)"
 "{"
   "float i=.1,t=1e-05;"
   "return vec2(o(f));"
 "}"
 "void main()"
 "{"
   "float f=iBlockOffset+(gl_FragCoord.x-.5+(gl_FragCoord.y-.5)*512.)/iSampleRate;"
   "vec2 i=iVolume*l(f),e=floor((.5+.5*i)*65536.),g=mod(e,256.)/255.,t=floor(e/256.)/255.;"
   "gl_FragColor=vec4(g.x,t.x,g.y,t.y);"
 "}";

#endif // SFX_H_
