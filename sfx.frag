#version 130

uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iVolume;

#define DEBUG_ONLY_DRUMS false
#define DEBUG_ONLY_MELO  false
#define DEBUG_ONLY_TRACK -1

#define PI radians(180.)
float _sin(float a) { return sin(2. * PI * mod(a,1.)); }
float _unisin(float a,float b) { return (.5*_sin(a) + .5*_sin((1.+b)*a)); }
float _sq(float a) { return sign(2.*fract(a) - 1.); }
float _squ(float a,float pwm) { return sign(2.*fract(a) - 1. + pwm); }
float _tri(float a) { return (4.*abs(fract(a)-.5) - 1.); }
float _saw(float a) { return (2.*fract(a) - 1.); }
float quant(float a,float div,float invdiv) { return floor(div*a+.5)*invdiv; }
float quanti(float a,float div) { return floor(div*a+.5)/div; }
float clip(float a) { return clamp(a,-1.,1.); }
float theta(float x) { return smoothstep(0., 0.01, x); }
float freqC1(float note){ return 32.7 * pow(2.,note/12.); }
float minus1hochN(int n) { return (1. - 2.*float(n % 2)); }

const float BPM = 32.; 
const float BPS = BPM/60.;
const float SPB = 60./BPM;

float doubleslope(float t, float a, float d, float s)
{
    return smoothstep(-.00001,a,t) - (1.-s) * smoothstep(0.,d,t-a);
}

float s_atan(float a) { return 2./PI * atan(a); }
float s_crzy(float amp) { return clamp( s_atan(amp) - 0.1*cos(0.9*amp*exp(amp)), -1., 1.); }
float squarey(float a, float edge) { return abs(a) < edge ? a : floor(4.*a+.5)*.25; } 

float TRISQ(float t, float f, int MAXN, float MIX, float INR, float NDECAY, float RES, float RES_Q)
{
    float ret = 0.;
    
    int Ninc = 4; // try this: leaving out harmonics...
    
    for(int N=0; N<=MAXN; N+=Ninc)
    {
        float mode     = 2.*float(N) + 1.;
        float inv_mode = 1./mode;         // avoid division? save table of Nmax <= 20 in some array or whatever
        float comp_TRI = (N % 2 == 1 ? -1. : 1.) * inv_mode*inv_mode;
        float comp_SQU = inv_mode;
        float filter_N = pow(1. + pow(float(N) * INR,2.*NDECAY),-.5) + RES * exp(-pow(float(N)*INR*RES_Q,2.));

        ret += (MIX * comp_TRI + (1.-MIX) * comp_SQU) * filter_N * _sin(mode * f * t);
    }
    
    return ret;
}

float QTRISQ(float t, float f, float QUANT, int MAXN, float MIX, float INR, float NDECAY, float RES, float RES_Q)
{
    return TRISQ(quant(t,QUANT,1./QUANT), f, MAXN, MIX, INR, NDECAY, RES, RES_Q);
}

float env_ADSR(float x, float L, float A, float D, float S, float R)
{
    float att = pow(x/A,8.);
    float dec = S + (1.-S) * exp(-(x-A)/D);
    float rel = (x < L-R) ? 1. : pow((L-x)/R,4.);

    return (x < A ? att : dec) * rel;
    
}

float macesaw(float t, float f, float CO, float Q, float det1, float det2, float res, float resQ)
{
    float s = 0.;
    float inv_CO = 1./CO;
    float inv_resQ = 1./resQ;
    float p = f*t;
        for(int N=1; N<=200; N++)
        {
            // saw
            float sawcomp = 2./PI * (1. - 2.*float(N % 2)) * 1./float(N);
            float filterN  = pow(1. + pow(float(N)*f*inv_CO,Q),-.5)
                     + res * exp(-pow((float(N)*f-CO)*inv_resQ,2.));
            
            if(abs(filterN*sawcomp) < 1e-6) break;
                
            if(det1 > 0. || det2 > 0.)
            {
                s += 0.33 * (_sin(float(N)*p) + _sin(float(N)*p*(1.+det1)) + _sin(float(N)*p*(1.+det2)));
            }
            else
            {
                s += filterN * sawcomp * _sin(float(N)*p);
            }
        }
    return s;
}

float maceskuh(float t, float f, float CO, float Q, float det1, float det2, float res, float resQ, float pw)
{
    float s = 0.;
    float inv_CO = 1./CO;
    float inv_resQ = 1./resQ;
    float p = f*t;
        for(int N=1; N<=200; N++)
        {
            // variable pulse wave: voll verrechnet, klingt aber geil =D
            float plscomp  = 1./float(N) * (1. + (2.*float(N%2)-1.)*_sin(pw)); 
            float filterN  = pow(1. + pow(float(N)*f*inv_CO,Q),-.5)
                     + res * exp(-pow((float(N)*f-CO)*inv_resQ,2.));
            
            if(abs(filterN*plscomp) < 1e-6) break;
                
            if(det1 > 0. || det2 > 0.)
            {
                s += 0.33 * (_sin(float(N)*p) + _sin(float(N)*p*(1.+det1)) + _sin(float(N)*p*(1.+det2)));
            }
            else
            {
                s += filterN * plscomp * _sin(float(N)*p);
            }
        }
    return 2.*s-1.;
}


float macesanderekuh(float t, float f, float CO, float Q, float det1, float det2, float res, float resQ, float pw)
{
    float s = 0.;
    float inv_CO = 1./CO;
    float inv_resQ = 1./resQ;
    float p = f*t;
        for(int N=1; N<=200; N++)
        {
            // varialbe pulse wave:
            float plscomp  = 1./(2.*PI*float(N)) * (minus1hochN(N)*_sin(pw*float(N)+.25) - 1.);
            float filterN  = pow(1. + pow(float(N)*f*inv_CO,Q),-.5)
                     + res * exp(-pow((float(N)*f-CO)*inv_resQ,2.));
            
            if(abs(filterN*plscomp) < 1e-6) break;
                
            if(det1 > 0. || det2 > 0.)
            {
                s += 0.33 * (_sin(float(N)*p) + _sin(float(N)*p*(1.+det1)) + _sin(float(N)*p*(1.+det2)));
            }
            else
            {
                s += filterN * plscomp * _sin(float(N)*p);
            }
        }
    return 2.*s-1.;
}


float freq_malformation(float t, float t_on, int vel, int Bsyn)
{
    t = t - min(t, t_on);
    
    float f = 80.;
    
    float fFM = 100.;
    float aFM = 0.01 * doubleslope(t, 0.8, 0.4, 0.5);
    float aFB = 0.000;
    
    float E = doubleslope(t, 0.2, 1., 0.);
    float r = _sin(t * f * (1. + aFM * _sin(t * fFM * (1. + aFB * _sin(t*fFM)))));
    return E * r;
}

float snare(float t, float t_on, float vel)
{
    // #define _tri(a) (4.*abs(fract(a)-.5) - 1.)
    t = t - min(t, t_on);
    float f1 = 6000.;
    float f2 = 800.;
    float f3 = 350.;
    float dec12 = 0.01;
    float dec23 = 0.01;
    float rel = 0.1;
    float snr = _tri(t * (f3 + (f1-f2)*smoothstep(-dec12,0.,-t)
                             + (f2-f3)*smoothstep(-dec12-dec23,-dec12,-t))) * smoothstep(-rel,-dec12-dec23,-t);
        
    //noise part
    float noise = 2. * fract(sin(t * 90.) * 45000.) * doubleslope(t,0.05,0.3,0.3);
   
    return vel * clamp(1.7 * (2. * snr + noise), -1.5, 1.5) * doubleslope(t,0.0,0.25,0.3);
}

float hut(float t, float t_on, float vel)
{
    t = t - min(t, t_on);
    float noise = fract(sin(t * 90.) * 45000.);
    noise = 1./(1.+noise);
    return vel * 2. * noise * doubleslope(t,0.,0.12,0.0);
    
    // might think of this one! - maybe tune length / pitch
    //float kick_blubb = (1.-exp(-1000.*t))*exp(-30.*t) * _sin((400.-200.*t)*t * _saw(4.*f*t));
}

float shake(float t, float t_on, float vel) // shaker is just some mod of hihat (hut)
{
    t = t - min(t, t_on);
    return vel * 0.5 * fract(sin(t * 90.) * 45000.) * doubleslope(t,0.03,0.15,0.15);
}

float hoskins_noise(float t) // thanks to https://www.shadertoy.com/view/4sjSW1 !
{
    float p = floor(t * (1500.0 * exp(-t*.100)));
    vec2 p2 = fract(vec2(p * 5.3983, p * 5.4427));
    p2 += dot(p2.yx, p2.xy + vec2(21.5351, 14.3137));
    return fract(p2.x * p2.y * 3.4337) * .5 * smoothstep(-.3,0.,-t);    
}

float facekick(float t, float t_on, float vel)
{
    t = t - min(t, t_on); // reset time to Bon event
    
    float f   = 50. + 150. * smoothstep(-0.12, 0., -t);
    float env = smoothstep(0.,0.015,t) * smoothstep(-0.08, 0., 0.16 - t);
    
    float kick_body = env * TRISQ(t, f, 3, 1., 0.8, 8., 4., 1.); // more heavy bass drum: increase reso parameters?
    
    float kick_click = 0.4 * step(t,0.03) * _sin(t*1100. * _saw(t*800.));
    
    float kick_blobb = (1.-exp(-1000.*t))*exp(-40.*t) * _sin((400.-200.*t)*t * _sin(1.*f*t));
    
    return vel * (kick_body + kick_blobb + 0.1*kick_click);
}

float hardkick(float t, float t_on, float vel)
{
    t = t - min(t, t_on); // reset time to Bon event
    
    float f   = 50. + 150. * smoothstep(-0.12, 0., -t);
    float env = smoothstep(0.,0.015,t) * smoothstep(-0.08, 0., 0.16 - t);
    
    float kick_body = env * TRISQ(t, f, 3, 1., 0.8, 8., 4., 1.); // more heavy bass drum: increase reso parameters?
    
    float kick_click = 0.4 * step(t,0.03) * _sin(t*1100. * _saw(t*800.));
    
    float kick_blobb = (1.-exp(-1000.*t))*exp(-40.*t) * _sin((400.-200.*t)*t * _sin(1.*f*t));
    
    return vel * (kick_body + kick_blobb + 0.1*kick_click);
}

float FM23(float _t, float B, float Boff, float note)
{
    // now model after DX reface
    // frequency ratio
    float FR2 = .9999;
    float FR1 = 6.;
    // level
    float LV2 = 0.5;
    float LV1 = 0.6;
    // feedback
    float FB2 = 0.1;
    float FB1 = 0.2;
    // actual algorithm (5?)
    float OP2 = _saw(_t * FR2 * (freqC1(note) + FB2 * _saw(_t * FR2 * freqC1(note))));
    float OP1 = LV1 * _tri(_t * FR1 * freqC1(note) + LV2 * OP2 + FB1 * _sin(_t * FR1 * freqC1(note)));

    return OP1;
}

float distsin(float t, float B, float Bon, float Boff, float note, int Bsyn)
{
    float Bprog = B-Bon;            // progress within Bar
    float Bproc = Bprog/(Boff-Bon); // relative progress
    float _t = SPB*(B - Bon); // reset time to Bon event
    float f = freqC1(note);

    float env = theta(B-Bon) * theta(Boff-B);
    float sound = clamp(1.1 * _sin(freqC1(note)*t), -0.999,0.999);

    if(Bsyn == -1) return 0.;
    
    if(Bsyn == 0)
        return env * sound; // test reasons: just give out something simple

    if(Bsyn == 23)
    {

        float kappa = 1.4;
    //sound = _unisin(t*freqC1(note),0.001)
      //  + _tri(t*2.*freqC1(note)*(1.+0.0001*_sin(0.5*Bproc))) * smoothstep(0.0,0.3, Bprog);

    //try some trumpet sound
    //float freq_LFO = 0.005 * smoothstep(0.,1.,_t) * _tri(_t);
    //sound = _saw(_t*freqC1(note)*(1.+freq_LFO) + 0.3 * smoothstep(0.,0.7,_t) * _saw(_t*3.01*freqC1(note)));
        
        float rel = .3;
        env = smoothstep(0., 0.001, B-Bon) * pow(smoothstep(rel, 0., B-Boff),2.5);

        env *= exp(-kappa*(B-Bon));
        
        note += 0.;

        float delta = 2.e-5 * (1.+.5*sin(2.*B));
        int filtN = 40;
        sound = 0.;
        
        for(int i = 0; i < filtN; i++)
            sound += 0.08 * FM23(_t + float(i)*delta, B, Boff, note);
        
        sound = .7*clamp(4.*sound,-1.,1.);
        
    }    
    
    if(Bsyn == 14) // test of mace-sq (matzeskuh)
    {
        env *= env_ADSR(Bprog,Boff-Bon,2.,0.,0.2,2.);
        
        env *= 0.5;
        
        float filterCO = 600. * env_ADSR(Bprog,Boff-Bon,2.,0.,1.,2.) + 40. * sqrt(f);

        sound += 0.3*macesanderekuh(t, f, filterCO, 30., 0.002, -0.01, 0.0, 0.0, 0.1);
        
        sound = clip(0.8*sound);
    }

    if(Bsyn == 15) // can I manage some horns, some day?
    {
        env *= env_ADSR(Bprog,Boff-Bon,3.,0.,0.2,2.);
        
        env *= 0.5;
        
        float filterCO = 100. * env_ADSR(Bprog,Boff-Bon,2.,0.,1.,2.) + 20. * sqrt(f);

        sound = 0.2*macesaw(t, f, filterCO, 50., -0.01-0.005*_sin(0.5*Bproc), 0.01-0.008*_sin(0.25*Bproc+.25), 0.0, 0.);

        sound += 0.3*maceskuh(t, f, filterCO, 30., 0.002, -0.01, 0.0, 0.0, 0.1+0.06*_sin(0.25*_t));

//        sound += 0.1*maceskuh(t, 2.*f, filterCO, 10., 0.004, -0.002, 0.0, 0.0, -0.1+0.03*_sin(2.*t));
        
        sound = clip(0.8*sound);
    }
    
    if(Bsyn == 16) // super saw pad
    {
        env *= env_ADSR(Bprog,Boff-Bon,1.5,2.,0.2,0.8);
        
        env *= 0.5;
        
        float filterQ = 100.;
        float filterCO = 200. + 100. * env_ADSR(Bprog,Boff-Bon,1.5,2.5,0.2,10.);

        sound = 0.9*macesaw(t, f, filterCO, filterQ, 0.010, 0.016, 0.0, 0.);
           // t --> quanti(t, 4096); // some nice wicked lo-bit shit
        
     
        sound = 0.4 * (0.5*_saw(t*f+.1+0.1*_sin(2.*t)) + _saw(t*f+.25) + _saw(t*f+0.2*_sin(0.5*t)));
        // lo-fi noise
        //sound = quanti(sound, 128.);
    }
    
    if(Bsyn == 17) // deftiges pad
    {
        env *= env_ADSR(Bprog,Boff-Bon,1.5,2.,0.2,0.8);
        
        env *= 0.5;
        
        float filterQ = 20.;
        float filterCO = 200. + 100. * env_ADSR(Bprog,Boff-Bon,1.5,2.5,0.2,10.);

        sound = 0.9*macesaw(t, f, filterCO, filterQ, 0.010, 0.020, 0.3, 3.);
           // t --> quanti(t, 4096); // some nice wicked lo-bit shit
        
        // lo-fi noise
        //sound = quanti(sound, 128.);
    }

    
    if(Bsyn == 10) // Alright, this is not mellow, but wayne.
    {
        sound = _tri(2.*f*t) + _tri(0.999*2.*f*_t*(1.+0.0001*_sin(4.*_t)));
        
        sound += 0.2 * _saw(5.01*f*t)
               + 0.2 * _saw(14.01*f*t);
        
        sound += _saw(t * (f + 1000. * smoothstep(-0.1,0.,-Bprog))) * smoothstep(-0.1,0.,-Bprog);
        
        sound *= 0.2;
    }
    
    // replace sound of Bsyn 0 by something squarey (with some automated PWM)
    if(Bsyn == 20)
    {
        sound = QTRISQ(t, freqC1(note), 1024., 16, 0.2, 8. + 20.*smoothstep(0.,0.3,Bproc), 1.2, 0.2, 0.1)
              + QTRISQ(t, 2.*freqC1(note), 1024., 16, 0.4, 8. + 20.*smoothstep(0.,0.1,Bproc), 1.2, 0.2, 0.1);
        
        sound = 0.8 * sound * exp(-4.0*Bprog);

        sound += 0.4*_sin(_t*f) * exp(-2.0*Bprog);
        sound += 0.4*_tri(_t*4.*f) * exp(-8.0*Bprog);
        
        sound += 0.6 * smoothstep(0.,1.2,Bproc) * (0.5 * _squ(_t*f,0.) + 0.5 * _squ(_t*f*.99,0.01));
        
        sound = clip(1.3*sound);
    }

    else if(Bsyn == 12)
    {
        float kappa = 1.4;
        float dist = 1.3;

        env *= exp(-kappa*(B-Bon));
        
        float freq_mod = 1. + 0.02 * exp(-10.*(B-Bon)); // edgy sound
        freq_mod = 1.;
            
        //sound = 0.7*sound + 0.3*clamp(2.*dist * _sin(freqC1(note+12)*(t-1e-3)), -0.999,0.999);
        
        sound  = 0.7 * _squ(freq_mod * freqC1(note)*t, -0.4 - 0.15*_sin(0.3*(B-Bon)));
        
        // add subbass
        sound += 0.5 * _unisin(freqC1(note)*t,0.05 * _sin(0.25*t));

        //reduce bit depth to 4
        //sound = quant(sound, 4., 0.25);
        
        // try something else, QTRISQ is my additive low-samplerate osci with improvized filter
        // QTRISQ(float t, float f, float QUANT, int MAXN, float MIX, float INR, float NDECAY, float RES, float RES_Q)
        sound = QTRISQ(t, freqC1(note), 1024., 16, 0.2, 8. + 20.*smoothstep(0.,0.3,Bproc), 1.2, 0.2, 0.1)
              + QTRISQ(t, 2.*freqC1(note), 1024., 16, 0.4, 8. + 20.*smoothstep(0.,0.1,Bproc), 1.2, 0.2, 0.1);
        
        float im = 2. + .5*_sin(0.21*_t);
        float y = _sin(4.*f*_t + im*_sin(0.25*f*Bproc) );
        
    }

    else if(Bsyn == 11)
    {
        float kappa = 1.4;
        float dist = 1.3;

        env *= exp(-kappa*(B-Bon));
        
        float freq_mod = 1. + 0.02 * exp(-10.*(B-Bon)); // edgy sound
        freq_mod = 1.;
            
        //sound = 0.7*sound + 0.3*clamp(2.*dist * _sin(freqC1(note+12)*(t-1e-3)), -0.999,0.999);
        
        sound  = 0.7 * _squ(freq_mod * freqC1(note)*t, -0.4 - 0.15*_sin(0.3*(B-Bon)));
        
        // add subbass
        sound += 0.5 * _unisin(freqC1(note)*t,0.05 * _sin(0.25*t));

        //reduce bit depth to 4
        //sound = quant(sound, 4., 0.25);

    }

    else if(Bsyn == 13) // inspired by some other shader. forgot which one. sorry.
    {
        float im = 2. + .5*_sin(0.21*_t);
        float y = _sin(4.*f*_t + im*_sin(0.25*f*Bproc) );
        
        sound = y;
    }

    
    else if(Bsyn == 4)
    {
        sound = QTRISQ(t, freqC1(note), 2048., 20, 0.2, 3. + 5.*smoothstep(0.,0.3,Bproc), 1.2, 0.2, 0.1)
              + QTRISQ(t, 2.*freqC1(note), 1024., 10, 0.4, 1. + 3.*smoothstep(0.,0.3,Bproc), 1.2, 0.2, 0.1);
        
        sound *= 0.3;
    }
    
    else if(Bsyn == 30)
    {
        float fir_d = 1.e-3;
        float fir_a = 2.e-2;
        sound = 0.;
        for(float fir = -fir_a; fir < fir_a; fir += fir_d)
        {
        sound += .8 * _sin(.5*f*(_t+fir))
              + TRISQ(t+fir, freqC1(note), 20, 0.2, .01 *smoothstep(0.,0.3,Bproc), 1.2, 0.2, 0.1)
              + .3*TRISQ(t+fir, 4.*freqC1(note), 10, .3, 1., 0.2, 2., 0.1);
        }
        sound = 0.5 * s_crzy(squarey(60.*sound,5.));
    }    
    
    else if(Bsyn == 1)
    {
        float dist = 1.55;

        env = smoothstep(0., 0.25, B-Bon) * (1.+0.01*sin(2.*PI*20.*t)); // * theta(Boff-B);

        float rel = .2;
        env *= smoothstep(rel, -rel, B-Boff);
       
        sound = 0.3*(clamp(dist * _sin(freqC1(note)*t)           , -0.999,0.999)
                   + clamp(dist * _sin(0.999*freqC1(note)*t+0.05), -0.999,0.999));

        sound += 0.7 * _unisin(freqC1(note)*t,0.05);    
    }
    
    else if(Bsyn == 2)
    {   
        float fFM = 0.33*f;
        float aFM = 0.3 * doubleslope(_t,1.,5.,0.01);
        float aFB = 0.02 * doubleslope(_t,2.,8.,0.00);
    
        env = doubleslope(B-Bon, 0.002, 2., 0.);
        sound = _sin(_t * f * (1. + aFM * _sin(_t * fFM * (1. + aFB * _sin(_t * fFM)))));
        
        //reduce bit depth to 16
        sound = quant(sound, 16., 0.0625);
        
        // rectify
        //sound = sign(sound);

        // try downsampling
        env = doubleslope(B-Bon, 0.01, 8., 0.);
        sound = 0.3 * _unisin(quant(_t,256.,0.004) * f * .25, 0.05);
    
        // check tuning
        sound += _sin(_t * f * 0.25);
    }    
    
    else if(Bsyn == 3)
    {
        float kappa = 1.4;
        env *= exp(-kappa*(B-Bon));
        sound = _unisin(t*freqC1(note),0.001)
              + _tri(t*2.*freqC1(note)*(1.+0.0001*_sin(0.5*Bproc))) * smoothstep(0.0,0.3, Bprog);

        //try some trumpet sound
        float freq_LFO = 0.005 * smoothstep(0.,1.,_t) * _tri(_t);
        sound = _saw(_t*freqC1(note)*(1.+freq_LFO) + 0.3 * smoothstep(0.,0.7,_t) * _saw(_t*3.01*freqC1(note)));
        
        // now model after DX reface
        // frequency ratio
        float FR2 = 0.998;
        float FR1 = 0.5;
           // level
        float LV2 = 0.5;
        float LV1 = 1.;
        // feedback
        float FB2 = 0.1;
        float FB1 = 0.6;
        // actual algorithm (5?)
        float OP2 = _saw(_t * FR2 * (freqC1(note) + FB2 * _saw(_t * FR2 * freqC1(note))));
        float OP1 = LV1 * _saw(_t * FR1 * freqC1(note) + LV2 * OP2 + FB1 * _sin(_t * FR1 * freqC1(note)));
        
        sound = 0.66 * OP1;
        
    }    
    
    return clamp(env,0.,1.) * clamp(sound, -1., 1.);
}


float mainSynth(float time)
{
    // START TRACK INFO
int NO_trks = 7;
int trk_sep[8] = int[8](0,3,9,12,14,16,18,20);
int trk_syn[7] = int[7](1,23,3,3,6,6,30);
float mod_on[20] = float[20](4.,6.,8.,4.25,5.25,6.25,7.25,8.25,8.75,8.,10.,11.,10.,11.,4.,6.,0.,2.,0.,2.);
float mod_off[20] = float[20](6.,8.,10.,4.75,5.75,6.75,7.75,8.75,10.25,10.,11.,12.,11.,12.,6.,8.,2.,4.,2.,4.);
int mod_ptn[20] = int[20](1,1,1,3,3,3,3,3,4,1,6,6,5,5,0,0,7,7,8,8);
float mod_transp[20] = float[20](0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.);
float inv_NO_tracks = .14285714285714285;
float max_mod_off = 4.;
int NO_ptns = 9;
int ptn_sep[10] = int[10](0,6,11,11,15,24,30,48,69,89);
float note_on[89] = float[89](0.,.125,.5,1.,1.25,1.5,0.,.5,1.,1.5,1.75,0.,.0625,.125,.1875,0.,.125,.1875,.25,.5625,.625,.6875,.75,1.,0.,.3125,.4375,.5,.8125,.9375,.03125,.0625,.125,.1875,.21875,.25,.3125,.375,.4375,.53125,.5625,.625,.6875,.71875,.75,.8125,.875,.9375,0.,0.,.125,.25,.25,.375,.5,.625,.75,.75,.875,1.,1.,1.125,1.25,1.25,1.375,1.5,1.625,1.75,1.875,0.,.125,.25,.375,.5,.625,.75,.875,1.,1.,1.125,1.125,1.25,1.375,1.375,1.5,1.625,1.625,1.75,1.875);
float note_off[89] = float[89](.125,.25,.625,1.125,1.375,1.625,.5,1.,1.5,1.75,2.,.03125,.125,.1875,.21875,.03125,.1875,.21875,.5625,.625,.6875,.71875,1.,1.25,.15625,.375,.5,.65625,.875,1.,.0625,.09375,.15625,.21875,.25,.3125,.375,.4375,.5,.5625,.59375,.65625,.71875,.75,.8125,.875,.9375,1.,.125,.125,.25,.375,.375,.5,.625,.75,.875,.875,1.,1.125,1.125,1.25,1.375,1.375,1.5,1.625,1.75,1.875,2.,.125,.25,.375,.5,.625,.75,.875,1.,1.125,1.125,1.25,1.25,1.375,1.5,1.5,1.625,1.75,1.75,1.875,2.);
float note_pitch[89] = float[89](24.,24.,25.,25.,27.,25.,24.,24.,23.,27.,26.,39.,39.,42.,38.,38.,42.,45.,47.,48.,44.,42.,42.,42.,44.,42.,41.,40.,39.,37.,20.,20.,20.,20.,20.,16.,26.,23.,25.,16.,16.,16.,16.,16.,20.,27.,13.,21.,34.,33.,36.,34.,33.,36.,37.,36.,34.,33.,37.,33.,23.,23.,33.,27.,25.,23.,25.,23.,25.,24.,24.,24.,25.,24.,24.,32.,24.,27.,34.,34.,27.,24.,20.,36.,24.,34.,22.,25.,24.);
float note_vel[89] = float[89](1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.,1.);
    

int drum_index = 6;
float drum_synths = 6.;

float max_release = 0.3;
    
float global_norm = .6;
float track_norm[7] = float[7](1.,1.,1.,1.,1.,.75,1.4);
    
    float r = 0.;
    float d = 0.;

    //which beat are we at?
    float BT = mod(BPS * time, max_mod_off); // mod for looping
    if(BT > max_mod_off) return r;

    // drum / sidechaining parameters
    float amt_drum = 0.3;
    float r_sidechain = 1.;
    float amt_sidechain = 0.99;
    float dec_sidechain = 0.6;

    float Bon = 0.;
    float Boff = 0.;

    for(int trk = 0; trk < NO_trks; trk++)
    {
        //if(DEBUG_ONLY_TRACK >=0 && trk != DEBUG_ONLY_TRACK) continue;

        int TLEN = trk_sep[trk+1] - trk_sep[trk];
        
        int _mod = TLEN;
        for(int i=0; i<TLEN; i++) if(BT < mod_off[(trk_sep[trk]+i)]) {_mod = i; break;}
        if(_mod == TLEN) continue;
        
        float B = BT - mod_on[trk_sep[trk]+_mod];
        
        int ptn = mod_ptn[trk_sep[trk]+_mod];
        int PLEN = ptn_sep[ptn+1] - ptn_sep[ptn];
        
        int _noteU = PLEN-1;
        for(int i=0; i<PLEN-1; i++) if(B < note_on[(ptn_sep[ptn]+i+1)]) {_noteU = i; break;}

        int _noteL = PLEN-1;
        for(int i=0; i<PLEN-1; i++) if(B <= note_off[(ptn_sep[ptn]+i)] + max_release) {_noteL = i; break;}
        
        for(int _note = _noteL; _note <= _noteU; _note++)
        {
            Bon    = note_on[(ptn_sep[ptn]+_note)];
            Boff   = note_off[(ptn_sep[ptn]+_note)];

            if(trk_syn[trk] == drum_index)
            {
                float Bdrum = mod(note_pitch[ptn_sep[ptn]+_note], drum_synths);
                float Bvel = 1.; //note_vel[(ptn_sep[ptn]+_note)] * pow(2.,mod_transp[_mod]/6.);

                float anticlick = 1.-exp(-1000.*(B-Bon));
                float _d = 0.;

                if(Bdrum < .5) // Sidechain
                {
                    r_sidechain = anticlick - amt_sidechain * theta(B-Bon) * smoothstep(-dec_sidechain,0.,Bon-B);
                }
                else if(Bdrum < 1.5) // Kick1
                {
                    _d = hardkick(B*SPB, Bon*SPB, Bvel) * theta(Boff-B);
                }
                else if(Bdrum < 2.5) // Kick2
                {
                    _d = facekick(B*SPB, Bon*SPB, Bvel) * theta(Boff-B);
                }
                else if(Bdrum < 3.5) // Snare1
                {
                    _d = snare(B*SPB, Bon*SPB, Bvel);
                }
                else if(Bdrum < 4.5) // HiHat
                {
                    _d = hut(B*SPB, Bon*SPB, Bvel);
                }                
                else if(Bdrum < 5.5) // Shake
                {
                    _d = shake(B*SPB, Bon*SPB, Bvel);
                }
                else if(Bdrum < 6.5) // ...
                {
                }          
                d += track_norm[trk] * _d;
            }
            else
            {
                r += track_norm[trk] * distsin(time, B, Bon, Boff,
                                               note_pitch[(ptn_sep[ptn]+_note)] + mod_transp[_mod], trk_syn[trk]);
            }

        }
    }

    d *= global_norm;
    r *= global_norm;

    if(DEBUG_ONLY_DRUMS) return 0.99 * clamp(d, -1., 1.);
    if(DEBUG_ONLY_MELO) return 0.99 * clamp(r, -1., 1.);

    r_sidechain = 1.;
    amt_drum = .5;

    float snd = s_atan((1.-amt_drum) * r_sidechain * r + amt_drum * d);

    return snd;
//    return sign(snd) * sqrt(abs(snd)); // eine von Matzes "besseren" Ideen
}

vec2 mainSound(float t)
{
    //maybe this works in enhancing the stereo feel
    float stereo_width = 0.1;
    float stereo_delay = 0.00001;
    
    //float comp_l = mainSynth(t) + stereo_width * mainSynth(t - stereo_delay);
    //float comp_r = mainSynth(t) + stereo_width * mainSynth(t + stereo_delay);
    
    //return vec2(comp_l * .99999, comp_r * .99999); 
    
    return vec2(mainSynth(t));
}


void main() 
{
   // compute time `t` based on the pixel we're about to write
   // the 512.0 means the texture is 512 pixels across so it's
   // using a 2 dimensional texture, 512 samples per row
   float t = iBlockOffset + ((gl_FragCoord.x-0.5) + (gl_FragCoord.y-0.5)*512.0)/iSampleRate;
    
//    t = mod(t, 4.5);
    
   // Get the 2 values for left and right channels
   vec2 y = iVolume * mainSound( t );

   // convert them from -1 to 1 to 0 to 65536
   vec2 v  = floor((0.5+0.5*y)*65536.0);

   // separate them into low and high bytes
   vec2 vl = mod(v,256.0)/255.0;
   vec2 vh = floor(v/256.0)/255.0;

   // write them out where 
   // RED   = channel 0 low byte
   // GREEN = channel 0 high byte
   // BLUE  = channel 1 low byte
   // ALPHA = channel 2 high byte
   gl_FragColor = vec4(vl.x,vh.x,vl.y,vh.y);
}

