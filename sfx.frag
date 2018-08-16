#version 130

uniform float iBlockOffset;
uniform float iSampleRate;
uniform float iVolume;

vec2 mainSound(float t)
{
    //maybe this works in enhancing the stereo feel
    float stereo_width = 0.1;
    float stereo_delay = 0.00001;
    
//     float comp_l = mainSynth(t) + stereo_width * mainSynth(t - stereo_delay);
//     float comp_r = mainSynth(t) + stereo_width * mainSynth(t + stereo_delay);
    
//     return vec2(comp_l * .99999, comp_r * .99999); 
    return vec2(tanh(sin(acos(-1.)*440.*t)));
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

