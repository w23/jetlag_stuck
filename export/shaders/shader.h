/* File generated with Shader Minifier 1.1.4
 * http://www.ctrl-alt-test.fr
 */
#ifndef SHADER_H_
# define SHADER_H_

const char *shader_glsl =
 "uniform float t;"
 "uniform sampler2D Tex;"
 "vec3 v=vec3(0.,.01,1.);"
 "vec2 i=vec2(1920.,1080.);"
 "float f(float i)"
 "{"
   "return fract(sin(i)*46347.4);"
 "}"
 "float s(vec2 i)"
 "{"
   "return f(dot(i,vec2(79.5325,31.4328)));"
 "}"
 "float n(vec2 i)"
 "{"
   "vec2 f=floor(i);"
   "i-=f;"
   "i*=i*(3.-2.*i);"
   "return mix(mix(s(f+v.xx),s(f+v.zx),i.x),mix(s(f+v.xz),s(f+v.zz),i.x),i.y);"
 "}"
 "float r(vec2 i)"
 "{"
   "if(t<1016.)"
     "return 0.;"
   "if(i.x<0.||i.x>32.)"
     "return 0.;"
   "if(i.y>2.8)"
     "return 0.;"
   "i-=2.;"
   "float v=floor(i.y/1.5);"
   "if(t<1152.)"
     "{"
       "if(v==-1.)"
         ";"
       "else"
         " if(v==-2.)"
           "v=t<1056.?0.:-2.-floor(t-1056.)/8.;"
         "else"
           " return 0.;"
     "}"
   "else"
     " v-=16.;"
   "i.y=v*1.5+mod(i.y,1.5);"
   "return texture2D(Tex,i/32.).x;"
 "}"
 "float m(vec3 i)"
 "{"
   "return min(i.y+2.,length(i)-2.);"
 "}"
 "vec3 x(vec3 i)"
 "{"
   "return normalize(vec3(m(i+v.yxx),m(i+v.xyx),m(i+v.xxy))-m(i));"
 "}"
 "float f(vec3 i,vec3 v,float f,float x,float y)"
 "{"
   "for(float e=0.;e<y;++e)"
     "{"
       "float s=m(i+v*f);"
       "f+=s;"
       "if(s<.001||f>x)"
         "break;"
     "}"
   "return f;"
 "}"
 "void main()"
 "{"
   "vec2 m=gl_FragCoord.xy/i*2.-1.;"
   "m.x*=i.x/i.y;"
   "vec3 e=vec3(0.),y,d,z,p;"
   "vec2 a;"
   "float c,b,l=fract(t+m.x),T=floor(t/16.);"
   "vec3 u=normalize(vec3(-.5,.5,.4));"
   "u=normalize(vec3(1.));"
   "float o=.4,g=1.+fract(t/32.);"
   "vec3 h=vec3(0.,-1.,0.),k=vec3(0.,1.9,5.);"
   "k.z+=4.*fract(t/32.);"
   "float w=0.,q=0.,F=0.,D=0.,C=0.;"
   "if(t<304.)"
     "q=0.,F=0.,D=0.,o=.01,g=2.,h=vec3(0.,0.,0.),k=vec3(0.,1.,5.+2.*fract(t/16.));"
   "else"
     " h=vec3(sin(T*3.)*4.,sin(T*4.),sin(T*5.)*4.),q=1.,F=0.,D=0.,o=.2,g=2.,k=vec3(sin(T)*10.,2.+2.*sin(T*7.),10.*cos(T*3.));"
   "if(t>362.)"
     "{"
       "q=1.;"
       "D=3.;"
       "F=0.;"
       "float Z=fract((t-362.)/128.);"
       "h=vec3(0.,0.,0.);"
       "k=vec3(0.,1.,5.+32.*Z);"
       "g=1.+4.*Z;"
     "}"
   "if(t>490.)"
     "F=1.,o=.4;"
   "if(t>580.)"
     "h=vec3(sin(T*3.)*4.,sin(T*4.),sin(T*5.)*4.),k=vec3(sin(T)*10.,2.+2.*sin(T*7.),10.*cos(T*3.));"
   "if(t>704.)"
     "F=2.,q=2.,o=.4;"
   "if(t>832.)"
     "{"
       "q=mod(floor(t/8.),4.);"
       "F=mod(floor(t/4.),5.);"
       "D=mod(floor(t/6.),4.);"
       "if(q==3.)"
         "D=1.;"
       "float Z=fract(t/16.),Y=sin(T),X=sin(T+3.);"
       "k=h+mix(vec3(cos(Y)*10.,2.+2.*sin(T*2.),sin(Y)*10.),vec3(cos(X)*10.,2.+2.*sin(T*3.),sin(X)*10.),Z*Z);"
       "k*=max(1.,8./length(k));"
       "h=vec3(sin(T*3.)*2.,sin(T*4.),sin(T*5.)*2.);"
       "C=sin(T*17.);"
     "}"
   "if(t>1016.)"
     "{"
       "h=vec3(5.,1.,0.);"
       "k.z=abs(k.z);"
       "k*=max(1.,8./length(k));"
       "float Z=fract(t/64.);"
       "g=1.+2.*Z;"
       "k=vec3(8.+10.*sin(T+t/16.),2.,10.);"
     "}"
   "if(t>1152.)"
     "q=0.,k=vec3(10.,2.,10.),C=0.;"
   "if(t>1184.)"
     "D=3.;"
   "if(t>1216.)"
     "F=1.;"
   "w+=length(h-k);"
   "for(float Z=0.;Z<32.;++Z)"
     "{"
       "d=normalize(cross(y=normalize(h-k),vec3(C,1.,0.)));"
       "vec3 Y=normalize(cross(d,y));"
       "mat3 X=mat3(d,Y,-y);"
       "if(t<266.)"
         "{"
           "d=X*normalize(vec3(m/g,-1.)*w-y);"
           "y=X*y+k;"
           "float W=100.,V=f(y,d,0.,W,40.);"
           "if(V<W)"
             "{"
               "z=y+d*V;"
               "vec3 U=x(z);"
               "if(T<2.)"
                 "e=vec3(1.);"
               "else"
                 " if(T<4.)"
                   "e=vec3(V/W);"
                 "else"
                   " if(T<6.)"
                     "e=z;"
                   "else"
                     " if(T<8.)"
                       "e=fract(z);"
                     "else"
                       " if(T<10.)"
                         "e=U;"
                       "else"
                         "{"
                           "e+=.8*vec3(max(0.,dot(U,u)));"
                           "if(T>14.)"
                             "e*=step(5.,f(z,u,.1,5.,20.));"
                           "if(T>12.)"
                             "e=vec3(.01)+.5*e;"
                         "}"
               "e*=32.;"
             "}"
           "break;"
         "}"
       "vec2 W=vec2(0.);"
       "float V=f(l+=y.x)*6.2831;"
       "y=vec3(vec2(cos(V),sin(V))*sqrt(f(l+=m.y))*o,0.);"
       "d=X*normalize(vec3((m+W)/g,-1.)*w-y);"
       "y=X*y+k;"
       "float U=1.,S=f(l+=z.x);"
       "vec3 R=clamp(abs(fract(vec3(3.,2.,1.)/3.+S)*6.-3.)-1.,0.,1.);"
       "for(float Q=0.;Q<6.;++Q)"
         "{"
           "vec3 P=vec3(0.),O=vec3(.8);"
           "float N=1.;"
           "vec2 M=vec2(1.,1.);"
           "b=1e+06;"
           "c=0.;"
           "if(d.y<0.)"
             "b=(-2.-y.y)/d.y,z=y+d*b,p=v.xzx,a=z.xz,c=1.;"
           "float L=dot(-y,d),K=L*L-dot(-y,-y)+4.;"
           "if(K>=0.)"
             "{"
               "K=sqrt(K);"
               "float J=L-K,I=L+K;"
               "if(J<0.)"
                 "J=1e+06;"
               "if(I<0.)"
                 "I=1e+06;"
               "float H=min(J,I);"
               "if(H<b)"
                 "b=H,z=y+d*b,p=U*normalize(z),c=2.;"
             "}"
           "float J=-y.z/d.z;"
           "if(J>0.&&J<b)"
             "{"
               "vec3 H=y+d*J;"
               "if(r(H.xy)>.5)"
                 "z=H,b=J,p=v.xxz,a=z.xy,c=3.;"
             "}"
           "if(c==0.)"
             "{"
               "O=vec3(0.);"
               "if(q==0.)"
                 "P=vec3(100.)*pow(max(0.,dot(d,u)),300.);"
               "else"
                 " if(q==1.)"
                   "{"
                     "vec2 H=d.xz*(10.-y.y)/d.y*.1;"
                     "float I=n(H)*.5+n(H*12.1)*.125+n(H*3.8)*.125+n(H*9.)*.0625;"
                     "P=.3*vec3(.3,.5,.9)+400.*vec3(.9,.6,.2)*pow(max(0.,dot(d,u)),400.)+vec3(smoothstep(.4,.6,I));"
                   "}"
                 "else"
                   " if(q==2.)"
                     "{"
                       "P=10.*vec3(step(abs(dot(d,u)+mod(t/16.,1.)-1.),.05));"
                       "vec3 I=normalize(vec3(-1.,.5,.1));"
                       "P+=.2*vec3(.4,.1,.3)*pow(max(0.,dot(d,I)),3.);"
                       "I=normalize(vec3(1.,.3,.3));"
                       "P+=.2*vec3(.1,.5,.4)*pow(max(0.,dot(d,I)),4.);"
                     "}"
             "}"
           "else"
             " if(c==2.)"
               "{"
                 "if(F==0.)"
                   ";"
                 "else"
                   " if(F==1.)"
                     "O=vec3(.8),M=vec2(.5,mix(.95,.8,S)),N=0.;"
                   "else"
                     " if(F==2.)"
                       "O=vec3(1.),M=vec2(1.,1.),N=.2;"
                     "else"
                       " if(F==3.)"
                         "P=p,O=vec3(1.),M=vec2(1.,1.),N=.4;"
                       "else"
                         " if(F==4.)"
                           "O=vec3(1.),M=vec2(0.,mix(.85,.8,S)),N=.4;"
               "}"
             "else"
               " if(c==1.)"
                 "{"
                   "if(D==0.)"
                     ";"
                   "else"
                     " if(D==1.)"
                       "{"
                         "O=vec3(.9);"
                         "N=.5;"
                         "a*=4.;"
                         "vec2 I=floor(a),H=fract(a)*2.-1.,G=vec2(floor(t/2.)+length(I),atan(I.x,I.y));"
                         "P+=vec3(1.,.6,.3)*step(.97,s(G));"
                         "P+=vec3(.3,.4,.9)*step(.98,s(G+1.));"
                         "P+=vec3(.3,.8,.3)*step(.98,s(G+2.));"
                         "P*=vec3(step(length(H),.8),step(length(H+vec2(.1,0.)),.8),step(length(H+vec2(0.,.1)),.8));"
                         "N=.001+s(I)*.3;"
                         "O*=step(abs(H.x),.9)*step(abs(H.y),.9);"
                         "P*=10.;"
                       "}"
                     "else"
                       " if(D==2.)"
                         "{"
                           "O=vec3(1.);"
                           "float I=.5*n(a)+.25*n(a*1.9)+.125*n(a*3.9)+.0625*n(a*9.);"
                           "N=.2+.2*I;"
                         "}"
                       "else"
                         " if(D==3.)"
                           "{"
                             "O=vec3(1.);"
                             "a*=.7;"
                             "float I=.5*n(a)+.25*n(a*1.9)+.125*n(a*3.9)+.0625*n(a*9.);"
                             "N=.01+.2*smoothstep(.4,.5,I);"
                           "}"
                         "else"
                           " if(D==4.)"
                             "O=vec3(1.),N=.3,P=vec3(step(abs(z.z+2.*floor(mod(t,8.))-8.),1.));"
                 "}"
               "else"
                 " P=vec3(4.),O=vec3(0.);"
           "e+=R*P;"
           "R*=O;"
           "if(all(lessThan(R,vec3(.001))))"
             "break;"
           "if(f(l+=z.y)>M.x)"
             "y=z-.01*p,d=normalize(refract(d,p,M.y)),U=-U;"
           "else"
             " y=z+.01*p,d=normalize(mix(reflect(d,p),vec3(f(l+=z.z),f(l+=d.x),f(l+=z.y))-.5,N)),d*=sign(dot(d,p));"
         "}"
     "}"
   "gl_FragColor=vec4(sqrt(e/32.),0.);"
 "}";

#endif // SHADER_H_
