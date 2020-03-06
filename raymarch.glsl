uniform vec2 R;
uniform float t;

const float PI=3.1415923;
const vec3 E=vec3(0.,.01,1.);

float hash1(float f){return fract(sin(f)*46347.423874);}
float hash2(vec2 v){return hash1(dot(v,vec2(79.53248,31.4328)));}

float noise2(vec2 v) {
	vec2 V=floor(v);v-=V;
	v*=v*(3.-2.*v);
	return mix(
		mix(hash2(V+E.xx), hash2(V+E.zx), v.x),
		mix(hash2(V+E.xz), hash2(V+E.zz), v.x), v.y);
}

mat3 RX(float a){float c=cos(a),s=sin(a);return mat3(1.,0.,0.,0.,c,s,0.,-s,c);}
mat3 RY(float a){float c=cos(a),s=sin(a);return mat3(c,0.,s,0.,1.,0.,-s,0.,c);}
mat3 RZ(float a){float c=cos(a),s=sin(a);return mat3(c,s,0.,-s,c,0.,0.,0.,1.);}

vec3 O, D, P, N;
vec2 UV;
int M;
float l;

void xsph(vec3 sc, float sr2, int mi, float ss) {
	vec3 v = sc - O;
	float b = dot(v, D);
	float c = dot(v, v) - sr2;
	float det2 = b * b - c;
	if (det2 < 0.) return;
	det2 = sqrt(det2);
	float t1 = b - det2, t2 = b + det2;
	if (t1 < 0.) t1 = 1e6;
	if (t2 < 0.) t2 = 1e6;
	float ls = min(t1, t2);
	if (ls < l) {
		l = ls;
		P = O + D * l;
		N = ss * normalize(P - sc);
		M = mi;
	}
}

vec3 hsv2rgb(vec3 c)
{
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

void main() {
	vec2 uv=(gl_FragCoord.xy/R)*2.-1.;uv.x*=R.x/R.y;

	//gl_FragColor=vec4(noise2(uv*10.));return;

	vec3 c=vec3(0.);
	float seed = t;

	////////////// TWEAK THESE
	float ls = .05, // lens size
				lf = 3.,// * sin(t*.1), // focus distance
				lS = 1.; // fov-ish
	//vec3 sz = vec3(3.5,2.,3.);
	float or = 4.;
	mat3 mv = mat3(1,0,0,0,1,0,0,0,1);
	//mv = RY(0.);
	//mv = RY(2.*sin(t*.01))*RY(3.*sin(t*.039));

	vec3 sundir = normalize(vec3(-.5,.5,.4));

	// PATHTRACER STARTS
	lS *= lf;
	float NS = 32.;
	//float T = t;
	sundir = normalize(sundir + vec3(0.,.5*sin(t/4.),0.));
	for (float s=0.;s<NS;++s) {
		float tt = t - s/NS/60.; // FIXME merge-with-previous-frame blur
		mv = RX(.2+.6*cos(tt*.037))*RY(2.+2.*sin(tt*.05));
		float a = hash1(seed+=uv.y)*2.*PI,
					r = s/NS;//*hash1(seed+=uv.y);
		O = vec3(vec2(cos(a),sin(a))*r*ls, or);
		vec3 at = vec3(uv*lS, O.z-lf);
		D = normalize(at - O) * mv; O *= mv;
		//O.z = 1.;
		vec3 kc = vec3(1.);
		float ins = 1.;
		float hue = hash1(seed += P.x);
		kc = hsv2rgb(vec3(hue,1.,1.));
		float mskymat = mod(floor(tt/8.),2.);
		float mballmat = mod(floor(tt/4.),4.);
		float mfloormat = mod(floor(tt/6.),3.);
		for (int i = 0; i < 6; ++i) {
			vec3 me = vec3(0.), ma = vec3(0.);
			float mr = 0.;
			vec2 mf = vec2(1., 1.);
			l = 30.;
			M = 0;

			if (D.y < 0.) {
				l = (-2. - O.y) / D.y;
				P = O + D * l;
				N = E.xzx;
				UV = P.xz;
				M = 1;
			}

			xsph(vec3(0.,-0.,1.), 4., 2, ins);

			if (M == 0) {
				if (mskymat == 0.) {
					vec2 skp = D.xz*(10.-O.y)/D.y * .1;
					float sk = noise2(skp)*.5 + noise2(skp*2.1)*.25 + noise2(skp*3.8)*.125;
					me = .3 * vec3(.3,.5,.9)
						+ 400.*vec3(.9,.6,.2) * pow(max(0., dot(D,sundir)), 400.)
						+ vec3(smoothstep(.4,.6,sk));
						;
				} else if (mskymat == 1.) {
					me = 10. * vec3(step(abs(dot(D,sundir)+mod(tt/16.,1.)-1.), .05));
					vec3 bd = normalize(vec3(-1., .5, .1));
					me += .1 * vec3(.4, .1, .3) * pow(max(0., dot(D,bd)), 3.);
					bd = normalize(vec3(1., .3, .3));
					me += .1 * vec3(.1, .5, .4) * pow(max(0., dot(D,bd)), 4.);
				}
			} else if (M == 2) {
				if (mballmat == 0.) {
					ma = vec3(.8);
					mf = vec2(.5, mix(.95,.9,hue)); // FIXME fresnel angle dependent?
					mr = .0;
				} else if(mballmat == 1.) {
					ma = vec3(1.);
					mf = vec2(1.,1.);
					mr = .2;
				} else if(mballmat == 2.) {
					me = N;
					ma = vec3(1.);
					mf = vec2(1.,1.);
					mr = .4;
				} else if(mballmat == 3.) {
					ma = vec3(1.);
					mf = vec2(.0, mix(.85,.8,hue)); // FIXME fresnel angle dependent?
					mr = .4;
				}
			} else {
				if (mfloormat == 0.) {
					ma = vec3(.9);
					mr = .5;
					UV *= 4.;
					vec2 cid = floor(UV);
					vec2 cc = fract(UV)*2. - 1.;

					vec2 rcid = vec2(floor(tt/2.) + length(cid), atan(cid.x, cid.y));

					me += vec3(1.,.6,.3) * step(.97,hash2(rcid));//+hash2(tpat));
					me += vec3(.3,.4,.9) * step(.98,hash2(rcid+1.));//+hash2(tpat));
					me += vec3(.3,.8,.3) * step(.98,hash2(rcid+2.));//+hash2(tpat));

					//me += vec3(1.) * step(mod(length(cid)-t*.5, 14.), 1.);

					me *= vec3(
						step(length(cc),.8),
						step(length(cc+vec2(.1,0.)),.8),
						step(length(cc+vec2(.0,.1)),.8)
					);

					//me = vec3(0.);

					//mr = .002 + step(.5,hash2(cid))*.2;// + .2 * smoothstep(.3, .6,noise2(UV*4.));
					mr = .001 + hash2(cid)*.3;// + .2 * smoothstep(.3, .6,noise2(UV*4.));
					ma *= step(abs(cc.x),.9)*step(abs(cc.y),.9);

					me *= 10.;
					//me = vec3(0.);
					//mr = 1.;
					//ma = vec3(.8);
				} else if (mfloormat == 1.) {
					ma = vec3(1.);
					float ns =
						.5 * noise2(UV)
						+ .25 * noise2(UV * 1.9)
						+ .125 * noise2(UV * 3.9)
						+ .0625 * noise2(UV * 9.);
					mr = .2 + .2 * ns;
				} else if (mfloormat == 2.) {
					ma = vec3(1.);
					UV *= .7;
					float ns =
						.5 * noise2(UV)
						+ .25 * noise2(UV * 1.9)
						+ .125 * noise2(UV * 3.9)
						+ .0625 * noise2(UV * 9.);
					mr = .01 + .2 * smoothstep(.4, .5, ns);
				}
			}

			//kc *= 1. - l/40.;
			c += kc * me;
			kc *= ma;

			if (all(lessThan(kc,vec3(.001)))) break;

			O = P;
			if (/*ins < 0 ||*/ hash1(seed+=P.x) > mf.x) {
				D = normalize(refract(D, N, mf.y));
				O -= .01 * N;
				ins = -ins;
			} else {
				O += .01 * N;
				D = normalize(mix(
					reflect(D, N),
					vec3(hash1(seed+=P.z),hash1(seed+=D.x),hash1(seed+=P.y))*2.-1., mr));
				D *= sign(dot(D, N));
			}
		}
	}

	c /= NS;

	gl_FragColor=vec4(sqrt(c), 0);
}
