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

vec3 hsv2rgb(vec3 c) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}

const float C_A = 434073., C_B = 497559., C_C = 397590., C_D = 498071.,C_E = 988959., C_F = 988945., C_G = 400790., C_H = 630681.,C_I = 467495., C_J = 467491., C_K = 611161., C_L = 69919.,C_M = 653721., C_N = 638361., C_O = 432534., C_P = 497425.,C_Q = 432606., C_R = 497497., C_S = 923271., C_T = 991778.,C_U = 629142., C_V = 629075., C_W = 646615., C_X = 628377.,C_Y = 628292., C_Z = 1016879., C_1 = 291919., C_2 = 493087.,C_3 = 495239., C_4 = 630408., C_5 = 988807., C_6 = 272278.,C_7 = 1016900., C_8 = 431766., C_9 = 433730., C_0 = 433590.,C_dot = 1024.;

float gB(float g, vec2 gp){
	return (gp.x<4.&&gp.y<5.) ? mod(floor(g / pow(2., gp.y*4. + gp.x)), 2.) : 0.;
}

#define PC(g) if(pc.x==lx){col=gB(g,pg);}lx+=1.

float diGlyph(in float di) {
	if (di == 0.) return C_0;
	if (di == 1.) return C_1;
	if (di == 2.) return C_2;
	if (di == 3.) return C_3;
	if (di == 4.) return C_4;
	if (di == 5.) return C_5;
	if (di == 6.) return C_6;
	if (di == 7.) return C_7;
	if (di == 8.) return C_8;
	if (di == 9.) return C_9;
	return C_E;
}
void printInt(in float num, in vec2 pg, in vec2 pc, inout float lx, inout float col) {
	/*if (num < 0.) {PC(C_N);
	num *= -1.;
	} else {PC(diGlyph(mod(floor(num/1000.),10.)));
}*/

	if (num >= 1000.) { PC(diGlyph(mod(floor(num/1000.),10.))); }
	if (num >= 100.) { PC(diGlyph(mod(floor(num/100.),10.))); }
	if (num >= 10.) { PC(diGlyph(mod(floor(num/10.),10.))); }
	PC(diGlyph(mod(floor(num),10.)));
}

float printText(vec2 p) {
#define PIXSZ 2.
	p = floor(p / PIXSZ);
	vec2 pc = floor(p / vec2(5.,6.));
	vec2 pg = mod(p, vec2(5.,6.));
	float lx = 1.;
	float col = 0.;

#define PUTN(n) printInt(n,pg,pc,lx,col)
	//float rnd = floor(hash1(floor(t/8.)) * 10.);
	float rnd = mod(floor(t/8.), 11.);
	if (pc.y == 0.) {
	PUTN(t);
	PC(0.);
	PC(0.);
	PC(0.);
	//PC(C_L); PC(C_E); PC(C_V); PC(C_E); PC(C_L); PC(0.);
	if (rnd == 0.) { PC(C_L);PC(C_O);PC(C_G);PC(C_I);PC(C_C);PC(C_O);PC(C_M);PC(C_A); }
	else if (rnd == 1.) { PC(C_F);PC(C_A);PC(C_R);PC(C_B);PC(C_R);PC(C_A);PC(C_U);PC(C_S);PC(C_C);PC(C_H); }
	else if (rnd == 2.) { PC(C_L);PC(C_J); }
	else if (rnd == 3.) { PC(C_P);PC(C_R);PC(C_I);PC(C_S);PC(C_M);PC(C_B);PC(C_E);PC(C_I);PC(C_N);PC(C_G);PC(C_S); }
	else if (rnd == 4.) { PC(C_A);PC(C_L);PC(C_C);PC(C_A);PC(C_T);PC(C_R);PC(C_A);PC(C_Z); }
	else if (rnd == 5.) { PC(C_C);PC(C_O);PC(C_N);PC(C_S);PC(C_P);PC(C_I);PC(C_R);PC(C_A);PC(C_C);PC(C_Y); }
	else if (rnd == 6.) { PC(C_Q);PC(C_U);PC(C_I);PC(C_T);PC(C_E); }
	else if (rnd == 7.) { PC(C_M);PC(C_E);PC(C_R);PC(C_C);PC(C_R);PC(C_U);PC(C_R);PC(C_Y); }
	else if (rnd == 8.) { PC(C_S);PC(C_A);PC(C_N);PC(C_D);PC(C_S); }
	else if (rnd == 9.) { PC(C_T);PC(C_I);PC(C_T);PC(C_A);PC(C_N); }
	else if (rnd == 10.) { PC(C_T);PC(C_H);PC(C_R);PC(C_O);PC(C_B); }
	} else if (pc.y == 1.) {
	}

	return col;
}

bool mask(vec2 p) {
	p.x += 32.;
	float g = printText(p);
	g *= step(hash2(floor((p+floor(vec2(.4,.7)*t))/PIXSZ)), .8);
	return g > 0.;
}

void main() {
	vec2 uv=(gl_FragCoord.xy/R)*2.-1.;uv.x*=R.x/R.y;

	//gl_FragColor=vec4(noise2(uv*10.));return;
	//gl_FragColor=vec4(printText(uv*80.));return;

	vec3 c=vec3(0.);
	float seed = t;

	////////////// TWEAK THESE
	float ls = .05, // lens size
				lf = 3.,// * sin(t*.1), // focus distance
				lS = 1.; // fov-ish
	//vec3 sz = vec3(3.5,2.,3.);
	float or = 8.;
	//mat3 mv = mat3(1,0,0,0,1,0,0,0,1);
	//mv = RY(0.);
	//mv = RY(2.*sin(t*.01))*RY(3.*sin(t*.039));

	vec3 sundir = normalize(vec3(-.5,.5,.4));

	// PATHTRACER STARTS
	lS *= lf;
	float NS = 64.;
	//float T = t;
	sundir = normalize(sundir + vec3(0.,.5*sin(t/4.),0.));
	for (float s=0.;s<NS;++s) {

		/*
		//float tt = t - s/NS/60.; // FIXME merge-with-previous-frame blur
		float tt = 221.;
		mv = RX(.2+.6*cos(tt*.037))*RY(2.+2.*sin(tt*.05));
		float a = hash1(seed+=uv.y)*2.*PI,
					r = s/NS;//*hash1(seed+=uv.y);
		O = vec3(vec2(cos(a),sin(a))*r*ls, or);
		vec3 at = vec3(uv*lS, O.z-lf);
		D = normalize(at - O) * mv; O *= mv;
		//O.z = 1.;
		O += vec3(5., 0., 0.);
		*/

		float tt = t;

		vec3 ca = vec3(0., -1., 0.);
		vec3 cp = vec3(0., 1.9, 5.);

		float yu = 0.;
		float dof = .4;
		float fov = 1. + 1. * fract(t/32.); cp.z += 4. * fract(t/32.);
		float lfoc = length(ca-cp);
		//lfoc = 1. + 24. * fract(t / 64.);

	float a = hash1(seed+=uv.x)*6.2831;
	O = normalize(ca-cp);
	D = normalize(cross(O, vec3(.3*yu,1.,0.)));
	vec3 up = normalize(cross(D, O));
	mat3 mv = mat3(D, up, -O);

	// TODO dolly zoom

	// TODO AA
	vec2 uvaa = vec2(0.);

	// TODO different aprerture forms
	O = vec3(vec2(cos(a),sin(a))*sqrt(hash1(seed+=uv.y))*dof,0.);
	D = mv*normalize(vec3((uv + uvaa) / fov, -1.)*lfoc - O);
	O = mv*O + cp;

		vec3 kc = vec3(1.);
		float ins = 1.;
		float hue = hash1(seed += P.x);
		kc = hsv2rgb(vec3(hue,1.,1.));
		float mskymat = mod(floor(tt/8.),4.);
		float mballmat = mod(floor(tt/4.),4.);
		float mfloormat = mod(floor(tt/6.),3.);

		mskymat = 1.;
		mballmat = 0.;
		mfloormat = 3.;

		for (int i = 0; i < 6; ++i) {
			vec3 me = vec3(0.), ma = vec3(0.);
			float mr = 0.;
			vec2 mf = vec2(1., 1.);
			l = 99.;
			M = 0;

			if (D.y < 0.) {
				l = (-2. - O.y) / D.y;
				P = O + D * l;
				N = E.xzx;
				UV = P.xz;
				M = 1;
			}

			// TODO buildup:
			// 1. sphere on a plane, no lighting, just color, no AAA ("fake" raymarching artifacts)
			// 2. length-based color
			// 3. normal
			// 4. simple lambert lighting
			// 5. ?????
			// 6. path trace for i in 1..N bounce w simple light source (sky? ball?), materials are very simple
			// 7. wet floor material
			// 8. ????? show more materials?
			// 9. first greeting
			// 10. materials caleidoscope and greetings

			// TODO add small sphere as light source
			// TODO transparency texture pattern
			xsph(vec3(0.,0.,0.), 4., 2, ins);

			// glyph plane
			float lp = -O.z / D.z;
			if (lp > 0. && lp < l) {
				vec3 p = O + D * lp;
				if (mask(p.xy*8.+vec2(0.,3.))) {
					P = p;
					l = lp;
					N = E.xxz;
					UV = P.xy;
					M = 4;
				}
			}

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
				} else {
					me = vec3(0.);
				}
			} else if (M == 4) {
				ma = vec3(0.);
				me = vec3(10.);
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
				} else if (mfloormat == 3.) {
					// https://circularlimited.bandcamp.com/track/disruption
					ma = vec3(1.);
					mr = .3;
					me = vec3(step(abs(P.z + 2. * floor(mod(t, 8.)) - 8.), 1.));
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
