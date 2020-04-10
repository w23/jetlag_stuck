#version 140
uniform vec2 R;
uniform float t;
uniform sampler2D text;

const float PI=3.1415923;
const vec3 E=vec3(0.,.01,1.);

// Is not random enough
/* uint rand_seed = 1u; */
/* float rand() { */
/* 	rand_seed = rand_seed * 1103515245u + 12345u; */
/* 	//rand_seed = (rand_seed * 1664525u + 1013904223u); */
/* 	//return float(rand_seed&0xffffu) / float(0xffffu); */
/* 	//return float(rand_seed>>16) / float(0xffffu); */
/* 	return float(rand_seed) / float(0xffffffffu); */
/* } */

float hash1(float f){return fract(sin(f)*46347.423874);}
float hash2(vec2 v){return hash1(dot(v,vec2(79.53248,31.4328)));}

float noise2(vec2 v) {
	vec2 V=floor(v);v-=V;
	v*=v*(3.-2.*v);
	return mix(
		mix(hash2(V+E.xx), hash2(V+E.zx), v.x),
		mix(hash2(V+E.xz), hash2(V+E.zz), v.x), v.y);
}

/* mat3 RX(float a){float c=cos(a),s=sin(a);return mat3(1.,0.,0.,0.,c,s,0.,-s,c);} */
/* mat3 RY(float a){float c=cos(a),s=sin(a);return mat3(c,0.,s,0.,1.,0.,-s,0.,c);} */
/* mat3 RZ(float a){float c=cos(a),s=sin(a);return mat3(c,s,0.,-s,c,0.,0.,0.,1.);} */

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

// FIXME: replace by simpler wavelength <-> sRGB function
vec3 hsv2rgb(vec3 c) {
    vec4 K = vec4(1.0, 2.0 / 3.0, 1.0 / 3.0, 3.0);
    vec3 p = abs(fract(c.xxx + K.xyz) * 6.0 - K.www);
    return c.z * mix(K.xxx, clamp(p - K.xxx, 0.0, 1.0), c.y);
}
bool mask(vec2 p) {
	return texture2D(text, p * .01).r > .5;
}

float w(vec3 p) { return min(p.y+2., length(p)-2.); }
vec3 wn(vec3 p) { return normalize(vec3(
	w(p+E.yxx),w(p+E.xyx),w(p+E.xxy))-w(p));
}

void main() {
	vec2 uv=(gl_FragCoord.xy/R)*2.-1.;uv.x*=R.x/R.y;

	//gl_FragColor = texture2D(text, uv); return;

	/* rand_seed = 0xfffffffu*uint(hash2(uv)+hash1(t)); */
	/* rand_seed = uint(gl_FragCoord.x + gl_FragCoord.y); */
	/* rand_seed += uint(t); */
	/* rand(); */
	/* rand_seed *= rand_seed; */
	//gl_FragColor=vec4(rand(), rand(), rand(), rand());return;
	//gl_FragColor=vec4(noise2(uv*10.));return;
	//gl_FragColor=vec4(printText(uv*80.));return;

	vec3 c=vec3(0.);
	float seed = t;

	//mat3 mv = mat3(1,0,0,0,1,0,0,0,1);

	vec3 sundir = normalize(vec3(-.5,.5,.4));

	int phase = int(t/32.);

	// PATHTRACER STARTS
	float NS = 16.;
	//float T = t;
	sundir = normalize(sundir + vec3(0.,.5*sin(t/4.),0.));
	for (float s=0.;s<NS;++s) {
		float tt = t;
		float yu = 0.;
		float dof = .4;
		float fov = 1. + 1. * fract(t/32.);

		vec3 ca = vec3(0., -1., 0.);
		vec3 cp = vec3(0., 1.9, 5.);
		cp.z += 4. * fract(t/32.);
		float lfoc = length(ca-cp);
		//lfoc = 1. + 24. * fract(t / 64.);
		if (phase < 6) {
			//ca = vec3(0.);
			//cp = vec3(0., 0., 5.);
			//fov = 2.;
		}

		float a = hash1(seed+=uv.x)*6.2831;
		O = normalize(ca-cp);
		D = normalize(cross(O, vec3(.3*yu,1.,0.)));
		vec3 up = normalize(cross(D, O));
		mat3 mv = mat3(D, up, -O);

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
		if (phase < 6) {
			//O = vec3(0.);
			D = mv*normalize(vec3(uv/fov, -1.)*lfoc - O);
			O = mv*O + cp;

			float d, l = 0., L = 5. + 10. * float(phase);
			float steps = 10. + 10. * float(phase);
			//float steps = 1. + mod(t, 10.);
			for (float i = 0.; i < steps; ++i) {
				d = w(P = O + D * l);
				l += d;
				if (d < .001 || l > L) break;
			}
			if (l < L) {
				if (phase == 0) c = vec3(1.);
				else if (phase == 1) c = vec3(l/L);
				else if (phase == 2) c = P;
				else if (phase == 3) c = fract(P);
				else if (phase == 4) c = wn(P);
				else c = vec3(max(0., dot(wn(P), sundir)));
				c *= NS;
			}
			break;
		}

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
		float mballmat = mod(floor(tt/4.)/*TODO beat sync, 4th beat is earlier*/,4.);
		float mfloormat = mod(floor(tt/6.),3.);

		//mskymat = 1.; mballmat = 0.; mfloormat = 3.;

		for (int i = 0; i < 6; ++i) {
			vec3 me = vec3(0.), ma = vec3(0.);
			float mr = 0.;
			vec2 mf = vec2(1., 1.);
			l = 1e6;
			M = 0;

			if (D.y < 0.) {
				l = (-2. - O.y) / D.y;
				P = O + D * l;
				N = E.xzx;
				UV = P.xz;
				M = 1;
			}

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
					//vec3(rand(), rand(), rand())*2.-1., mr));
				D *= sign(dot(D, N));
			}
		}
	}

	c /= NS;

	gl_FragColor=vec4(sqrt(c), 0);
}
