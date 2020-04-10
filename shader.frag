//#version 140
//uniform vec2 R;
uniform float t;
uniform sampler2D Tex;

//float PI=3.1415923;
vec3 E=vec3(0.,.01,1.);
vec2 R = vec2(1920., 1080.);

// Is not random enough
/* uint rand_seed = 1u; */
/* float rand() { */
/* 	rand_seed = rand_seed * 1103515245u + 12345u; */
/* 	//rand_seed = (rand_seed * 1664525u + 1013904223u); */
/* 	//return float(rand_seed&0xffffu) / float(0xffffu); */
/* 	//return float(rand_seed>>16) / float(0xffffu); */
/* 	return float(rand_seed) / float(0xffffffffu); */
/* } */

float hash1(float f){return fract(sin(f)*46347.4238);}
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

//float mtext = 0.;
float mask(vec2 p) {
	if (t < 1016.) return 0.;
	if (p.x < 0. || p.x > 32.) return 0.;
	if (p.y > 2.8 ) return 0.;
	//if (p.y > 1.8 || p.y < -1.) return 0.;
	p -= 2.;

		float line = floor(p.y / 1.5);
		//if (line >= 0.) return 0.;
		//if (line > 0. || line < -mtext) return 0.;
		//if (line > 0. || line < -1.) return 0.;
		if (t < 1152.) {
			if (line == -1.) {
			} else if (line == -2.) {
				line = (t < 1056.) ? 0. : -2. - floor(t-1056.)/8.;
			} else {
				return 0.;
			}
		} else {
			line -= 16.;
		}
		//line = 1.;
		//line -= (1. + mod(floor(t)/8., 16.));
		p.y = line * 1.5 + mod(p.y, 1.5);
		//p.y = (line - 1. + line * (1. + mod(floor(t)/8., 16.))) * 1.5 + mod(p.y, 1.5);
		return texture2D(Tex, p/32.).r;
}

float w(vec3 p) { return min(p.y+2., length(p)-2.); }
vec3 wn(vec3 p) { return normalize(vec3(
	w(p+E.yxx),w(p+E.xyx),w(p+E.xxy))-w(p));
}
float tr(vec3 O, vec3 D, float l, float L, float steps) {
	for (float i = 0.; i < steps; ++i) {
		float d = w(O + D * l);
		l += d;
		if (d < .001 || l > L) break;
	}
	return l;
}

void main() {
	vec2 uv=(gl_FragCoord.xy/R)*2.-1.;uv.x*=R.x/R.y;
	//gl_FragColor = texture2D(Tex, uv); return;

	vec3 c=vec3(0.);
	vec3 O, D, P, N;
	vec2 UV;
	float M;
	float l;
	float seed = fract(t + uv.x);
	float bar = floor(t/16.);
	//mat3 mv = mat3(1,0,0,0,1,0,0,0,1);

	vec3 sundir = normalize(vec3(-.5,.5,.4));

	//float T = t;
	//sundir = normalize(sundir + vec3(0.,.5*sin(t/4.),0.));
	sundir = normalize(vec3(1.));
	float dof = .4;
	float fov = 1. + 1. * fract(t/32.);

	vec3 ca = vec3(0., -1., 0.);
	vec3 cp = vec3(0., 1.9, 5.);
	cp.z += 4. * fract(t/32.);
	float lfoc = 0.;

	float mskymat = 0.;
	float mballmat = 0.;
	float mfloormat = 0;
	float yu = 0.;

	if (t < 304.) {
		//mskymat = 0.;
		mskymat = 0.;
		mballmat = 0.;
		mfloormat = 0.;
		dof = .01;
		fov = 2.;

		//ca = vec3(0., -.1, 0);
		ca = vec3(0., 0., 0.);
		cp = vec3(0., 1., 5. + 2. * fract(t/16.));
	} else {
		ca = vec3(
			sin(bar*3.) * 4.,
			sin(bar*4.),
			sin(bar*5.) * 4.);

		mskymat = 1.;
		mballmat = 0.;
		mfloormat = 0.;
		dof = .2;
		fov = 2.;

		//ca = vec3(0., sin(t*.1), 0.);
		cp = vec3(
			sin(bar) * 10.,
			2. + 2. * sin(bar*7.),
			10. * cos(bar*3.));
	}

	if (t > 362.) {
		mskymat = 1.;
		mfloormat = 3.;
		mballmat = 0.;

		float ph = fract((t-362.) / 128.);
		ca = vec3(0., 0., 0.);
		cp = vec3(0., 1., 5. + 32. * ph);
		fov = 1. + 4. * ph;
	}

	if (t > 490.) {
		mballmat = 1.;
		//mfloormat = 1.;
		//mskymat = 3.;
		dof = .4;
	}

	if (t > 580.) {
		ca = vec3(
			sin(bar*3.) * 4.,
			sin(bar*4.),
			sin(bar*5.) * 4.);
		cp = vec3(
			sin(bar) * 10.,
			2. + 2. * sin(bar*7.),
			10. * cos(bar*3.));
	}

	if (t > 704.) {
		mballmat = 2.;
		//mfloormat = 1.;
		mskymat = 2.;
		dof = .4;
	}

	if (t > 832.) {
		// FIXME proper material ranges
		mskymat = mod(floor(t/8.),4.);
		mballmat = mod(floor(t/4.)/*TODO beat sync, 4th beat is earlier*/,5.);
		mfloormat = mod(floor(t/6.),4.);
		if (mskymat == 3.)
			mfloormat = 1.;
		//mtext = 1.;

		float ph = fract(t/16.);
		float a = sin(bar), b = sin(bar+3.);
		/* cp = mix(cp, vec3( */
		/* 	sin(bar+4.) * 10., */
		/* 	2. + 2. * sin(bar*2.), */
		/* 	10. * cos(bar+5.)), ph * ph); */
		cp = ca + mix(
			vec3(cos(a)*10.,2+2.*sin(bar*2.),sin(a)*10.),
			vec3(cos(b)*10.,2+2.*sin(bar*3.),sin(b)*10.), ph * ph);
		cp *= max(1., 8. / length(cp));
		//cp *= min(1., length
		//a = mix(fract(a)*7., 7.*fract(b), ph * ph);
		//cp = vec3(cos(a)*10.,2+2.*sin(bar*2.),sin(a)*10.);
		ca = vec3(
			sin(bar*3.) * 2.,
			sin(bar*4.),
			sin(bar*5.) * 2.);
			yu = sin(bar*17.);
	}

	if (t > 1016.) {
		//mtext = 1.;
		ca = vec3(5., 1., 0.);
		cp.z = abs(cp.z);
		cp *= max(1., 8. / length(cp));

		float ph = fract(t/64.);
		fov = 1. + 2. * ph;
		cp = vec3(8. + 10. * sin(bar + t/16.), 2., 10.);
	}

	if (t > 1152.) {
		//ca = vec3(0.);
		mskymat = 0.;
		cp = vec3(10., 2., 10.);
		yu = 0.;
	}

	if (t > 1184.) {
		mfloormat = 3.;
	}

	if (t > 1216.) {
		mballmat = 1.;
	}

	/*
	if (t > 1232.) {
		mtext = 0.;
	}
	*/

	lfoc += length(ca-cp);

	// PATHTRACER STARTS
	for (float s=0.;s<32.;++s) {
		D = normalize(cross(O = normalize(ca-cp), vec3(yu,1.,0.)));
		vec3 up = normalize(cross(D, O));
		mat3 mv = mat3(D, up, -O);

		if (t < 266.) {
			D = mv*normalize(vec3(uv/fov, -1.)*lfoc - O);
			O = mv*O + cp;

			float L = 100., l = tr(O, D, 0., L, 40.);
			if (l < L) {
				P = O + D * l;
				vec3 n = wn(P);
				if (bar < 2) c = vec3(1.);
				else if (bar < 4) c = vec3(l/L);
				else if (bar < 6) c = P;
				else if (bar < 8) c = fract(P);
				else if (bar < 10) c = n;
				else {
					c += .8 * vec3(max(0., dot(n, sundir)));
					if (bar > 14)
						c *= step(5., tr(P, sundir, .1, 5., 20.));
					if (bar > 12)
						c = vec3(.01) + .5 *c;// + vec3(.5) * pow(max(0., dot(n,normalize(sundir-D))), 100.);
					//c += vec3(.01);
				}
				c *= 32.;
			}
			break;
		}

		// TODO dolly zoom
		//vec2 uvaa = (vec2(hash1(seed+=s), hash1(seed+=s)) - .5)/R;
		vec2 uvaa=vec2(0.);
		// TODO different aprerture forms
		float a = hash1(seed+=O.x)*6.2831;
		O = vec3(vec2(cos(a),sin(a))*sqrt(hash1(seed+=uv.y))*dof,0.);
		D = mv*normalize(vec3((uv + uvaa) / fov, -1.)*lfoc - O);
		O = mv*O + cp;

		float ins = 1.;
		float hue = hash1(seed += P.x);
	// FIXME: replace by simpler wavelength <-> sRGB function
		vec3 kc = clamp(abs(fract(vec3(3.,2.,1.)/3. + hue)*6.-3.)-1.,0.,1.);

		for (float i = 0.; i < 6.; ++i) {
			vec3 me = vec3(0.), ma = vec3(.8);
			float mr = 1.;
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

			// TODO add small sphere as light source?
			// TODO transparency texture pattern
			float b = dot(-O, D);
			float det2 = b * b - dot(-O, -O) + 4.;
			if (det2 >= 0.) {
				det2 = sqrt(det2);
				float t1 = b - det2, t2 = b + det2;
				if (t1 < 0.) t1 = 1e6;
				if (t2 < 0.) t2 = 1e6;
				float ls = min(t1, t2);
				if (ls < l) {
					l = ls;
					P = O + D * l;
					N = ins * normalize(P);
					M = 2;
				}
			}

			// text plane
			float lp = -O.z / D.z;
			if (lp > 0. && lp < l) {
				vec3 p = O + D * lp;
				if (mask(p.xy) > .5) {
					P = p;
					l = lp;
					N = E.xxz;
					UV = P.xy;
					M = 3;
				}
			}

			if (M == 0) { // SKY
				ma = vec3(0.);
				if (mskymat == 0.) {
						me = vec3(100.) * pow(max(0., dot(D,sundir)), 300.);
				} else if (mskymat == 1.) {
					vec2 skp = D.xz*(10.-O.y)/D.y * .1;
					float sk = noise2(skp)*.5 + noise2(skp*12.1)*.125 + noise2(skp*3.8)*.125 + noise2(skp*9.)*.0625;
					//sk =0.;
					me = .3 * vec3(.3,.5,.9)
						+ 400.*vec3(.9,.6,.2) * pow(max(0., dot(D,sundir)), 400.)
						+ vec3(smoothstep(.4,.6,sk));
						;
				} else if (mskymat == 2.) {
					me = 10. * vec3(step(abs(dot(D,sundir)+mod(t/16.,1.)-1.), .05));
					vec3 bd = normalize(vec3(-1., .5, .1));
					me += .2 * vec3(.4, .1, .3) * pow(max(0., dot(D,bd)), 3.);
					bd = normalize(vec3(1., .3, .3));
					me += .2 * vec3(.1, .5, .4) * pow(max(0., dot(D,bd)), 4.);
				}
			} else if (M == 2) { // BALL
				if (mballmat == 0.) {
				} else if (mballmat == 1.) {
					ma = vec3(.8);
					mf = vec2(.5, mix(.95,.8,hue)); // FIXME fresnel angle dependent?
					mr = .0;
				} else if(mballmat == 2.) {
					ma = vec3(1.);
					mf = vec2(1.,1.);
					mr = .2;
				} else if(mballmat == 3.) {
					me = N;
					ma = vec3(1.);
					mf = vec2(1.,1.);
					mr = .4;
				} else if(mballmat == 4.) {
					ma = vec3(1.);
					mf = vec2(.0, mix(.85,.8,hue)); // FIXME fresnel angle dependent?
					mr = .4;
				}
			} else if (M == 1) { // GROUND
				if (mfloormat == 0.) {
					//float r = length(P.xz);
					//me = vec3(10.) * step(mod(r - (t-266), 100.), 1.) * max(0., 1.-(t-266)/16);
				} else if (mfloormat == 1.) {
					ma = vec3(.9);
					mr = .5;
					UV *= 4.;
					vec2 cid = floor(UV);
					vec2 cc = fract(UV)*2. - 1.;

					vec2 rcid = vec2(floor(t/2.) + length(cid), atan(cid.x, cid.y));

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
				} else if (mfloormat == 2.) {
					ma = vec3(1.);
					float ns =
						.5 * noise2(UV)
						+ .25 * noise2(UV * 1.9)
						+ .125 * noise2(UV * 3.9)
						+ .0625 * noise2(UV * 9.);
					mr = .2 + .2 * ns;
				} else if (mfloormat == 3.) {
					ma = vec3(1.);
					UV *= .7;
					float ns =
						.5 * noise2(UV)
						+ .25 * noise2(UV * 1.9)
						+ .125 * noise2(UV * 3.9)
						+ .0625 * noise2(UV * 9.);
					mr = .01 + .2 * smoothstep(.4, .5, ns);
				} else if (mfloormat == 4.) {
					// https://circularlimited.bandcamp.com/track/disruption
					ma = vec3(1.);
					mr = .3;
					me = vec3(step(abs(P.z + 2. * floor(mod(t, 8.)) - 8.), 1.));
				}
			} else { // TEXT
				me = vec3(4.);
				ma = vec3(0.);
			}

			//kc *= 1. - l/40.;
			c += kc * me;
			kc *= ma;

			if (all(lessThan(kc,vec3(.001)))) break;

			if (/*ins < 0 ||*/ hash1(seed+=P.y) > mf.x) {
				O = P - .01 * N;
				D = normalize(refract(D, N, mf.y));
				ins = -ins;
			} else {
				O = P + .01 * N;
				/* if (hash1(seed+=D.z) > mr) { */
				/* 	D = reflect(D, N); */
				/* } else { */
				/* 	D = vec3(hash1(seed+=P.z),hash1(seed+=D.x),hash1(seed+=P.y))-.5; */
				/* } */
				/* D = normalize(D); */
				D = normalize(mix(
					reflect(D, N),
					vec3(hash1(seed+=P.z),hash1(seed+=D.x),hash1(seed+=P.y))-.5, mr));
				D *= sign(dot(D, N));
			}
		}
	}

	gl_FragColor=vec4(sqrt(c/32.), 0.);
}
