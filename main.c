
#include "py/dynruntime.h"

#include "qfplib.h"

#ifndef FIXEDPT_BITS
#define FIXEDPT_BITS	32
#endif

void raytracetest(int lightx, int lighty, int lightz, unsigned char* colorbuffer);

static mp_obj_t productf(mp_obj_fun_bc_t* self, size_t n_args, size_t n_kw, mp_obj_t* args) {
    // Check number of arguments is valid
    //mp_arg_check_num(n_args, n_kw, 4, 4, false);
    int x = mp_obj_get_int(args[0]);
    int y = mp_obj_get_int(args[1]);
    int z = mp_obj_get_int(args[2]);
    // Extract buffer pointer and verify typecode
    mp_buffer_info_t bufinfo;
    mp_get_buffer_raise(args[3], &bufinfo, MP_BUFFER_RW);
    //if (bufinfo.typecode != 'B') {
        //mp_raise_ValueError(MP_ERROR_TEXT("expecting byte array"));
    //}
    // Compute product, store result back in first element of array
    unsigned char* colorddata = bufinfo.buf;
    raytracetest(x, y, z, colorddata);

    return mp_const_none;
}

// This is the entry point and is called when the module is imported
mp_obj_t mpy_init(mp_obj_fun_bc_t* self, size_t n_args, size_t n_kw, mp_obj_t* args) {
    // This must be first, it sets up the globals dict and other things
    MP_DYNRUNTIME_INIT_ENTRY

    // The productf function uses the most general C argument interface
    mp_store_global(MP_QSTR_productf, MP_DYNRUNTIME_MAKE_FUNCTION(productf));

    // This must be last, it restores the globals dict
    MP_DYNRUNTIME_INIT_EXIT
}

//optee_os implementation
unsigned int __aeabi_uidivmod(unsigned n, unsigned p) {
	unsigned i = 1, q = 0;
	if (p == 0) {
		return 0;
	}

	while ((p >> 31) == 0) {
		i = i << 1;		/* count the max division steps */
		p = p << 1;     /* increase p until it has maximum size*/
	}

	while (i > 0) {
		q = q << 1;		/* write bit in q at index (size-1) */
		if (n >= p)
		{
			n -= p;
			q++;
		}
		p = p >> 1; 	/* decrease p */
		i = i >> 1; 	/* decrease remaining size in q */
	}
	return q;
}

inline void *memcpy(void *dest, const void *src, size_t n) {
    for (size_t i = 0; i < n; i++) {
        ((unsigned char*)dest)[i] = ((unsigned char*)src)[i];
    }
	return 0;
}

float __aeabi_i2f(int x) {
	return qfp_int2float(x);
}

typedef int fixedpt;
typedef	long long int fixedptd;

#ifndef FIXEDPT_WBITS
#define FIXEDPT_WBITS	16
#endif

#define FIXEDPT_FBITS	(FIXEDPT_BITS - FIXEDPT_WBITS)

#define fixedpt_add(A,B) ((A) + (B))
#define fixedpt_sub(A,B) ((A) - (B))

#define FIXEDPT_ONE	((fixedpt)((fixedpt)1 << FIXEDPT_FBITS))
#define FIXEDPT_ONE_HALF (FIXEDPT_ONE >> 1)
#define FIXEDPT_TWO	(FIXEDPT_ONE + FIXEDPT_ONE)


fixedpt float_to_fixed(float input)
{
    return (fixedpt)qfp_float2int(qfp_fmul(input, (65536.0f))); //(1 << FIXEDPT_WBITS));
}

/* Multiplies two fixedpt numbers, returns the result. */
inline fixedpt
fixedpt_mul(fixedpt A, fixedpt B)
{
    return (((fixedptd)A * (fixedptd)B) >> FIXEDPT_FBITS);
}

/* Divides two fixedpt numbers, returns the result. */
inline fixedpt
fixedpt_div(fixedpt A, fixedpt B)
{
    return (((fixedptd)A << FIXEDPT_FBITS) / (fixedptd)B);
}

/* Returns the square root of the given number, or -1 in case of error */
fixedpt fixedpt_sqrt(fixedpt A)
{
    int invert = 0;
    int iter = FIXEDPT_FBITS;
    int l, i;

    if (A < 0)
        return (-1);
    if (A == 0 || A == FIXEDPT_ONE)
        return (A);
    if (A < FIXEDPT_ONE && A > 6) {
        invert = 1;
        A = fixedpt_div(FIXEDPT_ONE, A);
    }
    if (A > FIXEDPT_ONE) {
        int s = A;

        iter = 0;
        while (s > 0) {
            s >>= 2;
            iter++;
        }
    }

    /* Newton's iterations */
    l = (A >> 1) + 1;
    for (i = 0; i < iter; i++)
        l = (l + fixedpt_div(A, l)) >> 1;
    if (invert)
        return (fixedpt_div(FIXEDPT_ONE, l));
    return (l);
}

struct fixed1616 {
    int value;
};

inline struct fixed1616 inttofp(int x) {
    struct fixed1616 a;
    a.value = x << FIXEDPT_FBITS;
    return a;
}
struct fixed1616 floattofp(float x) {
    struct fixed1616 a;
    a.value = float_to_fixed(x);
    return a;
}
inline int fptoint(struct fixed1616 a) {
    return a.value >> FIXEDPT_FBITS;
}
inline bool eq(struct fixed1616 a, struct fixed1616 b) {
    return a.value == b.value;
}
inline bool neq(struct fixed1616 a, struct fixed1616 b) {
    return a.value != b.value;
}
inline struct fixed1616 add(struct fixed1616 a, struct fixed1616 b) {
    struct fixed1616 result;
    result.value = fixedpt_add(a.value, b.value);
    return result;
}
struct fixed1616 sub(struct fixed1616 a, struct fixed1616 b) {
    struct fixed1616 result;
    result.value = fixedpt_sub(a.value, b.value);
    return result;
}
struct fixed1616 mulfp(struct fixed1616 a, struct fixed1616 b) {
    struct fixed1616 result;
    result.value = fixedpt_mul(a.value, b.value);
    return result;
}
struct fixed1616 div(struct fixed1616 a, struct fixed1616 b) {
    struct fixed1616 result;
    result.value = fixedpt_div(a.value, b.value);
    return result;
}
inline struct fixed1616 neg(struct fixed1616 a) {
    struct fixed1616 result;
    result.value = -a.value;
    return result;
}
struct fixed1616 sqrtfp(struct fixed1616 a) {
    struct fixed1616 result;
    result.value = fixedpt_sqrt(a.value);
    return result;
}

typedef struct fixed1616 fixed_16_16;

const fixed_16_16 FPZERO 	= {0};
const fixed_16_16 FPONE 	= {1 << FIXEDPT_WBITS};

struct vec3fp
{
    fixed_16_16 x;
    fixed_16_16 y;
    fixed_16_16 z;
};
struct vec3fp vecadd(struct vec3fp a, struct vec3fp b)
{
    struct vec3fp out;
    out.x = add(a.x, b.x);
    out.y = add(a.y, b.y);
    out.z = add(a.z, b.z);
    return out;
}
struct vec3fp vecsub(struct vec3fp a, struct vec3fp b)
{
    struct vec3fp out;
    out.x = sub(a.x, b.x);
    out.y = sub(a.y, b.y);
    out.z = sub(a.z, b.z);
    return out;
}
struct vec3fp vecmul(struct vec3fp a, struct vec3fp b)
{
    struct vec3fp out;
    out.x = mulfp(a.x, b.x);
    out.y = mulfp(a.y, b.y);
    out.z = mulfp(a.z, b.z);
    return out;
}
struct vec3fp vecmuls(struct vec3fp a, fixed_16_16 b)
{
    struct vec3fp out;
    out.x = mulfp(a.x, b);
    out.y = mulfp(a.y, b);
    out.z = mulfp(a.z, b);
    return out;
}
struct vec3fp vecneg(struct vec3fp a)
{
    struct vec3fp out;
    out = vecmuls(a, floattofp(-1.0f));
    return out;
}
 struct vec3fp vecdiv(struct vec3fp a, struct vec3fp b)
{
    struct vec3fp out;
    out.x = div(a.x, b.x);
    out.y = div(a.y, b.y);
    out.z = div(a.z, b.z);
    return out;
}
 struct vec3fp vecdivs(struct vec3fp a, fixed_16_16 b)
{
    struct vec3fp out;
    out.x = div(a.x, b);
    out.y = div(a.y, b);
    out.z = div(a.z, b);
    return out;
}
fixed_16_16 dot(struct vec3fp a, struct vec3fp b)
{
    return add(add(mulfp(a.x, b.x), mulfp(a.y, b.y)), mulfp(a.z, b.z));
}
fixed_16_16 mag(struct vec3fp a)
{
    return sqrtfp(dot(a, a));
}

struct vec3fp vec3(fixed_16_16 inx, fixed_16_16 iny, fixed_16_16 inz) {
    struct vec3fp out;
    out.x = inx;
    out.y = iny;
    out.z = inz;
    return out;
}

struct vec3fp getnormalized(struct vec3fp a)
{
    fixed_16_16 vectormag = mag(a);
    if (eq(vectormag, FPZERO)) {
        return vec3(FPZERO, FPZERO, FPONE);
    }
    else {
        return vecdivs(a, vectormag);
    }
}

struct Ray {
    struct vec3fp origin;
    struct vec3fp direction;
};
 struct Ray MakeRay (struct vec3fp origin, struct vec3fp direction) {
    struct Ray out;
    out.origin = origin;
    out.direction = getnormalized(direction);
    return out;
}
 struct vec3fp point_at_parameter(struct Ray r, fixed_16_16 t) {
    return vecadd(r.origin, vecmuls(r.direction, t));
}

struct Sphere {
    fixed_16_16 radius;
    struct vec3fp center;
};
 struct Sphere MakeSphere(struct vec3fp acenter, fixed_16_16 aradius) {
    struct Sphere s;
    s.center = acenter;
    s.radius = aradius;
    return s;
}
bool intersect(const struct Sphere s, const struct Ray ray, fixed_16_16* t) {
    struct vec3fp oc = vecsub(ray.origin, s.center);
    fixed_16_16 a = dot(ray.direction, ray.direction);
    fixed_16_16 b = mulfp(floattofp(2.0f), dot(oc, ray.direction));
    fixed_16_16 c = sub(dot(oc, oc), mulfp(s.radius, s.radius));
    fixed_16_16 discriminant = sub(mulfp(b, b), mulfp(mulfp(floattofp(4.0f), a), c));
    if (discriminant.value < 0) {
        return false;
    }
    else {
        // t = (-b - sqrt(discriminant)) / ((fixed_16_16)2.0f * a);
        *t = div((sub(neg(b), sqrtfp(discriminant))), (mulfp(floattofp(2.0f), a)));
        return true;
    }
}

struct Plane {
    struct vec3fp position;
    struct vec3fp normal;
};
 struct Plane MakePlane(struct vec3fp aposition, struct vec3fp anormal) {
    struct Plane p;
    p.position = aposition;
    p.normal = anormal;
    return p;
}
bool intersectp(struct Plane p, struct Ray ray, fixed_16_16* t)
{
    fixed_16_16 denom = dot(p.normal, ray.direction);
    if (denom.value > floattofp(0.0001f).value) {
        struct vec3fp p0l0 = vecsub(p.position, ray.origin);
        *t = div(dot(p0l0, p.normal), denom);
        return ((*t).value >= 0);
    }
    return false;
}

 bool visiblecheck(struct Ray ray, const struct Sphere* spheres, int numspheres, fixed_16_16 minlength) {
    bool bvis = true;
    for (int i = 0; i < numspheres; ++i) {
        fixed_16_16 t = FPZERO;
        if (intersect(spheres[i], ray, &t)) {
            if (t.value <= minlength.value) {
                bvis = false;
            }
        }
    }
    return bvis;
}


unsigned char find_closest(int x, int y, fixed_16_16 c0)
{
    const fixed_16_16 dither44[4][4] = {
        { div(floattofp(00.0f), floattofp(16.0f)), 
        div(floattofp(12.0f), floattofp(16.0f)), 
        div(floattofp(03.0f), floattofp(16.0f)), 
        div(floattofp(15.0f), floattofp(16.0f)) },
        { div(floattofp(08.0f), floattofp(16.0f)), 
        div(floattofp(04.0f), floattofp(16.0f)), 
        div(floattofp(11.0f), floattofp(16.0f)), 
        div(floattofp(07.0f), floattofp(16.0f)) },
        { div(floattofp(02.0f), floattofp(16.0f)), 
        div(floattofp(14.0f), floattofp(16.0f)), 
        div(floattofp(01.0f), floattofp(16.0f)), 
        div(floattofp(13.0f), floattofp(16.0f)) },
        { div(floattofp(10.0f), floattofp(16.0f)), 
        div(floattofp(06.0f), floattofp(16.0f)), 
        div(floattofp(09.0f), floattofp(16.0f)), 
        div(floattofp(05.0f), floattofp(16.0f)) } };

    fixed_16_16 limit = floattofp(0.0f);
    if (x < 4)
    {
        limit = (dither44[x][y]);
    }
    if (c0.value <= limit.value)
        return 0;
    return 1;
}

void dithershadeintodisplay(unsigned char* colorbuffer, unsigned char* internalbuffer) {
	const int width = 72;
	const int height = 40;  
    for (int x = 0; x < width; ++x) {
        for (int y = 0; y < height; ++y) {
            fixed_16_16 shadeestimate = floattofp((int)internalbuffer[x * height + y]);
            unsigned char color = find_closest(x % 4, y % 4, div(shadeestimate, floattofp(255.0f)));
            colorbuffer[x * height + y] = color;
        }
    }
}

void raytracetest(int lightx, int lighty, int lightz, unsigned char* colorbuffer) {
	const int width = 72;
	const int height = 40;  
	const fixed_16_16 aspectratio = div(floattofp(width), floattofp(height));
	const int starty = height - 1;
    const int endy = 0;
	unsigned char internalbuffer[width * height];
	
	const fixed_16_16 visbias = floattofp(0.001f);
    struct vec3fp lightpos;
    lightpos.x = mulfp(inttofp(lightx), floattofp(0.001f));
    lightpos.y = mulfp(inttofp(lighty), floattofp(0.001f));
    lightpos.z = mulfp(inttofp(lightz), floattofp(0.001f));
    struct Sphere lightsphere = MakeSphere(lightpos, floattofp(0.1f));
    const struct Sphere sphere0 = MakeSphere(vec3(floattofp(0.5f), FPZERO, floattofp(-1.0f)), floattofp(0.5f));
    const struct Sphere sphere1 = MakeSphere(vec3(floattofp(-0.8f), FPZERO, floattofp(-1.5f)), floattofp(0.5f));
    const struct Sphere spheres[] = { lightsphere, sphere0, sphere1 };
    const struct Sphere visspheres0[] = { sphere1 };
    const struct Sphere visspheres1[] = { sphere0 };
    const struct Sphere visspheres2[] = { sphere0, sphere1 };
    const struct Sphere* visspheres[] = { NULL, visspheres0, visspheres1, visspheres2 };
    const int visspheresnum[] = {
        0,
        sizeof(visspheres0) / sizeof(struct Sphere),
        sizeof(visspheres1) / sizeof(struct Sphere),
        sizeof(visspheres2) / sizeof(struct Sphere) };
    const int numspheres = sizeof(spheres) / sizeof(struct Sphere);

    for (int j = starty; j >= endy; --j) {
        for (int i = 0; i < width; ++i) {
            fixed_16_16 u = div(floattofp(i), floattofp(width));
            fixed_16_16 v = div(floattofp(j), floattofp(height));
            struct Ray ray = MakeRay(
                vec3(
                    FPZERO,
                    FPZERO,
                    FPZERO),
                vec3(
                    mulfp((sub(mulfp(floattofp(2.0f), u), floattofp(1.0f))), aspectratio),
                    sub(mulfp(floattofp(2.0f), v), floattofp(1.0f)),
                    floattofp(-1.0f))
            );

            int r = 0;
            fixed_16_16 t = FPZERO;
            fixed_16_16 mint = floattofp(1000.0f);
            for (int spherei = 0; spherei < numspheres; ++spherei) {
                if (intersect(spheres[spherei], ray, &t)) {
                    if (spherei == 0) {
                        if (t.value < mint.value) {
                            mint = t;
                            r = 255;
                        }
                    }
                    else {
                        if (t.value < mint.value) {
                            struct vec3fp point = point_at_parameter(ray, t);
                            struct vec3fp normal = getnormalized(vecsub(point, spheres[spherei].center));
                            struct vec3fp lightvec = vecneg(getnormalized(vecsub(point, lightpos)));
                            fixed_16_16 length = sqrtfp(dot(vecsub(point, lightpos), vecsub(point, lightpos)));
                            struct Ray visray = MakeRay(vecadd(point, vecmuls(normal, visbias)), lightvec);
                            mint = t;
                            if (visiblecheck(visray, visspheres[spherei], visspheresnum[spherei], length)) {
                                fixed_16_16 color = div((add(dot(lightvec, normal), floattofp(1.0f))), floattofp(2.0f));
                                r = fptoint(mulfp(floattofp(255.99f), color));
								
                            }
                        }
                    }
                }
            }
			const struct Plane plane0 = MakePlane(
			vec3(floattofp(0.0f), floattofp(1.0f), floattofp(0.0f)),
			vec3(floattofp(0.0f), floattofp(1.0f), floattofp(0.0f)));
			struct Plane planes[] = { plane0 };
			const int numplanes = sizeof(planes) / sizeof(struct Plane);
            for (int planei = 0; planei < numplanes; ++planei) {
                if (intersectp(planes[planei], ray, &t)) {
                    if (t.value < mint.value) {
                        struct vec3fp point = point_at_parameter(ray, t);
                        struct vec3fp normal = planes[planei].normal;
                        struct vec3fp lightvec = getnormalized(vecsub(point, lightpos));
                        fixed_16_16 length = sqrtfp(dot(vecsub(point, lightpos), vecsub(point, lightpos)));
                        struct Ray visray = MakeRay(vecadd(point, vecmuls(normal, visbias)), lightvec);
                        mint = t;
                        if (visiblecheck(visray, visspheres[3], visspheresnum[3], length)) {
                            fixed_16_16 color = div(add(dot(lightvec, normal), floattofp(1.0f)), floattofp(2.0f));
                            r = fptoint(mulfp(floattofp(255.99f), color));
                        }
                    }
                }
            }
			internalbuffer[i * height + j] = r;
        }
    }
	dithershadeintodisplay(colorbuffer, internalbuffer);
}


