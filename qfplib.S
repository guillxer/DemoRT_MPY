@ Copyright 2015 Mark Owen
@ http://www.quinapalus.com
@ E-mail: qfp@quinapalus.com
@
@ This file is free software: you can redistribute it and/or modify
@ it under the terms of version 2 of the GNU General Public License
@ as published by the Free Software Foundation.
@
@ This file is distributed in the hope that it will be useful,
@ but WITHOUT ANY WARRANTY; without even the implied warranty of
@ MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
@ GNU General Public License for more details.
@
@ You should have received a copy of the GNU General Public License
@ along with this file.  If not, see <http://www.gnu.org/licenses/> or
@ write to the Free Software Foundation, Inc., 51 Franklin Street,
@ Fifth Floor, Boston, MA  02110-1301, USA.


.syntax unified
.text
.thumb
.cpu cortex-m0

@ exported symbols

.global qfp_fadd
.global qfp_fsub
.global qfp_fmul
.global qfp_fcmp
.global qfp_float2int
.global qfp_float2fix
.global qfp_float2uint
.global qfp_float2ufix
.global qfp_int2float
.global qfp_fix2float
.global qfp_uint2float
.global qfp_ufix2float


@ exchange r0<->r1, r2<->r3
xchxy:
 push {r0,r2,r14}
 mov r0,r1
 mov r2,r3
 pop {r1,r3,r15}

@ IEEE single precision floats in r0,r1-> mantissae in r1,r0 exponents in r3,r2 *respectively*
@ trashes r4
unpackxy:
 push {r14}
 bl unpackx
 bl xchxy
 pop {r4}
 mov r14,r4

@ IEEE single in r0-> signed (two's complemennt) mantissa in r0 9Q23 (24 significant bits), signed exponent (bias removed) in r2
@ trashes r4; zero, denormal -> mantissa=+/-1, exponent=-380; Inf, NaN -> mantissa=+/-1, exponent=+640
unpackx:
 lsrs r2,r0,#23 @ save exponent and sign
 lsls r0,#9     @ extract mantissa
 lsrs r0,#9
 movs r4,#1
 lsls r4,#23
 orrs r0,r4     @ reinstate implied leading 1
 cmp r2,#255    @ test sign bit
 uxtb r2,r2     @ clear it
 bls 1f         @ branch on positive
 rsbs r0,#0     @ negate mantissa
1:
 subs r2,#1
 cmp r2,#254    @ zero/denormal/Inf/NaN?
 bhs 2f
 subs r2,#126   @ remove exponent bias: can now be -126..+127
 bx r14

2:              @ here with special-case values
 cmp r0,#0
 mov r0,r4      @ set mantissa to +1
 bpl 3f
 rsbs r0,#0     @ zero/denormal/Inf/NaN: mantissa=+/-1
3:
 subs r2,#126   @ zero/denormal: exponent -> -127; Inf, NaN: exponent -> 128
 lsls r2,#2     @ zero/denormal: exponent -> -508; Inf, NaN: exponent -> 512
 adds r2,#128   @ zero/denormal: exponent -> -380; Inf, NaN: exponent -> 640
 bx r14

@ normalise and pack signed mantissa in r0 nominally 3Q29, signed exponent in r2-> IEEE single in r0
@ trashes r4, preserves r1,r3
@ r5: "sticky bits", must be zero iff all result bits below r0 are zero for correct rounding
packx:
 lsrs r4,r0,#31 @ save sign bit
 lsls r4,r4,#31 @ sign now in b31
 bpl 2f         @ skip if positive
 cmp r5,#0
 beq 11f
 adds r0,#1     @ fiddle carry in to following rsb if sticky bits are non-zero
11:
 rsbs r0,#0     @ can now treat r0 as unsigned
packx0:
 bmi 3f         @ catch r0=0x80000000 case
2:
 subs r2,#1     @ normalisation loop
 adds r0,r0
 beq 1f         @ zero? special case
 bpl 2b         @ normalise so leading "1" in bit 31
3:
 adds r2,#129   @ (mis-)offset exponent
 bne 12f        @ special case: highest denormal can round to lowest normal
 adds r0,#0x80  @ in special case, need to add 256 to r0 for rounding
 bcs 4f         @ tripped carry? then have leading 1 in C as required
12:
 adds r0,#0x80  @ rounding
 bcs 4f         @ tripped carry? then have leading 1 in C as required (and result is even so can ignore sticky bits)
 cmp r5,#0
 beq 7f         @ sticky bits zero?
8:
 lsls r0,#1     @ remove leading 1
9:
 subs r2,#1     @ compensate exponent on this path
4:
 cmp r2,#254
 bge 5f         @ overflow?
 adds r2,#1     @ correct exponent offset
 ble 10f        @ denormal/underflow?
 lsrs r0,#9     @ align mantissa
 lsls r2,#23    @ align exponent
 orrs r0,r2     @ assemble exponent and mantissa
6:
 orrs r0,r4     @ apply sign
1:
 bx r14

5:
 movs r0,#0xff  @ create infinity
 lsls r0,#23
 b 6b

10:
 movs r0,#0     @ create zero
 bx r14

7:              @ sticky bit rounding case
 lsls r5,r0,#24 @ check bottom 8 bits of r0
 bne 8b         @ in rounding-tie case?
 lsrs r0,#9     @ ensure even result
 lsls r0,#10
 b 9b

@ unpack two arguments (r0,r1) and shift one down to have common exponent, returned in r2; note that arguments are exchanged
@ sticky bits shifted off bottom of smaller argument in r5
@ following code is unnecessarily general for fadd, but is shared with atan2
unpackxyalign:
 push {r14}
 bl unpackxy
 lsls r0,r0,#6  @ Q29
 lsls r1,r1,#6  @ Q29
 subs r4,r2,r3  @ calculate shift
 bge 1f         @ x>=y?
 mov r2,r3      @ no: take common exponent from y
 mov r5,r0      @ potential sticky bits from x
 rsbs r4,#0     @ make shift positive
 asrs r0,r4
 cmp r4,#32
 blo 2f
 movs r0,#0    @ large shift, so all bits are sticky and result is zero
 pop {r15}
1:
 mov r5,r1     @ common exponent from x; potential sticky bits from y
 asrs r1,r4
 cmp r4,#32
 blo 2f
 movs r1,#0    @ large shift, so all bits are sticky and result is zero
 pop {r15}
2:
 rsbs r4,#0
 adds r4,#32
 lsls r5,r4    @ extract sticky bits
 pop {r15}

.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_fsub:
 movs r2,#1    @ subtract: flip sign bit of second argument and fall through to fadd
 lsls r2,#31
 eors r1,r2
.thumb_func
qfp_fadd:
 push {r4,r5,r14}
 bl unpackxyalign
 adds r0,r1    @ do addition
 bne 2f        @ not in Inf-Inf case?
 cmp r2,#200
 blt 2f
 movs r0,#1
 lsls r0,#29   @ for Inf-Inf, set mantissa to +1 to prevent zero result
2:
packret:       @ common return point: "pack and return"
 bl packx
 pop {r4,r5,r15}

@ signed multiply r0 1Q23 by r1 4Q23, result in r0 7Q25, sticky bits in r5
@ trashes r3,r4
mul0:
 uxth r3,r0      @ Q23
 asrs r4,r1,#16  @ Q7
 muls r3,r4      @ L*H, Q30 signed
 asrs r4,r0,#16  @ Q7
 uxth r5,r1      @ Q23
 muls r4,r5      @ H*L, Q30 signed
 adds r3,r4      @ sum of middle partial products
 uxth r4,r0
 muls r4,r5      @ L*L, Q46 unsigned
 lsls r5,r4,#16  @ initialise sticky bits from low half of low partial product
 lsrs r4,#16     @ Q25
 adds r3,r4      @ add high half of low partial product to sum of middle partial products
                 @ (cannot generate carry by limits on input arguments)
 asrs r0,#16     @ Q7
 asrs r1,#16     @ Q7
 muls r0,r1      @ H*H, Q14 signed
 lsls r0,#11     @ high partial product Q25
 lsls r1,r3,#27  @ sticky
 orrs r5,r1      @ collect further sticky bits
 asrs r1,r3,#5   @ middle partial products Q25
 adds r0,r1      @ final result
 bx r14

.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_fcmp:
 movs r2,#1      @ initialise result
 lsls r3,r2,#31  @ r3=0x80000000
 tst r0,r3       @ check sign of first argument
 beq 1f
 subs r0,r3,r0   @ convert to 2's complement form for direct comparison
1:
 tst r1,r3       @ repeat for second argument
 beq 2f
 subs r1,r3,r1
2:
 subs r0,r1     @ perform comparison
 beq 4f         @ equal? return 0
 bgt 3f         @ r0>r1? return +1
 rsbs r2,#0     @ r0<r1: return -1
3:
 mov r0,r2
4:
 bx r14

.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_fmul:
 push {r4,r5,r14}
 bl unpackxy
 add r2,r3      @ sum exponents
 adds r2,#4     @ adjust exponent for pack which expects Q29
fmul0:
 bl mul0
 b packret


@ convert float to signed int, rounding towards -Inf, clamping
.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_float2int:
 movs r1,#0      @ fall through

@ convert float in r0 to signed fixed point in r0, clamping
.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_float2fix:
 push {r4,r14}
 bl unpackx
 add r2,r1       @ incorporate binary point position into exponent
 subs r2,#23     @ r2 is now amount of left shift required
 blt 1f          @ requires right shift?
 cmp r2,#7       @ overflow?
 ble 4f
3:               @ overflow
 asrs r1,r0,#31  @ +ve:0 -ve:0xffffffff
 mvns r1,r1      @ +ve:0xffffffff -ve:0
 movs r0,#1
 lsls r0,#31
5:
 eors r0,r1      @ +ve:0x7fffffff -ve:0x80000000 (unsigned path: 0xffffffff)
 pop {r4,r15}
1:
 rsbs r2,#0      @ right shift for r0, >0
 cmp r2,#32
 blt 2f          @ more than 32 bits of right shift?
 movs r2,#32
2:
 asrs r0,r0,r2
 pop {r4,r15}

@ unsigned version
.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_float2uint:
 movs r1,#0      @ fall through
.thumb_func
qfp_float2ufix:
 push {r4,r14}
 bl unpackx
 add r2,r1       @ incorporate binary point position into exponent
 movs r1,r0
 bmi 5b          @ negative? return zero
 subs r2,#23     @ r2 is now amount of left shift required
 blt 1b          @ requires right shift?
 mvns r1,r0      @ ready to return 0xffffffff
 cmp r2,#8       @ overflow?
 bgt 5b
4:
 lsls r0,r0,r2   @ result fits, left shifted
 pop {r4,r15}

@ convert signed int to float, rounding
.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_int2float:
 movs r1,#0      @ fall through

@ convert signed fix to float, rounding; number of r0 bits after point in r1
.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_fix2float:
 push {r4,r5,r14}
1:
 movs r2,#29
 subs r2,r1      @ fix exponent
packretns:       @ pack and return, sticky bits=0
 movs r5,#0
 b packret

@ unsigned version
.syntax unified
.text
.thumb
.cpu cortex-m0
.thumb_func
qfp_uint2float:
 movs r1,#0      @ fall through
.thumb_func
qfp_ufix2float:
 push {r4,r5,r14}
 cmp r0,#0
 bge 1b          @ treat <2^31 as signed
 movs r2,#30
 subs r2,r1      @ fix exponent
 lsls r5,r0,#31  @ one sticky bit
 lsrs r0,#1
 b packret


qfp_lib_end:
