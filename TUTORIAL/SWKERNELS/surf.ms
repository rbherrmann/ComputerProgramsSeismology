.ND
.nr PS 12
.nr VS 14
.nr LL 6.0i
.nr PO 1.25i
.ps 12
.vs 14
.ll 6.0i
.po 1.25i
.EQ
delim $$
define lam2 ' LAMBDA sup 2 '
define lamsi ' ( LAMBDA sup 2 ~+~ sigma sup 2 I ) '
define pcpb1 ' { partial c } over { partial beta sub 1 } '
define pcpbn ' { partial c } over { partial beta sub n } '
define pcpa1 ' { partial c } over { partial alpha sub 1 } '
define pcpan ' { partial c } over { partial alpha sub n } '
define pupb1 ' { partial U } over { partial beta sub 1 } '
define pupbn ' { partial U } over { partial beta sub n } '
define pupa1 ' { partial U } over { partial alpha sub 1 } '
define pupan ' { partial U } over { partial alpha sub n } '
define pupq1 ' { partial U } over { partial Q sub { beta sub 1 } sup -1 } '
define pupqn ' { partial U } over { partial Q sub { beta sub n } sup -1 } '
define pgpv ' { partial gamma } over { partial v } '
define pgpb1 ' { partial gamma } over { partial beta sub 1 } '
define pgpb ' { partial gamma } over { partial beta  } '
define pgpbn ' { partial gamma } over { partial beta sub n } '
define pgpq1 ' { partial gamma } over { partial Q sub { beta sub 1 } sup -1 } '
define pgpqn ' { partial gamma } over { partial Q sub { beta sub n } sup -1 } '
define pcpq1 ' { partial c } over { partial Q sub { beta sub 1 } sup -1 } '
define pcpqn ' { partial c } over { partial Q sub { beta sub n } sup -1 } '
define pgpq ' { partial gamma } over { partial Q sub v sup -1 } '
define om2c2 ' omega over { 2 c sub o sup 2 } '
define qb1 ' Q sub { beta sub 1 } sup -1 '
define qbn ' Q sub { beta sub n } sup -1 '
define qan ' Q { sub alpha sub n } sup -1 '
define qa1 ' Q { sub alpha sub 1 } sup -1 '
define a1 ' alpha sub 1 '
define an ' alpha sub n '
define b1 ' beta sub 1 '
define bn ' beta sub n '
define qbqa ' { Q sub beta } over { Q sub alpha } '
define qaiqbi ' { Q sub alpha sup -1 } over { Q sub beta sup -1 } '
define pilnqv ' { 1 over {pi Q sub v } } ln ( { omega } over { omega sub r } ) '
define piln ' { 1 over pi } ln ( { omega } over { omega sub r } ) '
define pilnq1 ' { 1 over {pi Q sub beta sub 1 } } ln ( { omega } over { omega sub r } ) '
define pilnqn ' { 1 over {pi Q sub beta sub n } } ln ( { omega } over { omega sub r } ) '
define piqv ' 1 over { pi Q sub v } '
define omc2 ' omega over { 2 c sub o sup 2 } '
define co ' c sub o '
define uo ' U sub o '
define uoco '   {U sub o over c sub o }  '
define uc '   {U  over c  }  '
define uo2co2 ' U sub o sup 2 over c sub o sup 2 '
define u2c2 ' U sup 2 over c sup 2 '
define ccoco ' left ( { c ^-^ c sub o } over { c sub o } right ) '
define 2gupo ' { 2 gamma U sub o } over { pi omega } '
define dcodt ' { dc sub o } over dt '
define pupv ' { partial U } over { partial v } '
define pupm ' { partial U } over { partial m } '
define puopv ' { partial U sub o } over { partial v } '
define uuo ' U over { U sub o } '
define pcopv ' { partial c sub o } over { partial v } '
define pcopm ' { partial c sub o } over { partial m } '
define pcpm ' { partial c } over { partial m } '
define pcpv ' { partial c } over { partial v } '
define pcpq ' { partial c } over { partial Q sub v sup -1 } '
define pupq ' { partial U } over { partial Q sub v sup -1  } '
define pcpa ' { partial c sub o } over { partial alpha } '
define pcpb ' { partial c sub o } over { partial beta } '
define pcpbc ' { partial c  } over { partial beta } '
define pupb ' { partial U sub o } over { partial beta } '
define pupbc ' { partial U  } over { partial beta } '
define qa ' Q sub alpha sup -1 '
define qb ' Q sub beta sup -1 '
define qv ' Q sub v sup -1 '
define sm ' sum left ( pcpb beta qb ~+~ pcpa alpha qa right ) '
define dtp ' { { 2 a } over { r sub i ^+^ r sub i-1 } } '
.EN
.LP
.ce 1
\f3\s+4SURFACE WAVE THEORY\s-4\f1
.LP
This  section discusses the surface-wave partial derivatives
required for the inversion.
.LP
\f3\s+21. Eigenvalue Problem\s-2\f1
.LP
\f2Love Waves\f1
.PP
Following Keilis-Borok \f2et al.\f1 (1989), the 
equation for the SH-wave eigenfunction, $ V sup (3) $, in a cylindrical coordinate system
is 
.sp
.EQ I (1)
d over { dz }  
left [ cpile { { V sup (3) } above { p sub { phi z } }} right ]
^=^
left [ { cpile { 0 above { mu k sup 2 ^-^ rho omega sup 2 }} ~
cpile { { 1 ^/^ mu } above 0 } } right  ] ^
left [ cpile { { V sup (3) } above { p sub { phi z } }} right ]
.EN
.sp
where $ k $ is the wavenumber, $ omega $ is the angular frequency,
$ mu $ is the rigidity and $ rho $ is density.  
The boundary conditions for surface waves are that $ p sub { phi z } (z = 0 ) ^=^0$, and $ V sup (3) (z) ^ -> ^ 0 $, as $ z ^ -> ^ inf $.
.PP
A discrete solution to the boundary conditions satisfies the equation
.sp
.EQ I (2)
L mark ^=^ omega sup 2 I sub 0 ^-^ k sup 2 I sub 1 ^-^ I sub 2
.EN
.sp
.EQ I
lineup  ^=^ omega sup 2 int from 0 to inf rho left ( V sup (3) right ) sup 2 dz
^ - ^ k sup 2 int from 0 to inf  mu left ( V sup (3) right ) sup 2 dz
^-^  int from 0 to inf mu left ( {dV sup (3) } over dz right  ) sup 2 dz
.EN
.EQ I
lineup ^=^ 0
.EN
.sp
where $ L  $ is the Lagrangian. 
The group velocity is determined from the relation
.sp
.EQ I (3)
U ^=^ { d omega } over dk ^=^ { k I sub 1 } over { omega I sub 0 }
.EN
.PP
If the medium is split into constant velocity-density layers, of
thickness $ d sub m $, density $ rho sub m $, and shear-wave
velocity $ beta sub m $, and the layer is bounded by
$ z^=^ z sub m - ^d sub m $ and $ z^=^ z sub m $, then
the partials of the phase velocity with shear velocity for 
fixed frequency and density, and
with density for fixed frequency and shear-wave velocity are
.sp
.EQ I (4)
left ({ partial c } over { partial beta } right ) sub m   ^=^  { beta sub m rho sub m } over
{ U I sub 0 }
 ^ int from { z = z sub m - ^d sub m } to { z sub m }
left [ left ( V sup (3) right ) sup 2 + left ( 1 over k { d V sup (3) } over dz right ) sup 2 right ] ~ dz
.EN
.sp
and
.EQ I (5)
left ({ partial c } over { partial rho }  right )  sub m ^=^  
{ { beta sub m } over { 2 rho sub m } } left ({ partial c } over { partial 
beta } right ) sub m 
^-^
{ c sup 2 } over { 2 U I sub 0 } 
 ^ int from { z = z sub m - ^d sub m } to { z sub m }
 left ( V sup (3) right ) sup 2  dz
.EN
.sp
If a layer boundary depth is perturbed by
an amount $ h $, while keeping the material properties and
the frequencies constant, the variation in phase velocity is
determined from
.sp
.EQ I (6)
{ partial c } over { partial h }
~=~
{ c sup 3 } over { 2 omega sup 2 I sub 1 }
left {
{
omega sup 2 ( V sup (3) ) sup 2 [ rho ] sub +- 
^-^
k sup 2 ( V sup (3) ) sup 2 [ mu ] sub +- 
^+^
[ mu left ({ d V sup (3) } over dz right ) sup 2 ] sub +- 
}
right }
.EN
.sp
where the symbol $ [ ^ ] sub +- $ indicates the computation of the jump in the
parameter across the interface, e.g., 
.sp
.EQ I (7)
 [ ^ chi ^ ] sub +- ^=^ chi ( h ^+^0) ^-^ chi ( h^-^0).
.EN
.sp
.LP
\f2Rayleigh Waves\f1
.PP
In a cylindrical coordinate system, the equation of motion governing the
eigenfunctions is
.sp
.EQ I (8)
{
d over { dz }  
left [ {  
 cpile {{ V sup (1) } above { V sup (2) } above { p sub zz } above { p sub { r z } }
 }}  right ]
^=^
left [ 
{
cpile { 0 above { -k } above { - rho omega sup 2 } above 0 }
~~
cpile { {k lambda ^/^ ( lambda + 2 mu ) } above 0 above 0 above {
- rho omega sup 2 ^+^ 4 k sup 2 mu ( lambda + mu )^/^ ( lambda + 2 mu ) }}
~~
cpile { { 1 ^/^ ( lambda + 2 mu )} above 0 above 0 above { -k lambda ^/^ ( lambda + 2 mu ) }}
~~
cpile { 0 above { 1 ^/^ mu } above k above 0 }
}
 right ]
~
left [ {  
 cpile {{ V sup (1) } above { V sup (2) } above { p sub zz } above { p sub { r z } }
 }}  right ]
}
.EN
.br
.LP
where the vertical component eigenfunction is $ V sup (1) $, and the radial
component eigenfunction is $ V sup (2) $. The ellipticity at the
free surface is $ epsilon ^=^ V sup (2) (0) ^/^ V sup (1) (0) $.
The boundary conditions for surface waves are that 
$ p sub { r z } (z = 0 ) ^=^0$, 
$ p sub { z z } (z = 0 ) ^=^0$, 
$ V sup (1) (z) ^ -> ^ 0 $, as $ z ^ -> ^ inf $.
and $ V sup (2) (z) ^ -> ^ 0 $, as $ z ^ -> ^ inf $.
.PP
A discrete solution to the boundary conditions satisfies the equation
.sp
.EQ I (9)
L ^=^ omega sup 2 I sub 0 ^-^ k sup 2 I sub 1 ^-^ 2 k I sub 2 ^-^ I sub 3 ^=^ 0
.EN
.sp
where $ L  $ is the Lagrangian,
and
.sp
.EQ I (10)
I sub 0 ^=^ int from 0 to inf rho 
left [ {(V sup (1) ) sup 2 ^+^ ( V sup (2) ) sup 2 } right ] dz
.EN
.br
.EQ I (11)
I sub 1 ^=^ int from 0 to inf 
left [ { mu {(V sup (1) ) sup 2 ^+^ 
( lambda + 2 mu )
( V sup (2) ) sup 2 }} right ] dz
.EN
.br
.EQ I (12)
I sub 2 ^=^ int from 0 to inf
left [ {
mu V sup (1) { d V sup (2) } over dz 
^-^
lambda V sup (2) { d V sup (1) } over dz
} right ]
dz
.EN
.br
.EQ I (13)
I sub 3 ^=^ int from 0 to inf
left [ {
( lambda + 2 mu ) left ( { d V sup (1) } over dz right ) sup 2
^+^
(  mu ) left ( { d V sup (2) } over dz right ) sup 2
} right ]
dz
.EN
.LP 
The group velocity is determined from the relation
.sp
.EQ I (14)
U ^=^ { d omega } over dk ^=^ {(^ k I sub 1 ^+^ I sub 2 ^) } ^/^ { omega I sub 0 }
.EN
.PP
If the medium is split into constant velocity-density layers, of
thickness $ d sub m $, density $ rho sub m $,
compressional-wave velocity $ alpha sub m $, and shear-wave
velocity $ beta sub m $, and the layer is bounded by
$ z^=^ z sub m - ^d sub m $ and $ z^=^ z sub m $, then
the partials of the phase velocity with shear velocity for 
fixed frequency, compressional velocity and density,
with compresisonal velocity for fixed frequency, shear velocity and
density, and
with density for fixed frequency and compresional and shear velocity are
.sp
.EQ I (15)
left ( { partial c } over { partial alpha } right ) sub m
^=^ left ( { alpha sub m rho sub m } over { U I sub 0 } right )
 ^ int from { z = z sub m - ^d sub m } to { z sub m }
left [ {
V sup (2) ^-^ 1 over k { d V sup (1) } over dz
 } right ] sup 2 dz
.EN
.br
.EQ I (16)
left ( { partial c } over { partial beta } right ) sub m
^=^ left ( { beta sub m rho sub m } over { U I sub 0 } right )
 ^ int from { z = z sub m - ^d sub m } to { z sub m }
left [ {
left ( { V sup (1) ^+^ 1 over k { d V sup (2) } over dz } right ) sup 2
^+^ 4 over k V sup (2) { d V sup (1) } over dz
} right ]
dz
.EN
.br
.EQ I (17)
left ( { partial c } over { partial rho } right ) sub m
^=^
1 over { 2 rho }
left [ { 
alpha left ({ { partial c } over { partial alpha }} right ) sub m 
^+^ 
beta  left ({ { partial c } over { partial beta  }} right ) sub m
} right ]
^-^ { c sup 2 } over { 2 U I sub 0 }
 ^ int from { z = z sub m - ^d sub m } to { z sub m }
left [ { ( V sup (1) ) sup 2 ^+^ (  V sup (2) ) sup 2 } right ]
dz
.EN
.sp
If a layer boundary depth is perturbed by
an amount $ h $, while keeping the material properties and
the frequencies constant, the variation in phase velocity is
determined from
.sp
.EQ I (18)
left ( {{ partial c } over { partial h } } right )
 ~=~ mark
{ c sup 3 } over { 2 omega ( { omega  I sub 1  ^+^ c I sub 2} ) } ~ .
.EN
.EQ I
lineup ~~~ left {  omega sup 2 [ rho ( (V sup (2) ) sup 2 ^+^ ( V sup (1) ) sup 2 ] sub +-
~ -k sup 2 [ ( mu v sup 1 ) sup 2 ] sub +- ~-~ k sup 2 [ ( lambda ^+^ 2 mu ) ( V sup (2) ) sup 2 ] sub +-
.EN
.EQ I
lineup ~~~~~~  
left "" -  [ ( lambda ^+^ 2 mu ) left ( { d V sup (1) } over dz right ) ] sub +-
~ [ mu left ( { d V sup (2) } over dz  right ) ] sub +- 
~
right }
.EN
.sp
where the symbol $ [ ^ ] sub +- $ indicates the computation of the jump in the
parameter across the interface, e.g., 
.sp
.EQ I (19)
 [ ^ chi ^ ] sub +- ^=^ chi ( h ^+^0) ^-^ chi ( h^-^0).
.EN
Equation (18) differs slightly from that in Keilis-Borok (1989) in the grouping of the
square brackets, in order to implement the correct partials for a water-solid
layer boundary, for which the radial eigenfunction is not continuous.
.LP
\f3\s+22. Numerical Partial Derivative Computation and Causality\s-2\f1 
.LP
The subroutines used in the programs \f2surf96\f1 and \f2joint96\f1  are essentially the same as those in PROGRAMS.330/VOLIII/src of
\f3Computer Programs in Seismology.\f1 Since the inversion programs
requires partial derivatives with respect to layer parameters, a detailed
discussion of the relations used is presented.
.LP
The program \f3srfdis96\f1 computes the phase velocity values for a given mode,
wave type and frequency. When group velocities are computed, the phase velocities
are also computed at the  two periods $ ( 1 ~ +- ~ h ) T $ rather than at the single period
$ T . $ The parameter $ h $, tpically 0.005, is that given on LINE 1 of the control file
\f2sobs.d.\f1 or \f2jobs.d\f1.  Denote the phase velocity output of \fBsrfdis96\f1 $ c sub o . $
.LP
The programs \f3srfdrr96\f1 and \fBsrfdrl96\f1 are essentially the same as
\f3sregn96\f1 and \fBslegn96\f1 except that the front end of the
program rearranges the output for use by the inversion program. Note that the  programs \f3sregn96\f1 and \fBslegn96\f1 provide the phase velocity partial derivatives but not the group velocity partial derivatives. In addition, although these two programs do provide the causal phase velocity and anelastic attenuation coefficient when the medium is described by a causal Q, the partial derivatives are only for infinite Q.
.LP
Given the
partial derivatives of phase velocity with respect
to layer shear velocity, $ { partial c sub o } over { partial beta } , $ compressional
velocity, $  {partial c sub o } over { partial alpha } , $ and group velocity, $ U sub o $, the expressions for group velocity are derived by following Rodi \f2et al\f1 (1975).
.LP
By definition
.EQ I
U ^=^ c over { 1 + { T over c } { dc over dT } } ^= ^ f(c, dc/dT, omega , m)
.EN
where $m$ is a model parameter. By applying the chain rule of differentiation, 
.EQ I
{ partial U } over { partial m } ^=^ { partial U } over { partial c } ^ { partial c } over { partial m } ^+^ { partial U } over { partial left ( { dc over dT } right ) } ^ {{  partial "" } over { partial m } left ( { dc over dT } right ) }  
^=^
{ partial U } over { partial c } ^ { partial c } over { partial m } ^+^ { partial U } over { partial left ({ dc over dT } right ) } ^ { partial "" } over { partial T }{ partial c } over { partial m } 
.EN
where the interchange of the order of partial differentiation is permitted since the $ c ( omega , m ) $ is a continuous function.
After some simple algebra, one obtains the followiong expression:
.sp
.EQ I (20)
pupm mark ~=~ uc ^ left ( { 2 - uc } right ) pcpm + u2c2 omega { partial "" } over { partial omega } pcpm
.EN
.EQ I
lineup ~=~ uc ^ left ( { 2 - uc } right ) pcpm - u2c2 T { partial "" } over { partial T } pcpm
.EN
.sp
The required partial with respect to period of the phase velocity partial is \f2numerically\f1
computed using the relation
.sp
.EQ I
{ partial ""} over { partial T }  pcopv ~=~ T left ( { ( pcopv ) sub { T+hT } -
( pcopv ) sub { T-hT } } over 2hT right )
.EN
.sp
In these expressions, the parameter $ m  $ can take on the values $ alpha $,  $  beta $, layer thickness or inverse Q.
and the subscript $ o $ represents the value in the parameter in the purely elastic
model.
.LP
If a causal Q is introduced, with a Futterman causality tied to a reference
angular frequency $ omega sub r , $ the elastic parameters will equal the
anelastic values only when the angular frequency $ omega $ equals $ omega sub r . $
This Q operator causes the medium velocity to be complex, e.g.,
.EQ I
v sub causal ^=^ v ^ left [ { 1 + { 1 }  over { pi Q sub v } ln left ({ omega over { omega sub r }  } right ) + i 1 over { 2 Q sub v  } }  right ]
.EN
The causal phase velocity $ c $ is given by applying the first terms of a Taylor series expansion:
.EQ I
c sub complex ^=^ c sub o 
^+^ { partial c sub o } over { partial alpha  }   ^  { partial alpha sub causal } over { partial { Q sub alpha sup -1 } } ^ Q  sub alpha sup -1 
^+^ { partial c sub o } over { partial beta  }   ^  { partial beta sub causal } over { partial { Q sub beta sup -1 }} ^  Q sub beta sup -1 
.EN
from which we obtain, assuming that the perturbation is small,
.sp
.EQ I (21)
c ~=~ co ~+~ piln sm
.EN
.sp
and the value of the spatial anelastic attenuation factor
.sp
.EQ I (22)
gamma ~=~ omc2 sm
.EN
In both (21) and (22) the summation extends over all all layers and the halfspace of the model.
To obtain this we expanded the propagation term $ e sup { - i omega r / c sub complex } $ to form $ e sup { - i omega r / c sub o } ^ e sup { - gamma r } $
.sp
These values are typically output by the programs \f3sregn96\f1
and \f3slegn96\f1 if causality is required. However, those programs
do not output the causal partial derivatives or group velocities since
these are not required for synthetic seismogram construction. The correct
causal relations for these parameters follow. A  $c $ or $ U $  without the subscript
$ o $   represents the causal value. 
Note that 
a partial with respect to $ Q sub beta sup -1 $ may involve a partial
of $ co $ with respect to $ alpha . $ To first order,
.sp
.EQ I (23)
pcpv ~=~ pcopv left ( 1 ~+~  { pilnqv }   right )
.EN
.sp
.EQ I (24)
pcpq = piln  ~ left ( pcopv ^ v  right )
.EN
.sp
.sp
.EQ I (25)
pgpq ~=~ om2c2   pcopv v 
.EN
for $ v ^=^ alpha $ or $ v ^=^ beta $.
An expression for $ pgpv $ cannot be obtained by simple first orper perturbation theory, since the change in $gamma $# depends on changes in $ { partial c } over { partial v }$, which are second order effects.
.sp
To obtain the expression for the causal group velocity, 
we express
.EQ I
U = uo ^+^ DELTA U ^=^ uo 
^+^ { partial U } over { partial Q sub alpha sup -1 } ^ DELTA Q sub alpha sup -1
^+^ { partial U } over { partial Q sub beta sup -1 } ^ DELTA Q sub beta sup -1
.EN
and use (20). From this we obtain
.EQ I (26)
U = uo left (  ~1 ~+~ ( 2 - uoco ) ccoco ~+~ 2gupo right )
.EN
.sp
where it is assumed that the higher order terms in $ ccoco  $
are negligible.
.sp
The partial derivatives of the causal group velocity 
.EQ I (27)
pupv mark = puopv left ( uuo - uoco ccoco + 2gupo  right ) 
.EN
.EQ I
lineup ~+~pcopv  uoco left ( { ~ -2 2gupo ~+~ uoco piqv }
.EN
.EQ I
lineup ~+~ left "" { ( 2 - uoco ) ^ [ { pilnqv - ccoco } ~~ ] + uoco ccoco  } right  )
.EN
.sp
.EQ I (28)
pupq mark ~=~ uoco (2 - uoco ) pcpq
.EN
.EQ I
lineup ~+~  { 1 over pi } uo2co2 ^ pcopv v 
.EN
.LP
An iterative linear inversion is performed because of the non-linear nature of
the problem. At any stage there is a current model which is used to
predict the observations and also a lack of fit. For simplicity, the compresional-wave velocity $ alpha $ is not directly inverted, rather it is tied to
the shear-wave velocity after determining the chane in the shear-wave velocity by one of the \f2assumptions\f1 made when running
\f3surf96\f1 or \f3joint96\f1.
In addition, the ratio $qbqa$ is fixed and the expression use in the inversion programs would look like
.EQ I
pcpq = piln left ( pcpb beta ~+~ pcpa alpha qbqa right )
.EN
rather than the single term in (24).
.LP
We may thus express the difference between observed and predicted values
in a linear model of changes in shear velocity model and changes in
the inverse Q model. Since the effect of compressional-wave Q may not be negligible,
the compresional-wave Q is related to the shear-wave Q by a ratio
$ Q sub alpha / Q sub beta $ for the layers.
Given these assumptions, the inversion can take on
two characters, non-causal and causal. To keep the notation general, the
partial derivatives with respect to layer velocity can either be 
causal or non-causal.
.LP
The differences in observed and predicted phase velocities are modeled as
.sp
.EQ I (29)
c sub obs - c sub pred mark ~=~ pcpb1 DELTA b1 + ... + pcpbn DELTA bn
.EN
.sp
.EQ I
lineup ~~~ ~+~ pcpq1 DELTA qb1 + ... + pcpqn DELTA qbn
.EN
.sp
.LP
The difference between observed and predicted group velocities
is given by
.sp
.EQ I (30)
U sub obs - U sub pred mark ~=~ pupb1 DELTA b1 + ... + pupbn DELTA bn
.EN
.sp
.EQ I
lineup ~~~ ~+~ pupq1 DELTA qb1 + ... + pupqn DELTA qbn
.EN
.sp
.LP
The difference between observed and gamma values predicted by
the current model is
.sp
.EQ I (31)
gamma sub obs - gamma sub pred mark ~=~ 
.EN
.sp
.EQ
lineup  pgpq1 DELTA qb1 + ... + pgpqn DELTA qbn
.EN
.sp
.LP
The inversion programs give the user the
choice of an (1) non-causal inversion, (2) a causal uncoupled
solution or (3) a causal-coupled solution.
The meaning of this is easily given in the following table:
.sp
.TS
center;
c c c c.
TERM	Non-Causal	Causal	Causal
		Uncoupled	Coupled
_
.sp
Phase Vel	$ c sub o $	(2)	(2)
.sp
Group Vel	$ U sub o $	(27)	(27)
.sp
$ pcpbc $	$ pcpb $	(23)	(23)
.sp
$ pupbc $	$ pupb $	(28)	(28)
.sp
$ pcpq $	0	0	(24)
.sp
$ pupq $	0	0	(28)
.sp
$ pgpq $	(7)	(7)	(25)
.sp
_
.sp
.TE
Zero entries in the table, indicate that the
respective partials are set to zero. Otherwise, the
partials  are given by the
function or  indicated equation.
.LP
.ne 3
\f3\s+23. Sphericity Corrections\s-2\f1 
.PP
Sphericity corrections are based on the work of Schwab and Knopoff (1972).
The concept is to convert the spherical earth model into
a flat earth model, compute the dispersion, and then
adjust the dispersion from the flat earth
model to make the spherical earth dispersion.
The sphericity correction for Love waves is exact, but that for the 
Rayleigh waves are approximate, both because of the nature
of the problem and also because the effects of gravitation are ignored.
.PP
In the discussion that follows, the subscript $ s $ represents the
spherical earth parameter, and the subscript $ f $ represents the
value used in the flat earth computations. If \f3obs.d\f1 file
indicates the use of a spherical earth model, then the program
\f3srfdis\f1 performs an earth flattening approximation, and computes
the dispersion for the equivalent flat earth model.
.LP
\f2Love Waves\f1
.PP
Let $ r $ be the radial distance from the center of the sphere, and
let the surface be given by $ r ^=^ a $. Also let $ z $ be the
depth from the free surface in the equivalent flat earth model.
The transformation used is
.sp
.EQ I (32)
z ~=~ a ^ ln ( a^/^r)
.EN
Given a spherical layer bounded by $r sub i$ and $r sub i-1$ , with $ r sub i-1 ^>^ r sub i $,
The thickness of the spherical  $i$'th layer is
is
.EQ I (33)
( h sub i ) sub s ~=~ r sub i-1 ^-^ r sub i
.EN
and the thickness of the transformed flat layer model is
.EQ I (34)
( h sub i ) sub f ~=~ 
 a ^ ln ( a^/^r sub i )
~-~ a ^ ln ( a^/^r sub i-1 )
.EN
The transformed   shear-wave  velocity and density in the equivalent flat layer model
are given by
.EQ I (35)
{ ( beta sub i ) sub f } bar
~=~
{  ( beta sub i ) sub s  } ~ dtp
.EN
.EQ I (36)
{ ( rho sub i ) sub f } bar
~=~
{  ( rho sub i ) sub s  } left ( { dtp  } right ) sup { -5}
.EN
The Love wave equation is from Schwab and Knopoff (1972).
Given this flat earth model, flat earth phase, $ c sub f $, and
group, $ U sub f $, are computed, as well as the
partial derivatives with respect to velocity, density, and layer
thickness. The program \f3srfdrl96\f1 computes the corresponding
spherical model values through the relations:
.EQ I (37)
c sub s ( omega ) ~=~  c sub f left [ 1 ^+^ (  3 c sub f ^/^ 2 a omega ) sup 2 right ] ^ sup {- 1^/^2}
.EN
.EQ I (38)
U sub s ( omega ) ~=~  U sub f left [ 1 ^+^ (  3 c sub f ^/^ 2 a omega ) sup 2 right ] sup {1^/^2}
.EN
To obtain the partials, the chain rule of differentiation is used:
.EQ I
{ partial c sub s } over { partial p sub s }
~=~ 
{ partial c sub s } over { partial c sub f }
~
{ partial c sub f } over { partial p sub f }
~
{ partial p sub f } over { partial p sub s }
.EN 
where $p$ is $ beta $ or $ rho $. The resulting expressions are
.EQ I (39)
left ( {{ partial c sub s } over { partial beta sub s } } right ) sub i
~=~
left ( 1 ^+^ ( 3 c sub f ^/^ 2 a omega ) sup 2 right ) sup { -3^/^2 }
~
left ( {{ partial c sub f } over { partial beta sub f } } right ) sub i
~
dtp
.EN
.EQ I (40)
left ( {{ partial U sub s } over { partial beta sub s } } right ) sub i
mark 
~=~
left {
left ( 1 ^+^ ( 3 c sub f ^/^ 2 a omega ) sup 2 right ) sup { 1^/^2 }
~
left ( {{ partial U sub f } over { partial beta sub f } } right ) sub i
.EN
.EQ I
lineup ~~~ ~+~
left ""
U sub f c sub f 
left ( {{ partial c sub f } over { partial beta sub f } } right ) sub i
( 3 ^/^ 2 a omega ) sup 2
left ( 1 ^+^ ( 3 c sub f ^/^ 2 a omega ) sup 2 right ) sup { -1^/^2 }
right }
~
dtp
.EN
.EQ I (41)
left ( {{ partial c sub s } over { partial h sub s } } right ) sub i
~=~
left ( 1 ^+^ ( 3 c sub f ^/^ 2 a omega ) sup 2 right ) sup { -3^/^2 }
~
left ( {{ partial c sub f } over { partial h sub f } } right ) sub i
~
{ a } over {  r sub i  }
.EN
.EQ I (42)
left ( {{ partial U sub s } over { partial h sub s } } right ) sub i
mark ~=~
left {
left ( 1 ^+^ ( 3 c sub f ^/^ 2 a omega ) sup 2 right ) sup { 1^/^2 }
~
left ( {{ partial U sub f } over { partial h sub f } } right ) sub i
.EN
.EQ
lineup
~+~
left ""
U sub f c sub f 
left ( {{ partial c sub f } over { partial h sub f } } right ) sub i
( 3 ^/^ 2 a omega ) sup 2
left ( 1 ^+^ ( 3 c sub f ^/^ 2 a omega ) sup 2 right ) sup { -1^/^2 }
right }
{ a   } over {  r sub i  }
.EN
.LP
\f2Rayleigh Waves\f1
.PP
Let $ r $ be the radial distance from the center of the sphere, and
let the surface be given by $ r ^=^ a $. Also let $ z $ be the
depth from the free surface in the equivalent flat earth model.
The transformation used is
.sp
.EQ I (43)
z ~=~ a ^ ln ( a^/^r)
.EN
The thickness of the spherical layer bounded by $ r sub i-1 ~ > r sub i $, is
is
.EQ I (44)
( h sub i ) sub s ~=~ r sub i-1 ^-^ r sub i
.EN
and the thickness of the transformed flat layer model is
.EQ I (45)
( d sub i ) sub f ~=~ 
 a ^ ln ( a^/^r sub i )
~-~ a ^ ln ( a^/^r sub i-1 )
.EN
The mean compressional- and shear-wave velocities  and density in the transformed flat layer model
are given by
.EQ I (46)
{ ( alpha sub i ) sub f } bar
~=~
{  ( alpha sub i ) sub s  } ~ dtp
.EN
.EQ I (47)
{ ( beta sub i ) sub f } bar
~=~
{  ( beta sub i ) sub s  } ~ dtp
.EN
.EQ I (48)
{ ( rho sub i ) sub f } bar
~=~
{  ( rho sub i ) sub s  } left ( { dtp } right ) sup { -2.275}
.EN
The exponent for the Rayleigh wave was determined empirically, e.g.,
\f2http://www.eas.slu.edu/eqc/eqc_cps/TUTORIAL/SPHERICITY/index.html\f1.
Given this flat earth model, flat earth phase, $ c sub f $, and
group, $ U sub f $, are computed, as well as the
partial derivatives with respect to velocity, density, and layer
thickness. The program \f3srfdrl\f1 computes the corresponding
spherical model values through the relations:
.EQ I (49)
c sub s ( omega ) ~=~  c sub f left [ 1 ^+^ ( c sub f ^/^ 2 a omega ) sup 2 right ] sup {-1^/^2}
.EN
.EQ I (50)
U sub s ( omega ) ~=~  U sub f left [ 1 ^+^ ( c sub f ^/^ 2 a omega ) sup 2 right ] sup {1^/^2}
.EN
To obtain the partials, the chain rule of differentiation is used:
.EQ I
{ partial c sub s } over { partial p sub s }
~=~ 
{ partial c sub s } over { partial c sub f }
~
{ partial c sub f } over { partial p sub f }
~
{ partial p sub f } over { partial p sub s }
.EN 
where $p$ is $ beta $ or $ rho $.
.EQ I (51)
left ( {{ partial c sub s } over { partial beta sub s } } right ) sub i
~=~
left ( 1 ^+^ (  c sub f ^/^ 2 a omega ) sup 2 right ) sup { -3^/^2 }
~
left ( {{ partial c sub f } over { partial beta sub f } } right ) sub i
~
dtp
.EN
.EQ I (52)
left ( {{ partial U sub s } over { partial beta sub s } } right ) sub i
mark
~=~
left {
left ( 1 ^+^ (  c sub f ^/^ 2 a omega ) sup 2 right ) sup { 1^/^2 }
~
left ( {{ partial U sub f } over { partial beta sub f } } right ) sub i
.EN
.EQ I
lineup
~~~ ~+~
left ""
U sub f c sub f 
left ( {{ partial c sub f } over { partial beta sub f } } right ) sub i
( 1^/^ 2 a omega ) sup 2
left ( 1 ^+^ (  c sub f ^/^ 2 a omega ) sup 2 right ) sup { -1^/^2 }
right }
~
dtp
.EN
.EQ I (53)
left ( {{ partial c sub s } over { partial h sub s } } right ) sub i
~=~
left ( 1 ^+^ (  c sub f ^/^ 2 a omega ) sup 2 right ) sup { -3^/^2 }
~
left ( {{ partial c sub f } over { partial h sub f } } right ) sub i
~
{ a } over {  r sub i  }
.EN
.EQ I (54)
left ( {{ partial U sub s } over { partial h sub s } } right ) sub i
mark ~=~
left {
left ( 1 ^+^ (  c sub f ^/^ 2 a omega ) sup 2 right ) sup { 1^/^2 }
~
left ( {{ partial U sub f } over { partial h sub f } } right ) sub i
.EN
.EQ I
lineup
~~~ ~+~
left ""
U sub f c sub f 
left ( {{ partial c sub f } over { partial h sub f } } right ) sub i
( 1 ^/^ 2 a omega ) sup 2
left ( 1 ^+^ (  c sub f ^/^ 2 a omega ) sup 2 right ) sup { -1^/^2 }
right }
{ a   } over {  r sub i  }
.EN
.LP
\f3\s+4References\s-4\f1

.LP
.in +3m
.sp
.ti -3m
Keilis-Borok, V. I., A. L. Levshin, T. B. Yanovskaya,
A. V. Lander, B. G. Bukchin, M. P. Barmin, L. I. Ratnikova,
and E. N. Its (1989).
\f2Seismic surface waves in a laterally inhomogeneous earth\f1,
Kluwer Academic Publishers, Dordrecht.
.sp
.ti -3m 
Rodi, W. L., P. Glover, T. M. C. Li and S. S. Alexander (1975). A fast, accurate method for computing group-velocity partial derivatives for Rayleigh and Love modes, 
.I
Bull. Seism. Soc. Am.
.R
.B 65,
1105-1114.
.sp
.ti -3m
Schwab, F. A., and L. Knopoff (1972). Fast surface wave and free mode
computations,
in
\f2Methods in Computational Physics, Vol. 11: Seismology: Surface Waves and Earth Oscillations,
.R
B. A. Bolt, ed, Academic Press, New York, pp. 87-180.
.sp
.ti -3m
Wang, C.-Y. (1981). Wave Theory for Seismogram Synthesis,
\f2Ph. D. Dissertation,\f1 Saint Louis University, St. Louis, MO.
.in -3m
.LP
.EQ
delim $$
.EN
.ce 1
\s+4\f3DEFINITIONS\s-4\f1
.LP
\f3Gamma Values\f1 - a measure of the non-geometrical
attenuation of a signal due to anelastic processes which is of the form $ exp ( - gamma r $). If a 
two-station technique is used, and if surface-wave propagation is
assumed, then the interstation $ gamma $ is defined by the
relation:
.EQ
gamma ~=~    ln left ( {
 { A sub 1 r sub 1 sup { 1 ^/^ 2 } }
^/^ 
{ A sub 2 r sub 2 sup { 1 ^/^ 2 } }
} right )
/ ^ ( r sub 2 ^-^ r sub 1 )
.EN
where $ A sub j $ is the instrument
corrected surface-wave spectral amplitude observed at distance $ r sub j $.
Both spectral amplitude observations are measured at the same frequency. 
To obtain clean estimates of the spectral amplitude estimates, multiple
filter analysis or phase match filter techniques can be used.
.sp
The units are $ roman km sup -1 $, $ roman m sup -1 $ or $ roman ft sup -1 $.
These must be in the same units as the velocities input.
.LP
\f3Group Velocity\f1 - velocity of energy propagation. This can be estimated
graphically or by using multiple filter techniques, which
bandpass filter the surface wave, and compute the envelope of the
resulting function. The group velocity is the obtained by dividing the 
epicentral distance by the time of arrival of the envelope maximum.  
.sp
The units are km/sec or m/sec or sec.
.LP
\f3Phase Velocity\f1 - velocity of a given phase. Usually measured
from phase spectrum at a given frequency. If source and instrument phase
is known, a single station technique can be used. For
data acquired along the same azimuth, two station interstation
phase is easily obtained if the instrument phase is known.
If more than two stations are available, a stacking technique
can be used to reduce the problems of spatial aliasing.
.sp
The units are km/sec or m/sec or sec.

