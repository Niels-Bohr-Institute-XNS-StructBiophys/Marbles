'''*****************************************************************************
Copyright (C) 2020  Niels Bohr Institute

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
*****************************************************************************'''

try:
    from pymol.cgo import *
    from pymol import cmd
    from pymol.vfont import plain
    import __main__
    __main__.pymol_argv = [ 'pymol', '-e' ]
    
except ModuleNotFoundError:
    print("Usage: pymol plotdisc.py -- path/fit.wit path/model.pdb")
    exit(-1)

import sys

def parse():

    if len(sys.argv) == 3:
        f = sys.argv[1]
        s = sys.argv[2]
    else:
        print("Usage: pymol plotdisc.py -- path/fit.wit path/model.pdb")
        exit(-1)

    return f, s
#-------------------------------------------------------------------------------

def parse_wif( wif_file ):

    params = {}
    with open( wif_file, 'r' ) as f:
        lines = f.readlines()

    for l in lines:

        if 'Major semiaxis' in l:
            params['mjs'] = float( l.split(' ')[-2] )
        if 'Minor semiaxis' in l:
            params['mis'] = float( l.split(' ')[-2] )
        if 'Width of belt' in l:
            params['wb']  = float( l.split(' ')[-2] )
        if 'Height of bilayer' in l:
            params['hb']  = float( l.split(' ')[-2] ) / 2.
        if 'Height of hydrophobic bilayer' in l:
            params['hhb']  = float( l.split(' ')[-2] ) / 2.
        if 'HeightOfBelt' in l:
            n = [ x for x in l.split(' ') if x != '' ]
            params['hbb']  = float( n[1] ) / 2.
        if 'Height of methyl layer' in l:
            params['hm']  = float( l.split(' ')[-2] ) / 2.

    return params
#-------------------------------------------------------------------------------

def signOfFloat(f):
        if f < 0: return -1
        if f > 0: return 1
        return 0
#-------------------------------------------------------------------------------

def sqC(v, n):
        return signOfFloat(math.cos(v)) *  math.pow(math.fabs(math.cos(v)), n)
#-------------------------------------------------------------------------------

def sqCT(v, n, alpha):
        return alpha + sqC(v, n)
#-------------------------------------------------------------------------------

def sqS(v, n):
        return signOfFloat(math.sin(v)) * math.pow(math.fabs(math.sin(v)), n)
#-------------------------------------------------------------------------------

def sqEllipsoid(x, y, z, a1, a2, a3, u, v, n, e):
        m=n
        x = a1 * sqC(u, n) * sqC(v, e) + x
        y = a2 * sqC(u, n) * sqS(v, e) + y
        z = a3 * sqS(u, m) + z
        nx = sqC(u, 2 -n) * sqC(v, 2 - e)  / (a1)
        ny = sqC(u, 2 -n) * sqS(v, 2 - e)  / (a2)
        nz = sqS(u, 2 -m)   / (a3)
        return x, y, z, nx, ny, nz
#-------------------------------------------------------------------------------

def makeSuperQuadricEllipsoid(x, y, z, a1, a2, a3, n = 0, e = 1, u1 = -1.5708, u2 = 1.5708, v1 = -3.14159,
                              v2 = 3.14159, u_segs = 30, v_segs = 30, color = [0.5, 0.5, 0.5], alpha = 0.5 ):

        r, g, b = color
        r /= 255
        g /= 255
        b /= 255

        # Calculate delta variables */
        dU = (u2 - u1) / u_segs
        dV = (v2 - v1) / v_segs

        o = [BEGIN, TRIANGLES]

        U = u1
        for Y in range(0, u_segs):
                # Initialize variables for loop */
                V = v1
                for X in range(0, v_segs):
                        # VERTEX #1 */
                        x1, y1, z1, n1x, n1y, n1z = sqEllipsoid(x, y, z, a1, a2, a3, U, V, n, e)
                        x2, y2, z2, n2x, n2y, n2z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V, n, e)
                        x3, y3, z3, n3x, n3y, n3z = sqEllipsoid(x, y, z, a1, a2, a3, U + dU, V + dV, n, e)
                        x4, y4, z4, n4x, n4y, n4z = sqEllipsoid(x, y, z, a1, a2, a3, U, V + dV, n, e)

                        o.extend([ALPHA, alpha, COLOR, r, g, b, NORMAL, n1x, n1y, n1z, VERTEX, x1, y1, z1])
                        o.extend([ALPHA, alpha, COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([ALPHA, alpha, COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])
                        o.extend([ALPHA, alpha, COLOR, r, g, b, NORMAL, n2x, n2y, n2z, VERTEX, x2, y2, z2])
                        o.extend([ALPHA, alpha, COLOR, r, g, b, NORMAL, n3x, n3y, n3z, VERTEX, x3, y3, z3])
                        o.extend([ALPHA, alpha, COLOR, r, g, b, NORMAL, n4x, n4y, n4z, VERTEX, x4, y4, z4])

                        # Update variables for next loop */
                        V += dV
                # Update variables for next loop */
                U += dU
        o.append(END)
        return o
#-------------------------------------------------------------------------------

f, s = parse()
params = parse_wif( f )

major_wbelt  = params['mjs'] + params['wb']
minor_wbelt  = params['mis'] + params['wb']

cmd.load( s, 'model' )
cmd.show('spheres')
cmd.set('sphere_scale', 1.5)
cmd.set('sphere_color', 'grey60')

cmd.load_cgo( makeSuperQuadricEllipsoid(0, 0, 0, params['mjs'],       params['mis'],       params['hb'],  color = [213., 94., 0.] ),   "NHeads" )
cmd.load_cgo( makeSuperQuadricEllipsoid(0, 0, 0, params['mjs'] + 0.1, params['mis'] + 0.1, params['hhb'], color = [230., 159., 0.] ),  "NTails" )
cmd.load_cgo( makeSuperQuadricEllipsoid(0, 0, 0, params['mjs'] + 0.2, params['mis'] + 0.2, params['hm'],  color = [240., 228., 66.] ), "NMethyls" )
cmd.load_cgo( makeSuperQuadricEllipsoid(0, 0, 0, major_wbelt,         minor_wbelt,         params['hbb'], color = [0., 158., 115.] ),  "NBelt" )
cmd.zoom('model_', 50)
