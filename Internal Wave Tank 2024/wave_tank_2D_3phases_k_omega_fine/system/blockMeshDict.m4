/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  v2312                                 |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    `format'      ascii;
    class       dictionary;
    object      blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
changecom(//)changequote([,])

define(calc, [esyscmd(perl -e 'printf ($1)')])

define(VCOUNT, 0)

define(vlabel, [[// ]Vertex $1 = VCOUNT define($1, VCOUNT)define([VCOUNT], incr(VCOUNT))])


define(hex2D, hex ($1b $2b $3b $4b $1t $2t $3t $4t))

define(quad2D, ($1b $2b $3b $4b))

define(topquad2D, ($1t $2t $3t $4t))

define(backQuad, ($1b $4b $3b $2b))

scale 0.01; // in cm

// defining dimensions

define(length, 100.0)
define(height, 20.0)
define(obj_h,13.6)
define(obj_wbase,3.9)
define(obj_wtop,1.3)
define(L1, 60.0)
define(Nx1, 120 ) // = 
define(Nx2, 73 ) // = 
define(Ny1, 28 ) // = 
define(Ny2, 13 ) // = 
define(NxObj, 3)
// x 
define(X00, 0.0)
define(X01, L1)
define(X02, calc(L1+ obj_wbase))
define(X03, length)
define(X10, 0.0)
define(X11, calc( L1 + (obj_wbase-obj_wtop)/2 ) )
define(X12, calc(X11 + obj_wtop) )
define(X13, length)
define(X20, 0.0)
define(X21, X11)
define(X22, X12)
define(X23, length)
// y
define(Y00, 0.0)
define(Y01, 0.0)
define(Y02, 0.0)
define(Y03, 0.0)
define(Y10, obj_h)
define(Y11, obj_h)
define(Y12, obj_h)
define(Y13, obj_h)
define(Y20, height)
define(Y21, height)
define(Y22, height)
define(Y23, height)


define(origox, 0.0)
define(origoy, 0.0)


define(nrcellsx, 10)
define(nrcellsy, 30)

define(Z1, 0.0)
define(Z2, 1.0)

// defining vertices
vertices
(
    (X00 Y00 Z1) vlabel(a1b)
    (X10 Y10 Z1) vlabel(a2b)
    (X11 Y11 Z1) vlabel(a3b)
    (X01 Y01 Z1) vlabel(a4b)
    (X00 Y00 Z2) vlabel(a1t)
    (X10 Y10 Z2) vlabel(a2t)
    (X11 Y11 Z2) vlabel(a3t)
    (X01 Y01 Z2) vlabel(a4t)

	(X02 Y02 Z1) vlabel(b1b)
	(X12 Y12 Z1) vlabel(b2b)
	(X13 Y13 Z1) vlabel(b3b)
	(X03 Y03 Z1) vlabel(b4b)
	(X02 Y02 Z2) vlabel(b1t)
	(X12 Y12 Z2) vlabel(b2t)
	(X13 Y13 Z2) vlabel(b3t)
	(X03 Y03 Z2) vlabel(b4t)

	// (X12 Y12 Z1) vlabel(c1b)
	(X22 Y22 Z1) vlabel(c2b)
	(X23 Y23 Z1) vlabel(c3b)
	// (X13 Y13 Z1) vlabel(c4b)
	//(X12 Y12 Z2) vlabel(c1t)
	(X22 Y22 Z2) vlabel(c2t)
	(X23 Y23 Z2) vlabel(c3t)
	//(X13 Y13 Z2) vlabel(c4t)
	
	//(X11 Y11 Z1) vlabel(d1b)
	//(X21 Y21 Z1) vlabel(d2b)
	//(X22 Y22 Z1) vlabel(d3b)
	//(X12 Y12 Z1) vlabel(d4b)
	//(X11 Y11 Z2) vlabel(d1t)
	//(X21 Y21 Z2) vlabel(d2t)
	//(X22 Y22 Z2) vlabel(d3t)
	//(X12 Y12 Z2) vlabel(d4t)
	
	//(X10 Y10 Z1) vlabel(e1b)
	(X20 Y20 Z1) vlabel(e2b)
	(X21 Y21 Z1) vlabel(e3b)
	//(X11 Y11 Z1) vlabel(e4b)
	//(X10 Y10 Z2) vlabel(e1t)
	(X20 Y20 Z2) vlabel(e2t)
	(X21 Y21 Z2) vlabel(e3t)
	//(X11 Y11 Z2) vlabel(e4t)
	
	
	
);

blocks
(
    // block A
    //hex2D(a1, a2, a3, a4)
    hex (a1b a4b a3b a2b a1t a4t a3t a2t)
    (Nx1 Ny1 1)
    simpleGrading (1 1 1)
    
    // block B
    //hex2D(b1, b2, b3, b4)
    hex (b1b b4b b3b b2b b1t b4t b3t b2t)
    (Nx2 Ny1 1)
    simpleGrading (1 1 1)
    
    // block C
    hex (b2b b3b c3b c2b b2t b3t c3t c2t)
    (Nx2 Ny2 1)
    simpleGrading (1 1 1)
    
    // block D
    hex (a3b b2b c2b e3b a3t b2t c2t e3t)
    (NxObj Ny2 1)
    simpleGrading (1 1 1)
    
    // block E
    hex (a2b a3b e3b e2b a2t a3t e3t e2t)
    (Nx1 Ny2 1)
    simpleGrading (1 1 1)
);

boundary
(
	rightWall
	{
		type wall;
		faces
		(
		    (b3b b3t b4t b4b)
		    (c3b c3t b3t b3b)
		);
	}
	leftWall
	{
		type wall;
		faces
		(
		    (a1b a1t a2t a2b)
		    (a2b a2t e2t e2b)
		    
		);
	}
	
	lowerWall
	{
		type wall;
		faces
		(
		    (a1b a4b a4t a1t)
		    (a3b a3t a4t a4b)
		    (a3b b2b b2t a3t)
		    (b1b b1t b2t b2b)
		    (b1b b1t b4t b4b)
		);
	}
	atmosphere
	{
		type patch;
		faces
		(
		    (e3b e3t e2t e2b)
		    (e3b c2b c2t e3t)
		    (c3b c3t c2t c2b)
		);
	}

	
);

mergePatchPairs
(
);

