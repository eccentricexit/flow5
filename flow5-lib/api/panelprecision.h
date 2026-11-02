/****************************************************************************

    flow5 application
    Copyright (C) 2025 André Deperrois 
    
    This file is part of flow5.

    flow5 is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License,
    or (at your option) any later version.

    flow5 is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty
    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
    See the GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with flow5.
    If not, see <https://www.gnu.org/licenses/>.


*****************************************************************************/



#pragma once


// All panel methods
#define DISTANCEPRECISION   1.e-6
#define INPLANEPRECISION    1.e-6
// A field point lies in a triangle's plane, for the doublet's velocity and potential and the source's potential and velocity (Panel3),
// when its distance to the plane is at the rounding of the coordinates: a point in the plane of a triangle of a 6 m wing is off it by
// 0 to 2e-21 m, and the coordinates of an aircraft below 10 m round at 1e-15 m.
#define PLANEROUNDING       1.e-12
#define MINREFAREA          1.e-10
#define MINREFDIMENSION     1.0e-6

//specific to M. Carley's method
#define VERTEXANGLEPRECISION  1.e-6   // radians
#define SIDELENGTHPRECISION   1.e-6
#define INTEGRALPRECISION     1.e-14


