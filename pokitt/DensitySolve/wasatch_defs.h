#ifndef WASATCH_DEFS_H
#define WASATCH_DEFS_H

/*
 * The MIT License
 *
 * Copyright (c) 1997-2016 The University of Utah
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

// Define for having POKITT (which also means CANTERA exists) will be found below 
// (if it was specified during configure).

/* 
   For reference, (if Wasatch is enabled, then) you are linking versus this SpatialOps (W3P):

   -I/Users/josh/development/builds/uintah/opt/Wasatch3P/install/SpatialOps/include

   -L/Users/josh/development/builds/uintah/opt/Wasatch3P/install/SpatialOps/lib/spatialops

   (Note, the above info is used to tell when the user configures
   versus a new SpatialOps library, so please leave it here.)
*/

#define HAVE_POKITT 1

#endif // WASATCH_DEFS_H
