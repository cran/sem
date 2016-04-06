/*
 * =====================================================================================
 *
 *       Filename:  csem.h
 *
 *    Description:  Header files for csem 
 *
 *        Version:  1.0
 *        Created:  Tue 27 Dec 2011 00:36:32 EST
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Zhenghua Nie (ZHN), zhenghua.nie@gmail.com
 *        Company:  McMaster University
 *
 *    Copyright (C) 2011 Zhenghua Nie. All Rights Reserved.
 *    This code is published under GNU GENERAL PUBLIC LICENSE.
 *
 *    This program is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License,  or
 *    (at your option) any later version.
 *      
 *    This program is distributed WITHOUT ANY WARRANTY. See the
 *    GNU General Public License for more details.
 *           
 *    If you do not have a copy of the GNU General Public License,  
 *    write to the Free Software Foundation, Inc., 
 *    59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 *              
 *        
 * =====================================================================================
 */

#ifndef __CSEM_HPP__
#define __CSEM_HPP__

#include <iomanip> 
#include <iostream> 
#include <fstream> 
#include <sstream>
#include <stdio.h> 
#include <string.h> 
#include <string> 
#include <assert.h> 
#include <sys/types.h>     
#include <sys/stat.h>    
#include <math.h> 
#include <cmath> 
#include <limits>

#include <R.h>
#include <Rdefines.h>
#include <R_ext/Utils.h>
#include <float.h>		/* for DBL_MAX */
#include <Rmath.h>
#include <R_ext/Applic.h>
#include <R_ext/RS.h>	       	/* for Memcpy */
#include <R_ext/Parse.h>
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>

#include "csemnlm.h"

#endif
