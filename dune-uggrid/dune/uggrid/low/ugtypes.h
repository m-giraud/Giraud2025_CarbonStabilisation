// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LGPL-2.1-or-later

/** \file
    \brief Globally set data types

    The basic types (normally 'short', 'int', 'float' and 'double') are
    replaced by 'SHORT', 'INT', 'FLOAT' and 'DOUBLE'. The type 'DOUBLE'
    is used for all Cartesian coordinates of the (x,y[,z])-directions of
    the grids.
 */

#ifndef UGTYPES_H
#define UGTYPES_H

#include "namespace.h"

START_UG_NAMESPACE

/* standard types */

typedef short SHORT;

/* these types are used for several bitfields. I'd guess that it needs
   at least 32 bits... */
typedef int INT;
typedef unsigned int UINT;

typedef float FLOAT;
typedef double DOUBLE;

END_UG_NAMESPACE

#endif
