// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LGPL-2.1-or-later

/** \file
    \brief Provide the usual preprocessor-defines for the dimension and complain
    if the dimension was set incorrectly.
 */


#ifndef DIMENSION_H
#define DIMENSION_H


#ifdef UG_DIM_2
#ifdef UG_DIM_3
#error ****    define EITHER dimension UG_DIM_2 OR UG_DIM_3       ****
#endif
#define DIM 2
#define DIM_OF_BND 1
#endif

#ifdef UG_DIM_3
#define DIM 3
#define DIM_OF_BND 2
#endif

#ifndef UG_DIM_2
#ifndef UG_DIM_3
#error ****    define at least dimension two OR three        ****
#endif
#endif

#endif
