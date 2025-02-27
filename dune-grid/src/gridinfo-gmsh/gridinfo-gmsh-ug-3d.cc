// SPDX-FileCopyrightText: Copyright (C) DUNE Project contributors, see file COPYING in module root
// SPDX-License-Identifier: LicenseRef-GPL-2.0-only-with-DUNE-exception
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string>

#include <dune/grid/uggrid.hh>

const std::string programName = "dune-gridinfo-gmsh-ug-3d";
typedef Dune::UGGrid<3> Grid;

#include <dune/grid/utility/gridinfo-gmsh-main.hh>
