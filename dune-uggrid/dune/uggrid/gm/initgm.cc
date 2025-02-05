// SPDX-FileCopyrightText: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LGPL-2.1-or-later
// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/****************************************************************************/
/*																			*/
/* File:	  initgm.c														*/
/*																			*/
/* Purpose:   call the init routines of the grid manager module                         */
/*																			*/
/* Author:	  Henrik Rentz-Reichert                                                                                 */
/*			  Institut fuer Computeranwendungen III                                                 */
/*			  Universitaet Stuttgart										*/
/*			  Pfaffenwaldring 27											*/
/*			  70569 Stuttgart												*/
/*			  email: ug@ica3.uni-stuttgart.de						        */
/*																			*/
/* History:   27.02.95 begin, ug version 3.0								*/
/*																			*/
/* Remarks:                                                                                                                             */
/*																			*/
/****************************************************************************/


/****************************************************************************/
/*																			*/
/* include files															*/
/*			  system include files											*/
/*			  application include files                                                                     */
/*																			*/
/****************************************************************************/

/* ANSI-C includes */
#include <config.h>
#include <cstdio>

/* low module */
#include <dune/uggrid/low/misc.h>
#include <dune/uggrid/low/ugstruct.h>
#include <dune/uggrid/low/ugtypes.h>

/* gm module */
#include "gm.h"
#include "algebra.h"
#include "cw.h"
#include "ugm.h"
#include "ugio.h"
#include "elements.h"
#include "refine.h"
#include "rm.h"

/* own header */
#include "initgm.h"


USING_UG_NAMESPACE
USING_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                          */
/* definition of variables global to this source file only (static!)        */
/*                                                                          */
/****************************************************************************/

/****************************************************************************/
/*
   InitGm - Call the inits for the grid manager module

   SYNOPSIS:
   INT InitGm ();

   PARAMETERS:
   .  void

   DESCRIPTION:
   This function calls the inits for the grid manager module.

   RETURN VALUE:
   INT
   .n     0 if ok
   .n     1 if some error occurred.
 */
/****************************************************************************/

INT NS_DIM_PREFIX InitGm ()
{
  INT err;

  /* cw.c */
  if ((err=InitCW())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* elements.c */
  if ((err=PreInitElementTypes())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* ugm.c */
  if ((err=InitUGManager())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* ugio.c */
  if ((err=InitUgio())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* rm.c */
  if ((err=InitRuleManager())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  /* set config variables for the script */
  if (SetStringValue("conf:dim",(DOUBLE)DIM))
    return(__LINE__);

  return (0);
}


INT NS_DIM_PREFIX ExitGm()
{
  INT err;

  /* ugm.c */
  if ((err=ExitUGManager())!=0)
  {
    SetHiWrd(err,__LINE__);
    return (err);
  }

  return 0;
}
