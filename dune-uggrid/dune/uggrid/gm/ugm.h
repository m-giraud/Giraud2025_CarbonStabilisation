// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
// SPDX-FileCopyrightInfo: Copyright © DUNE Project contributors, see file LICENSE.md in module root
// SPDX-License-Identifier: LGPL-2.1-or-later
/*! \file ugm.h
 * \ingroup gm
 */


/****************************************************************************/
/*                                                                          */
/* File:      ugm.h                                                         */
/*                                                                          */
/* Purpose:   unstructured grid manager header file                         */
/*            internal interface in grid manager module                     */
/*                                                                          */
/* Author:    Peter Bastian                                                                                                 */
/*            Interdisziplinaeres Zentrum fuer Wissenschaftliches Rechnen   */
/*            Universitaet Heidelberg                                                                               */
/*            Im Neuenheimer Feld 368                                                                               */
/*                        6900 Heidelberg                                                                                               */
/*                        internet: ug@ica3.uni-stuttgart.de                                            */
/*                                                                                                                                                      */
/* History:   09.03.92 begin, ug version 2.0                                                            */
/*                        13.12.94 begin, ug version 3.0                                                                */
/*                                                                                                                                                      */
/* Remarks:                                                                                                                             */
/*                                                                                                                                                      */
/****************************************************************************/


/****************************************************************************/
/*                                                                                                                                                      */
/* auto include mechanism and other include files                                                       */
/*                                                                                                                                                      */
/****************************************************************************/

#ifndef __UGM__
#define __UGM__

#include <dune/uggrid/low/namespace.h>
#include <dune/uggrid/low/ugtypes.h>

#include "gm.h"
#include "dlmgr.h"

START_UGDIM_NAMESPACE

/****************************************************************************/
/*                                                                                                                                                      */
/* defines in the following order                                                                                       */
/*                                                                                                                                                      */
/*                compile time constants defining static data size (i.e. arrays)        */
/*                other constants                                                                                                       */
/*                macros                                                                                                                        */
/*                                                                                                                                                      */
/****************************************************************************/

#define MAX_PAR_DIST    1.0E-6          /* max.dist between different parameter */

/****************************************************************************/
/*                                                                                                                                                      */
/* data structures exported by the corresponding source file                            */
/*                                                                                                                                                      */
/****************************************************************************/


/****************************************************************************/
/*                                                                                                                                                      */
/* definition of exported global variables                                                                      */
/*                                                                                                                                                      */
/****************************************************************************/

/****************************************************************************/
/*                                                                                                                                                      */
/* function declarations                                                                                                        */
/*                                                                          */
/****************************************************************************/

/* init */
INT              InitUGManager                  (void);
INT              ExitUGManager ();

/* object handling */
GM_OBJECTS              GetFreeOBJT                    (void);
INT              ReleaseOBJT                    (GM_OBJECTS);

/* create basic objects */
/** \todo Commented out because no definitions exist for these declarations

   #if defined(UG_DIM_2)
   int      GetElemLink (NODE *from, NODE *to, ELEMENT *theElement);
   ELEMENT *NbElem     (const ELEMENT *theElement, int i);
   void     Set_NbElem (ELEMENT *theElement, int i, ELEMENT *Neighbor);
   #endif
 */

#ifdef ModelP
EDGE * CreateEdge (GRID *theGrid, ELEMENT *theElement, INT i, bool with_vector);
#endif
ELEMENT * CreateElement          (GRID *theGrid, INT tag, INT objtype,
                                  NODE **nodes, ELEMENT *Father, bool with_vector);
INT         CreateSonElementSide    (GRID *theGrid, ELEMENT *theElement,
                                     INT side, ELEMENT *theSon, INT son_side);

GRID            *CreateNewLevel                 (MULTIGRID *theMG);

/* dispose basic objects */
INT              DisposeElement                 (GRID *theGrid, ELEMENT *theElement);
INT              DisposeTopLevel                (MULTIGRID *theMG);
INT              DisposeNode                    (GRID *theGrid, NODE *theNode);

/* miscellaneous */
INT              FindNeighborElement    (const ELEMENT *theElement, INT Side, ELEMENT **theNeighbor, INT *NeighborSide);
INT             CheckOrientation                (INT n, VERTEX **vertices);
INT             CheckOrientationInGrid  (GRID *theGrid);

NODE        *CreateSonNode          (GRID *theGrid, NODE *FatherNode);
NODE            *CreateMidNode                  (GRID *theGrid, ELEMENT *theElement, VERTEX *theVertex, INT edge);
NODE            *GetCenterNode                  (const ELEMENT *theElement);
NODE        *CreateCenterNode       (GRID *theGrid, ELEMENT *theElement, VERTEX *theVertex);

#ifdef UG_DIM_3
NODE            *CreateSideNode                 (GRID *theGrid, ELEMENT *theElement, VERTEX *theVertex, INT side);
NODE            *GetSideNode                    (const ELEMENT *theElement, INT side);
#endif
INT          GetSideIDFromScratch   (ELEMENT *theElement, NODE *theNode);
NODE        *GetMidNode             (const ELEMENT *theElement, INT edge);
INT                     GetNodeContext                  (const ELEMENT *theElement, NODE **theElementContext);
void            GetNbSideByNodes                (ELEMENT *theNeighbor, INT *nbside, ELEMENT *theElement, INT side);


void *GetMemoryForObject (MULTIGRID *mg, INT size, INT type);
INT PutFreeObject (MULTIGRID *mg, void *object, INT size, GM_OBJECTS type);

/* determination of node classes */
INT             ClearNodeClasses                        (GRID *theGrid);
INT             SeedNodeClasses                     (ELEMENT *theElement);
INT             PropagateNodeClasses            (GRID *theGrid);
INT             ClearNextNodeClasses            (GRID *theGrid);
INT             SeedNextNodeClasses             (ELEMENT *theElement);
INT             PropagateNextNodeClasses        (GRID *theGrid);
INT             MaxNextNodeClass                        (const ELEMENT *theElement);
INT             MinNodeClass                            (const ELEMENT *theElement);
INT             MinNextNodeClass                        (const ELEMENT *theElement);

END_UGDIM_NAMESPACE

#endif
