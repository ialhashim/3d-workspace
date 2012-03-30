//=============================================================================
// Copyright (C) 2001-2005 by Computer Graphics Group, RWTH Aachen
// Copyright (C) 2011 by Graphics & Geometry Group, Bielefeld University
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Library General Public License
// as published by the Free Software Foundation, version 2.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
//=============================================================================


#ifndef MESH_IO_H
#define MESH_IO_H


//== INCLUDES =================================================================

#include <string>
#include "Surface_mesh.h"


//=============================================================================

// For Windows warnings..
#ifdef WIN32
    inline FILE* myfopen(const char * name, const char * mode)
    {FILE *f; fopen_s(&f, name, mode); return f;}
    #define fopen(x,y) myfopen(x,y)
    #define sscanf sscanf_s
    #define fscanf fscanf_s
#endif

bool read_mesh(Surface_mesh& mesh, const std::string& filename);
bool read_off(Surface_mesh& mesh, const std::string& filename);
bool read_obj(Surface_mesh& mesh, const std::string& filename);
bool read_stl(Surface_mesh& mesh, const std::string& filename);

bool write_mesh(const Surface_mesh& mesh, const std::string& filename);
bool write_off(const Surface_mesh& mesh, const std::string& filename);


//=============================================================================
#endif
//=============================================================================
