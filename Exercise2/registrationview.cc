//=============================================================================
//
//   Code framework for the lecture
//
//   "Surface Representation and Geometric Modeling"
//
//   Mark Pauly, Mario Botsch, Balint Miklos, and Hao Li
//
//   Copyright (C) 2007 by  Applied Geometry Group and
//                          Computer Graphics Laboratory, ETH Zurich
//
//-----------------------------------------------------------------------------
//
//                                License
//
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 51 Franklin Street, Fifth Floor,
//   Boston, MA  02110-1301, USA.
//
//=============================================================================

#include <cstdio>

#include "RegistrationViewer.hh"

using namespace std;

int main(int argc, char **argv)
{
    glutInit(&argc, argv);

    RegistrationViewer window("Registration Viewer", 512, 512);

    if (argc>2)
    {
        std::vector<std::string> filenames;
        for(int i = 2; i < argc; i++) filenames.push_back( std::string( argv[i]) );

        window.set_output(argv[1]);

        if( window.open_meshes(filenames) )
        {
            glutMainLoop();
        }
        else
        {
            printf("Could not load all files\n");
        }
    }
    else
    {
        printf("Usage: %s <ouptput-points> <meshes>*\n", argv[0]);
    }
}
