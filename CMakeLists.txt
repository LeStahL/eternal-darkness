# Eternal Darkness by Team210
# Copyright (C) 2018  Alexander Kraus <nr4@z10.info>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

cmake_minimum_required(VERSION 2.8)

project(eternal-darkness)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

if(WIN32)
    add_definitions(/Gy /O1 /fp:fast /GR- /GS- /MT)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:demo /SUBSYSTEM:Windows  /OPT:REF  /OPT:ICF=10 /VERBOSE /NODEFAULTLIB") 
    add_executable(debug debug.c)
    add_executable(etd etd.c)
elseif(UNIX)
    add_definitions(-Oz)
    add_executable(etd etd-linux.c)
endif()


if(WIN32)
    target_link_libraries(etd opengl32.lib Winmm.lib)
    target_link_libraries(debug opengl32.lib Winmm.lib User32.lib ucrt.lib msvcrt.lib)
elseif(UNIX)
    target_link_libraries(etd X11 GL asound m pthread dl)
endif()
