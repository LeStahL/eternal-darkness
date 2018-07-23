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

include_directories("C:\\Users\\NR4\\Downloads\\SDL2-devel-2.0.8-VC\\SDL2-2.0.8\\include")
add_definitions(/Gy /O1 /fp:fast /GR- /QIfist)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /subsystem:windows /NODEFAULTLIB /OPT:REF /ENTRY:demo /OPT:ICF=10 /VERBOSE")
add_executable(etd etd.c)

# TODO: remove ucrt.lib if DEBUG is undefined in-source
target_link_libraries(etd "C:\\Users\\NR4\\Downloads\\SDL2-2.0.8\\SDL2-2.0.8\\build\\Release\\SDL2-static.lib" version.lib imm32.lib winmm.lib  opengl32.lib libcmt.lib ucrt.lib)
