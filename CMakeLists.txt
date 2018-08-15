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
# STRING (REGEX REPLACE "/RTC(su|[1su])" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
add_definitions(/Gy /O1 /fp:fast /GR- /GS-)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /ENTRY:demo /SUBSYSTEM:Windows  /OPT:REF  /OPT:ICF=10 /VERBOSE /NODEFAULTLIB") # /ENTRY:demo /NODEFAULTLIB
add_executable(etd etd.c)

# TODO: remove ucrt.lib if DEBUG is undefined in-source
target_link_libraries(etd opengl32.lib ucrt.lib User32.lib) #ucrt.lib  libcmt.lib version.lib imm32.lib winmm.lib
