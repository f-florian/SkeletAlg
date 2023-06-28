#!/bin/bash

filename=$1
paramName="length"
meshPath=~/Documents/dev/meshes/

(gmsh -2 -o ${meshPath}${filename}00.msh -format msh4 ${filename}.geo -setnumber $paramName 0.9797 > /dev/null && echo "meshing 00 ended") &
(gmsh -2 -o ${meshPath}${filename}01.msh -format msh4 ${filename}.geo -setnumber $paramName 0.5656 > /dev/null && echo "meshing 01 ended") &
(gmsh -2 -o ${meshPath}${filename}02.msh -format msh4 ${filename}.geo -setnumber $paramName 0.3703 > /dev/null && echo "meshing 02 ended") &
(gmsh -2 -o ${meshPath}${filename}03.msh -format msh4 ${filename}.geo -setnumber $paramName 0.3098 > /dev/null && echo "meshing 03 ended") &
(gmsh -2 -o ${meshPath}${filename}04.msh -format msh4 ${filename}.geo -setnumber $paramName 0.2190 > /dev/null && echo "meshing 04 ended") &
(gmsh -2 -o ${meshPath}${filename}05.msh -format msh4 ${filename}.geo -setnumber $paramName 0.1788 > /dev/null && echo "meshing 05 ended") &
(gmsh -2 -o ${meshPath}${filename}06.msh -format msh4 ${filename}.geo -setnumber $paramName 0.1549 > /dev/null && echo "meshing 06 ended") &
(gmsh -2 -o ${meshPath}${filename}07.msh -format msh4 ${filename}.geo -setnumber $paramName 0.1385 > /dev/null && echo "meshing 07 ended") &
(gmsh -2 -o ${meshPath}${filename}08.msh -format msh4 ${filename}.geo -setnumber $paramName 0.1171 > /dev/null && echo "meshing 08 ended") &
(gmsh -2 -o ${meshPath}${filename}09.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0979 > /dev/null && echo "meshing 09 ended") &
(gmsh -2 -o ${meshPath}${filename}10.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0692 > /dev/null && echo "meshing 10 ended") &
(gmsh -2 -o ${meshPath}${filename}11.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0565 > /dev/null && echo "meshing 11 ended") &
(gmsh -2 -o ${meshPath}${filename}12.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0438 > /dev/null && echo "meshing 12 ended") &
(gmsh -2 -o ${meshPath}${filename}13.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0370 > /dev/null && echo "meshing 13 ended") &
(gmsh -2 -o ${meshPath}${filename}14.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0309 > /dev/null && echo "meshing 14 ended") &
(gmsh -2 -o ${meshPath}${filename}15.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0178 > /dev/null && echo "meshing 15 ended") &
(gmsh -2 -o ${meshPath}${filename}16.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0138 > /dev/null && echo "meshing 16 ended") &
(gmsh -2 -o ${meshPath}${filename}17.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0117 > /dev/null && echo "meshing 17 ended") &
(gmsh -2 -o ${meshPath}${filename}18.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0097 > /dev/null && echo "meshing 18 ended") &
(gmsh -2 -o ${meshPath}${filename}19.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0069 > /dev/null && echo "meshing 19 ended") &
(gmsh -2 -o ${meshPath}${filename}20.msh -format msh4 ${filename}.geo -setnumber $paramName 0.0056 > /dev/null && echo "meshing 20 ended") &
(echo "")&
