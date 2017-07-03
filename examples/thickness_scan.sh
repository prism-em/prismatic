#! /usr/bin/env bash

base_output_name="thickness_scan"
base_output_ext=".mrc"
atom_filename="../SI100.XYZ"
tileX=3
tileY=3

for tileZ in {1..17}
do
	prismatic -i ${atom_filename} -t ${tileX} ${tileY} ${tileZ} -o base_output_name${tileZ}base_output_ext
done