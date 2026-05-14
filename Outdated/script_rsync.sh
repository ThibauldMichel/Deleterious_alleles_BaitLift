#!/bin/bash

rsync -azhv  hyphy_iqtree/*  hyphy_FastTreeAstral/  --exclude gene_trees 

