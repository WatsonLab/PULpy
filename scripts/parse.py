#!/usr/bin/env python

from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from reportlab.lib.units import cm

import csv
import re
import sys


#gd_diagram = GenomeDiagram.Diagram("test", x=0.01, y=0.01, yt=0.001, yb=0.001, track_size=0.9)

width = int(sys.argv[3])
x = float(50.0/float(width))
gd_diagram = GenomeDiagram.Diagram("test",  track_size=0.9, x=x)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", scale_ticks=False, scale_largeticks=0)
gd_feature_set = gd_track_for_features.new_set()

# open PUL data
pul_names = {}
pul_file = open(sys.argv[1], mode="r")
pul_data = csv.reader(pul_file, delimiter="\t")
for row in pul_data:
	pul_names[row[3]] = row[4]

fmin = 1000000
fmax = 0

my_file = open(sys.argv[2], mode='r')
parsed_data = csv.reader(my_file, delimiter="\t")
for row in parsed_data:
	con  = row[0]
	star = int(row[3])
	end  = int(row[4])
	stra = row[6]
	gene = row[8]

	prot = con + "_" + gene.split(';')[0].split('_')[1]

	#print prot
	if prot in pul_names.keys():

		if star < fmin:
			fmin = int(star)

		if end > fmax:
			fmax = int(end)

		label_angle = 0
		strand = "+1"
		if stra == "-":
			label_angle=180
			strand = "-1"

		print int(star)
		print int(end)
		print int(strand)
		print pul_names[prot]

		defcolor = colors.grey

		if pul_names[prot] == "susC":
			defcolor = colors.purple

		if pul_names[prot] == "susD":
                        defcolor = colors.orange

		match = re.search("GH", pul_names[prot])
		if match:
			defcolor = colors.pink

		match = re.search("CBM", pul_names[prot])
                if match:
                        defcolor = colors.pink

		match = re.search("PL", pul_names[prot])
                if match:
                        defcolor = colors.pink

		feature = SeqFeature(FeatureLocation(start=int(star), end=int(end), strand=int(strand)))
		gd_feature_set.add_feature(feature, color=defcolor, name=pul_names[prot], sigil="ARROW",
                                   label=True, label_size = 14,
                                   label_color=colors.blue, label_angle=label_angle)
	else:
		continue
	


print fmin
print fmax

pad = int(int(fmax-fmin+1) / 10)
pad = 100

width = int(sys.argv[3])

gd_diagram.draw(format="linear", pagesize=(width,150), fragments=1,
                start=int(fmin)-pad, end=int(fmax)+pad)

gd_diagram.write(sys.argv[4] + "linear.pdf", "PDF")
gd_diagram.write(sys.argv[4] + "linear.jpg", "JPG")
#gd_diagram.write("plasmid_linear_nice.eps", "EPS")
#gd_diagram.write("plasmid_linear_nice.svg", "SVG")
#gd_diagram.write("plasmid_linear_nice.png", "PNG")

