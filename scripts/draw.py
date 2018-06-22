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

width = int(sys.argv[3])
x = float(50.0/float(width))
gd_diagram = GenomeDiagram.Diagram("test",  track_size=0.9, x=x)
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features", scale_ticks=False, scale_largeticks=0)
gd_feature_set = gd_track_for_features.new_set()

pulid = sys.argv[2]

# open PUL data
fmin = 1000000000
fmax = 0

my_file = open(sys.argv[1], mode='r')
parsed_data = csv.reader(my_file, delimiter="\t")
for row in parsed_data:
	con  = row[3]
	pul  = row[1]

	if pul == "pulid":
		continue

	star = int(row[4])
	end  = int(row[5])
	stra = row[6]
	sus  = row[9]
	gene = row[10]

	if sus == "":
		sus = gene

	if not pul == pulid:
		continue

	#print prot

	print(star)
	print(end)
	print(sus)

	if star < fmin:
		fmin = int(star)

	if end > fmax:
		fmax = int(end)

	label_angle = 0
	strand = "+1"
	if stra == "-":
		label_angle=180
		strand = "-1"

	defcolor = colors.grey

	if sus == "susC":
		defcolor = colors.purple

	if sus == "susD":
                defcolor = colors.orange

	match = re.search("GH", sus)
	if match:
		defcolor = colors.pink

	match = re.search("CBM", sus)
	if match:
		defcolor = colors.pink

	match = re.search("PL", sus)
	if match:
		defcolor = colors.pink

	feature = SeqFeature(FeatureLocation(start=int(star), end=int(end), strand=int(strand)))
	gd_feature_set.add_feature(feature, color=defcolor, name=sus, sigil="ARROW",
		label=True, label_size = 14,
                label_color=colors.blue, label_angle=label_angle)
	


print(fmin)
print(fmax)


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

