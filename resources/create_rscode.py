

#.dic を適当にパースして、
# mandatory が no の場合は Some で入れる。



import re;
import os;
import sys;

with open("mmcif_pdbx_v50.dic") as fin:
	current_label = "";
	for ll in fin:
		mat = re.search("^save_(_[^\s]+)",ll);
		if mat:
			if len(current_label) > 0:
				print(current_label);
			current_label = mat.group(1);
		
		mat = re.search("^ +_item.mandatory_code +([^\s]+)",ll);
		if mat:
			gcode = mat.group(1);
			if gcode == "no":
				print("pub "+current_label+":"+"Some(String),");
			elif gcode == "yes":
				print("pub "+current_label+":String,");
			else:
				raise Exception(gcode+" is not handled.");
			current_label = "";
				