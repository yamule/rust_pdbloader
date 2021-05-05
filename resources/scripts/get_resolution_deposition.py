import re;
import sys;
import os;


#https://www.wwpdb.org/documentation/file-format-content/format33/sect2.html
"""

Overview

The HEADER record uniquely identifies a PDB entry through the idCode field. This record also provides a classification for the entry. Finally, it contains the date when the coordinates were deposited to the PDB archive.

Record Format

COLUMNS       DATA  TYPE     FIELD             DEFINITION
------------------------------------------------------------------------------------
 1 -  6       Record name    "HEADER"
11 - 50       String(40)     classification    Classifies the molecule(s).
51 - 59       Date           depDate           Deposition date. This is the date the
                                               coordinates  were received at the PDB.
63 - 66       IDcode         idCode            This identifier is unique within the PDB.



とのことなので、DEP で判断する。
#2.5 未満、dep 2020 Jan 以前、残基数 50-600、非天然アミノ酸 5 %未満
#まず大雑把にフィルタリングする。DNA、RNA は除く。
#5CPJ に "STRUCTURAL_PROTEIN/DNA となぜかダブルクォーテーションが入っている。。。手作業で削除
#release date が DL 日以降で、とれていない構造があった。
"""

#with open("small.dat")as fin:
with open("meta_info_lines.dat")as fin:
	depdate = "";
	classcode = "";
	pdbcode = "";
	resolution = "";
	for ll in fin:
		ll = re.sub("[\s]*$","",ll);
		mat = re.search("^HEADER",ll);
		if mat:
			if len(pdbcode) != 0:
				raise Exception("??? {}\n{}".format(pdbcode,ll));
			classcode = re.sub("^[\s]+","",re.sub("[\s]+$","",ll[10:50]));
			depdate = re.sub("[\s]","",ll[50:59]);
			pdbcode = re.sub("[\s]","",ll[62:66]);
			
		
		mat = re.search("^REMARK +2 +RESOLUTION. +([^\s]+) +ANGSTROM",ll);
		if mat:
			resolution = mat.group(1);
		if "REMARK   2 RESOLUTION. NOT APPLICABLE." in ll:
			resolution = "-";
		
		if ll == "END":
			passflag = True;
			
			if classcode == "DNA" or classcode == "RNA":
				passflag = False;
			try:
				if resolution == "-":
					passflag = False;
				elif float(resolution) >= 2.5:
					passflag = False;
			except Exception as e:
				sys.stderr.write(str(e));
				sys.stderr.write(resolution+"??\n");
				sys.stderr.write("{}\t{}\t{}\t{}\t{}\n".format(pdbcode,classcode,depdate,resolution,passflag));
				passflag = "?";
				
			dd = re.split("-",depdate);
			if dd[-1] == "21" or dd[-1] == "20":
				passflag = False;
			else:
				try:
					chk = int(dd[-1]);
				except Exception as e:
					sys.stderr.write(str(e));
					sys.stderr.write(dd[-1]+"??\n");
					sys.stderr.write("{}\t{}\t{}\t{}\t{}\n".format(pdbcode,classcode,depdate,resolution,passflag));
					passflag = "?";
					
			classcode = re.sub("[\s]","_",classcode);
			print("{}\t{}\t{}\t{}\t{}".format(pdbcode,classcode,depdate,resolution,passflag));
			depdate = "";
			pdbcode = "";
			resolution = "";
			classcode = "";
		
		
		
		
		