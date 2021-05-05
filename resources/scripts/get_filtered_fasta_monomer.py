import re;
import sys;
import os;

#get_resolution_deposition.py での出力を resotable.dat として保存しておく
#pdb_prop_check.py での出力を pdb_prop.dat として保存しておく

passed = {};
with open("D:\\dummy\\vbox_share\\bioo\\database\\resotable.dat") as fin:
	for ll in fin:
		ptt = re.split("\t",re.sub("[\r\n]","",ll));
		if ptt[-1] == "True":
			passed[ptt[0].upper()] = 100;
		else:
			if ptt[-1] != "False":
				raise Exception("???");

def parse_line(ll):
	ret = {};
	ptt = re.split("\t",re.sub("[\r\n]","",ll));
	for pp in ptt:
		mat = re.search("([^\\:]+)\:(.+)",pp);
		if mat:
			k = mat.group(1);
			v = mat.group(2);
		else:
			k = "?";
			v = pp;
		if not k in ret:
			ret[k] = [];
		ret[k].append(v);
	return ret;
passed2 = {};
with open("pdb_prop.dat") as fin:
	for ll in fin:
		pt_hs = parse_line(ll);
		
		if pt_hs["multi_chain"][0] == "False" and \
		pt_hs["hetflag"][0] == "False" and \
		pt_hs["abnormal_aa"][0] == "False" and \
		pt_hs["chain_break"][0] == "False" and \
		pt_hs["has_missing"][0] == "False":
			pdbid = "";
			mat = re.search("pdb(....)\.ent",pt_hs["?"][0]);
			if mat:
				pdbid = mat.group(1).upper();
			mat = re.search("([0-9]...)\.pdb",pt_hs["?"][0]);
			if mat:
				pdbid = mat.group(1).upper();
			if pdbid == "":
				raise Exception("Cannot parse pdbid. {}\n".format(ptt[0]));
			if len(pt_hs["aa"][0]) < 50 or len(pt_hs["aa"][0]) > 600:
				continue;
			if pdbid in passed:
				print(">"+pdbid+"_"+pt_hs["?"][1]+"\n"+pt_hs["aa"][0]+"\n");
				
				