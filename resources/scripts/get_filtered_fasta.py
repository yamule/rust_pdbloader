import re;
import sys;
import os;

#get_resolution_deposition.py での出力を resotable.dat として保存しておく
#pdb_prop_check.py での出力を pdb_prop.dat として保存しておく

passed = {};
monomerflag = False;
resotablefile = "D:\\dummy\\vbox_share\\bioo\\database\\resotable.dat";
pdbpropfile = "results/pdb_prop.dat";
with open(resotablefile) as fin:
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
	ret["file"] =  [ptt.pop(0)];
	ret["chain"] =  [ptt.pop(0)];
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
	
ngchain = {};#Chain のうち一つでも NG である場合
with open(pdbpropfile) as fin:
	for ll in fin:
		pt_hs = parse_line(ll);
		
		if monomerflag and pt_hs["multi_chain"][0] != "False":
			ngchain[pt_hs["file"][0]] = 100;
			continue;
			
		if pt_hs["hetflag"][0] == "False" and \
		pt_hs["hetflag_allchain"][0] == "False" and \
		pt_hs["abnormal_aa"][0] == "False" and \
		pt_hs["chain_break"][0] == "False" and \
		pt_hs["has_middle_missing_atom"][0] == "False" and \
		pt_hs["alt_flag"][0] == "False" and \
		len(pt_hs["aminoacids"][0]) >= 50 and \
		len(pt_hs["aminoacids"][0]) <= 600 :
			pdbid = "";
			mat = re.search("pdb(....)\.ent",pt_hs["file"][0]);
			if mat:
				pdbid = mat.group(1).upper();
			else:
				mat = re.search("([0-9]...)\.pdb",pt_hs["file"][0]);
				if mat:
					pdbid = mat.group(1).upper();
			if pdbid == "":
				raise Exception("Cannot parse pdbid. {}\n".format(pt_hs["file"][0]));
			if not pdbid in passed:
				ngchain[pt_hs["file"][0]] = 100;
		else:
			ngchain[pt_hs["file"][0]] = 100;

with open(pdbpropfile) as fin:
	for ll in fin:
		pt_hs = parse_line(ll);
		
		if pt_hs["file"][0] in ngchain:
			continue;
		mat = re.search("pdb(....)\.ent",pt_hs["file"][0]);
		if mat:
			pdbid = mat.group(1).upper();
		else:
			mat = re.search("([0-9]...)\.pdb",pt_hs["file"][0]);
			if mat:
				pdbid = mat.group(1).upper();
		if re.search("[a-z]",pt_hs["chain"][0]):
			pt_hs["chain"][0] = "lower_"+pt_hs["chain"][0];
		print(">"+pdbid+"_"+pt_hs["chain"][0]+"\n"+pt_hs["aminoacids"][0]+"\n");
