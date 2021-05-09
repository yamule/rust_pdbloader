
# coding: utf-8

# In[24]:


import re;
import os;
import sys;
import math;
aa_3_1 = {
"ALA":"A",
"ARG":"R",
"ASN":"N",
"ASP":"D",
"CYS":"C",
"GLN":"Q",
"GLU":"E",
"GLY":"G",
"HIS":"H",
"ILE":"I",
"LEU":"L",
"LYS":"K",
"MET":"M",
"PHE":"F",
"PRO":"P",
"SER":"S",
"THR":"T",
"TRP":"W",
"TYR":"Y",
"VAL":"V"};

atom_check = {
"ALA":["C", "CA", "CB", "N", "O"],
"ARG":["C", "CA", "CB", "CD", "CG", "CZ", "N", "NE", "NH1", "NH2", "O"],
"ASN":["C", "CA", "CB", "CG", "N", "ND2", "O", "OD1"],
"ASP":["C", "CA", "CB", "CG", "N", "O", "OD1", "OD2"],
"CYS":["C", "CA", "CB", "N", "O", "SG"],
"GLN":["C", "CA", "CB", "CD", "CG", "N", "NE2", "O", "OE1"],
"GLU":["C", "CA", "CB", "CD", "CG", "N", "O", "OE1", "OE2"],
"GLY":["C", "CA", "N", "O"],
"HIS":["C", "CA", "CB", "CD2", "CE1", "CG", "N", "ND1", "NE2", "O"],
"ILE":["C", "CA", "CB", "CD1", "CG1", "CG2", "N", "O"],
"LEU":["C", "CA", "CB", "CD1", "CD2", "CG", "N", "O"],
"LYS":["C", "CA", "CB", "CD", "CE", "CG", "N", "NZ", "O"],
"MET":["C", "CA", "CB", "CE", "CG", "N", "O", "SD"],
"PHE":["C", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "CZ", "N", "O"],
"PRO":["C", "CA", "CB", "CD", "CG", "N", "O"],
"SER":["C", "CA", "CB", "N", "O", "OG"],
"THR":["C", "CA", "CB", "CG2", "N", "O", "OG1"],
"TRP":["C", "CA", "CB", "CD1", "CD2", "CE2", "CE3", "CG", "CH2", "CZ2", "CZ3", "N", "NE1", "O"],
"TYR":["C", "CA", "CB", "CD1", "CD2", "CE1", "CE2", "CG", "CZ", "N", "O", "OH"],
"VAL":["C", "CA", "CB", "CG1", "CG2", "N", "O"]
};



class PDBData:
	def __init__(self):
		self.chains = [];
	@staticmethod
	def load(infilename):
		alllines = [];
		if re.search("\.gz",infilename):
			import gzip
			with gzip.open(infilename, 'rb') as f:
				alllines = re.split("[\r\n]+",f.read().decode('utf-8'));
		else:
			with open(infilename) as fin:
				for ll in fin:
					alllines.append(ll);
		return PDBData.load_from_lines(alllines);
	@staticmethod
	def load_from_lines(alllines):
		ret = PDBData();
		atoms = [];
		ligands = [];
		chain_ended_all = False;
		chain_ended = {};
		chains = [];
		chains_hs = {};
		chains_hs_ligand = {};
		chains_missing = {};
		mresflag = False;
		matomflag = False;
		for ll in alllines:
		
			ll = re.sub("[\r\n]$","",ll);
			
			head = ll[0:6] ;
			if head == "ATOM  " or head == "HETATM":
				att = PDBAtom(ll);
				if not att.chain_id in chains_hs:
					chains.append(att.chain_id);
					chains_hs[att.chain_id] = [];
					chains_hs_ligand[att.chain_id] = [];

				if att.chain_id in chain_ended or chain_ended_all:
					chains_hs_ligand[att.chain_id].append(att);
				else:
					chains_hs[att.chain_id].append(att);
			if head == "TER   ":
				if len(ll) > 21 and ll[21]:
					chain_ended[ll[21]] = 100;
				else:
					chain_ended_all = True;
			if head == "ENDMDL":
				break;
		for cc in chains:
			if not cc in chains_missing:
				chains_missing[cc] = [];
			ret.chains.append(PDBChain(cc,chains_hs[cc],chains_hs_ligand[cc],chains_missing[cc]));
		return ret;
	def save(self,outfilename):
		with open(outfilename,"w") as fout:
			for cc in self.chains:
				for aa in cc.atoms:
					fout.write(aa.make_line());
					fout.write("\n");
				if len(cc.ligands) > 0:
					strr = cc.atoms[-1].make_line()[6:30];
					
					fout.write("TER   "+strr);
					fout.write("\n");
					for aa in cc.ligands:
						fout.write(aa.make_line());
						fout.write("\n");
	


# In[26]:


class PDBChain:
	def __init__(self,name,atoms,ligands,missingres):
		self.name = name;
		self.atoms = atoms;
		self.ligands = ligands;
		self.missing = missingres;
		
	def list_residues(self):
		res_printed = {};
		ret = [];
		for aa in self.atoms:
			rcode = aa.residue_name+"#"+str(aa.residue_pos)+"#"+aa.insertion_code;
			if not rcode in res_printed:
				ret.append(aa.residue_name);
				res_printed[rcode] = 100;
		return ret;
		
	def assemble_residues(self):
		res_printed = {};
		ret = [];
		for aa in self.atoms:
			rcode = aa.residue_name+"#"+str(aa.residue_pos)+"#"+aa.insertion_code;
			if not rcode in res_printed:
				ret.append(rcode);
				res_printed[rcode] = [];
			res_printed[rcode].append(aa);
		rret = [];
		for rcode in list(ret):
			rret.append(res_printed[rcode]);
		return rret;
		
	def get_fasta_seq(self):
		aalist = self.list_residues();
		ret = [];
		for aa in aalist:
			if aa in aa_3_1:
				ret.append(aa_3_1[aa]);
			else:
				ret.append("X");
		return "".join(ret);
	def get_atom_lines(self):
		ret = [];
		for aa in self.atoms:
			ret.push(aa.make_line());
		return ret;
class PDBAtom:
	def __init__(self,line):
		line += "                                    ";
		self.head = line[0:6];
		self.serial_number = int(line[6:11]);
		self.atom_name = re.sub(" ","",line[12:16]);
		self.alt_loc = line[16];
		self.residue_name = re.sub("[\s]","",line[17:20]);
		self.chain_id = line[21];
		self.residue_pos = int(line[22:26]);
		self.insertion_code = line[26];
		self.x = float(line[30:38]);
		self.y = float(line[38:46]);
		self.z = float(line[46:54]);
		self.occupancy = line[54:60];
		self.bfactor = line[60:66];
		self.element = line[76:78];
		self.charge = line[79:80];
	def get_atom_label(self):
		return self.chain_id+"#"+self.residue_name+"#"+str(self.residue_pos)+"#"+self.insertion_code+"#"+self.atom_name;
	def distance(self,target):
		xx = self.x-target.x;
		yy = self.y-target.y;
		zz = self.z-target.z;
		rr = xx*xx+yy*yy+zz*zz;
		if rr == 0.0:
			return 0.0;
		return math.sqrt(rr);
		
	
	def make_line(self):
		xx = "{:>.3f}".format(self.x);
		yy = "{:>.3f}".format(self.y);
		zz = "{:>.3f}".format(self.z);
		
		if len(xx) > 8:
			raise Exception("string overflow".format(self.x));
		if len(yy) > 8:
			raise Exception("string overflow".format(self.y));
		if len(zz) > 8:
			raise Exception("string overflow".format(self.z));
		atomname = self.atom_name;
		if self.head == "ATOM  " and len(self.atom_name) < 4:
			atomname = " "+atomname;
		ret = "{head}{serial_number:>5} {atom_name:<4}{alt_loc}{residue_name:>3} {chain_id}{residue_pos:>4}{insertion_code}   {xx:>8}{yy:>8}{zz:>8}{occupancy:>6}{bfactor:>6}		  {element:>2} {charge:2}".format(		
			head = self.head,
			serial_number = self.serial_number,
			atom_name = atomname,
			alt_loc = self.alt_loc,
			residue_name = self.residue_name,
			chain_id = self.chain_id,
			residue_pos = self.residue_pos,
			insertion_code = self.insertion_code,
			xx = xx,
			yy = yy,
			zz = zz,
			occupancy = self.occupancy,
			bfactor = self.bfactor,
			element = self.element,
			charge = self.charge);
		return ret;

pdbb = PDBData.load(sys.argv[1]);
for cc in list(pdbb.chains):
	alen = len(cc.atoms);
	multi_chain = False;
	hetflag = False;
	abnormal_aa = False;
	chain_break = False;
	has_middle_missing_atom = False;
	has_terminal_missing_atom = False;
	has_alt = False;
	nonmissing_aa = 0;
	aas = [];
	for aa in list(cc.atoms):
		if "HETATM" in aa.head:
			hetflag = aa.residue_name;
		if aa.alt_loc != "" and aa.alt_loc != " ":
			has_alt = True;
	rres = cc.assemble_residues();
	pplen = len(rres);
	for (ii,rr) in enumerate(list(rres)):
		has_missing_atom = False;
		if not rr[0].residue_name in aa_3_1:
			abnormal_aa = rr[0].residue_name;
		else:
			coveredatom = {};
			for aa in list(rr):
				a = aa.atom_name;
				if a == "OXT":
					a = "O";
				#if a not in atom_check[rr[0].residue_name]:
				#	sys.stderr.write(a+" "+rr[0].residue_name+"\n");
				coveredatom[a] = 100;
			for aa in list(atom_check[rr[0].residue_name]):
				if aa not in coveredatom:
					has_missing_atom = True;
					if ii == 0 or ii == pplen -1:
						has_terminal_missing_atom = rr[0].residue_name+":"+aa;
					else:
						has_middle_missing_atom = rr[0].residue_name+":"+aa;
		
		if not has_missing_atom:				
			if not rr[0].residue_name in aa_3_1:
				aas.append("X");
			else:
				aas.append(aa_3_1[rr[0].residue_name]);
	
	for aa in list(cc.ligands):
		if aa.residue_name == "HOH":
			continue;
		if "HETATM" in aa.head:
			hetflag = aa.residue_name;
			
	if len(pdbb.chains) > 1:
		names = [];
		for pp in list(pdbb.chains):
			names.append(pp.name);
		multi_chain = ";".join(names);
		
	for ii in range(alen-1):
		cdiff = cc.atoms[ii+1].residue_pos - cc.atoms[ii].residue_pos;
		if cdiff != 0 and cdiff != 1:
			chain_break = str(cc.atoms[ii].residue_pos)+"->"+str(cc.atoms[ii+1].residue_pos);
			break;
	nonmissing_aa = len(aas);
	aminoacids = "".join(aas);
	res = "multi_chain:{}\thetflag:{}\tabnormal_aa:{}\tchain_break:{}\thas_middle_missing_atom:{}\thas_terminal_missing_atom:{}\thas_alt:{}\tnonmissing_aa:{}\taminoacids:{}".format(
	       multi_chain,    hetflag,    abnormal_aa,    chain_break,    has_middle_missing_atom,    has_terminal_missing_atom,    has_alt,    nonmissing_aa,    aminoacids)
	print(sys.argv[1]+"\t"+cc.name+"\t"+res);
	