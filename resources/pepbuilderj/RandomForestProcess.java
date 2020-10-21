package pepbuilderj;


import java.util.*;
import java.util.regex.*;
import java.io.*;




public class RandomForestProcess{
	ArrayList<RRFTree> trees;
	Pattern labelPat = Pattern.compile("^[$#][\\s]*([^\\s]+)");
	//Pattern lastNodePat = Pattern.compile("([^\\s=><]+)[\\s]*([><=]+)[\\s]*([^\\s=><:]+)[\\s]*:[\\s]*([^\\s=:;]+)[\\s]*");
	Pattern lastNodePat = Pattern.compile("([^\\s=><]+)[\\s]*([><=]+)[\\s]*([^\\s=><:]+)[\\s]*:[\\s]*([^\\s]+)[\\s]*");
	Pattern nodePat = Pattern.compile("([^\\s=><]+)[\\s]*([><=]+)[\\s]*([^\\s=><]+)[\\s]*");
	Hashtable<String,String> exLabel = new Hashtable<String,String>();//label for extra information
	
	/*shift_one was used for R result*/
	
	public void loadFromFile(String filename,boolean shift_one){
		ArrayList<String> al = new ArrayList<String>();
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(filename)));
			String line = "";
			while((line = br.readLine()) != null){
				al.add(line);
			}
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
			
		}
		loadFromString_Weka(al,shift_one);
		
	}
	public void loadFromString_Weka(ArrayList<String> pal,boolean shift_one){
		
		String currentlabel = "";
		Hashtable<String,StringBuffer> labeled = new Hashtable<String,StringBuffer>();
		
		for(int ii = 0;ii < pal.size();ii++){
			String str = pal.get(ii);
			Matcher mat = labelPat.matcher(str);
			if(mat.find()){
				if(currentlabel.equals("trees")){
					labeled.get(currentlabel).append(ii);
				}
				currentlabel = mat.group(1);
				
				if(currentlabel.equals("trees")){
					labeled.put(currentlabel,new StringBuffer());
					labeled.get(currentlabel).append(ii+1);
					labeled.get(currentlabel).append(",");
				}
			}else if(str.length() > 0){
				if(currentlabel.equals("trees")){
				}else{
					if(labeled.get(currentlabel) == null){
						labeled.put(currentlabel,new StringBuffer());
					}
					labeled.get(currentlabel).append(str+" \n");
				}
			}
		}
		
		if(currentlabel.equals("trees")){
			labeled.get(currentlabel).append(pal.size());
		}
		//pal.clear();
		Hashtable<String,Integer> labelmap = new Hashtable<String,Integer>();
		Pattern valpat = Pattern.compile("([^\\s=]+)[\\s=]+([^\\s=]+)");
		int placeshift = 0;
		Enumeration<String> el = labeled.keys();
		while(el.hasMoreElements()){
			String enn = el.nextElement();
			if(!enn.equals("place_shift") && !enn.equals("map") && !enn.equals("trees")){
				exLabel.put(enn,labeled.get(enn).toString());
			}
			
		}
		
		
		if(labeled.get("place_shift") != null){
			
			String lin = labeled.get("place_shift").toString().replaceAll("[\\s]","");
			placeshift = Integer.parseInt(lin);
		}
		if(labeled.get("map") != null){
			String[] lin = labeled.get("map").toString().split("[\r\n]+");;
			for(int ii = 0;ii < lin.length;ii++){
				Matcher vmat = valpat.matcher(lin[ii]);
				if(vmat.find()){
					labelmap.put(vmat.group(1),Integer.parseInt(vmat.group(2)));
					//System.out.println("mapping;;"+vmat.group(1)+"->"+vmat.group(2));
				}
				
			}
			
		}
		
		trees = new ArrayList<RRFTree>();
		if(labeled.get("trees") != null){
			ArrayList<String> buff = new ArrayList<String>();
			boolean flag = false;
			String[] n = labeled.get("trees").toString().split(",");
			int st = Integer.parseInt(n[0]);
			int en = Integer.parseInt(n[1]);
			for(int ii = st;ii < en;ii++){
				Matcher nm = nodePat.matcher(pal.get(ii));
				if(nm.find()){
					buff.add(pal.get(ii));
				}
				if(pal.get(ii).indexOf("==") == 0 || pal.get(ii).indexOf("--") == 0 ){
					if(buff.size()>0){
						trees.add(mapAllNodes_Weka(buff,labelmap,placeshift));
					}
					buff.clear();
				}
			}
			
			if(buff.size()>0){
				trees.add(mapAllNodes_Weka(buff,labelmap,placeshift));
			}
		}
	}
	
	
	public int countVerticalLines(String str){
		String[] pt = str.split("[\\s]+");
		int ret = 0;
		for(int ii = 0;ii < pt.length;ii++){
			if(pt[ii].equals("|")){
				ret++;
			}else if(pt[ii].replaceAll("[\\s]","").length() > 0){
				break;
			}
		}
		return ret;
	}
	
	public RRFTree mapAllNodes_Weka(ArrayList<String> str,Hashtable<String,Integer> labelmap,int pshift){
		/*
		RRFNode[] childNodes;
		RRFNode parentNode;
		int splitVarNum;
		boolean leafFlag;
		double splitThreshold;
		String predValue = null;
		boolean includeEquals_Up = true;
		*/
		RRFNode root = new RRFNode();
		ArrayList<String> al = new ArrayList<String>();
		for(int ii = 0;ii < str.size();ii++){
			Matcher pmat = nodePat.matcher(str.get(ii));
			if(pmat.find()){
				al.add(str.get(ii));
			//	System.out.println(al.get(ii));
			}
			
		}
		mapNodes(al,labelmap,root);
		ArrayList<RRFNode> nodes = getAllDescendants(root);
		RRFTree ret = new RRFTree();
		ret.nodes = new RRFNode[nodes.size()+1];
		ret.nodes[0] = root;
		ret.nodes[0].setSplitThreshold(ret.nodes[0].splitThreshold/Math.pow(10,pshift));
		for(int ii = 0;ii < nodes.size();ii++){
			ret.nodes[ii+1] = nodes.get(ii);
			ret.nodes[ii+1].setSplitThreshold(ret.nodes[ii+1].splitThreshold/Math.pow(10,pshift));
		}
		return ret;
	}
	
	public void mapNodes(ArrayList<String> al,Hashtable<String,Integer> labelmap,RRFNode parent){
		ArrayList<String> nextstr = new ArrayList<String>();
		if(al.size() == 0)return;
		al.add("");
		int vnum = countVerticalLines(al.get(0));
		//if(vnum == 0){
		//	return;
		//}
		
		RRFNode[] childtmp = new RRFNode[2];
		childtmp[0] = generateNodeFromString(al.get(0),labelmap,parent);
		int linnum = 1;
		ArrayList<String> clin = new ArrayList<String>();
		while(true){
			if(countVerticalLines(al.get(linnum)) <= vnum){
				if(clin.size() > 0){
					mapNodes(clin,labelmap,childtmp[0]);
				}
				childtmp[1] = generateNodeFromString(al.get(linnum),labelmap,parent);
				linnum++;
				break;
			}else{
				Matcher pmat = nodePat.matcher(al.get(linnum));
				if(pmat.find()){
					clin.add(al.get(linnum));
				}
			}
			linnum++;
		}
		
		clin = new ArrayList<String>();
		while(true){
			if(countVerticalLines(al.get(linnum)) <= vnum){
				if(clin.size() > 0){
					mapNodes(clin,labelmap,childtmp[1]);
				}
				clin.clear();
				break;
			}else{
				Matcher pmat = nodePat.matcher(al.get(linnum));
				if(pmat.find()){
					clin.add(al.get(linnum));
				}
			}
			linnum++;
		}
		if(al.get(0).indexOf("<")>-1){
			parent.setChildren(childtmp);
		}else if(al.get(0).indexOf(">")>-1){
			RRFNode[] rn = new RRFNode[2];
			rn[0] = childtmp[1];
			rn[1] = childtmp[0];
			parent.setChildren(rn);
		}else{
			System.out.println("Line cannot be parsed.\n"+al.get(0));
		}
		
		
		
		
	}
	
	public ArrayList<RRFNode> getAllDescendants(RRFNode rf){
		ArrayList<RRFNode> al = new ArrayList<RRFNode>();
		if(rf.childNodes != null){
			if(rf.childNodes[0] != null){
				al.add(rf.childNodes[0]);
				al.addAll(getAllDescendants(rf.childNodes[0]));
			}
			if(rf.childNodes[1] != null){
				al.add(rf.childNodes[1]);
				al.addAll(getAllDescendants(rf.childNodes[1]));
			}
		}
		return al;
		
	}
	
	
	public RRFNode generateNodeFromString(String str,Hashtable<String,Integer> labelmap,RRFNode parent){
		Matcher mat = nodePat.matcher(str);
		RRFNode ret = new RRFNode();
		if(mat.find()){
			Matcher lmat =lastNodePat.matcher(str);
			if(lmat.find()){
				//System.out.println(lmat.group(1));
				//System.out.println(labelmap.get(lmat.group(1)));
				//System.out.println(lmat.group(3));
				parent.setSplitVarNum(labelmap.get(lmat.group(1)));
				if(lmat.group(2).equals("<=")){
					parent.setIncludeEquals_Up(false);
				}
				parent.setSplitThreshold(Double.parseDouble(lmat.group(3)));
				ret.setPredValue(lmat.group(4));
				ret.asLeaf(true);
			}else{
				parent.setSplitVarNum(labelmap.get(mat.group(1)));
				if(mat.group(2).equals("<=")){
					parent.setIncludeEquals_Up(false);
				}
				parent.setSplitThreshold(Double.parseDouble(mat.group(3)));
			}
		}
		return ret;
		
	}
	/*
	public void setSplitVarNum(int v){
		splitVarNum = v;//when it is leaf, this should be -1
		if(v == -1){
			leafFlag = true;
		}else{
			leafFlag = false;
		}
		
	}
	public void setChildren(RRFNode[] c){
		childNodes = new RRFNode[c.length];
		for(int ii = 0;ii < c.length;ii++){
			childNodes[ii] = c[ii];
		}
	}
	public void setPredValue(String s){
		if(s == null)return;
		
		predValue = s;
		
		
	}
	
	public void setSplitThreshold(double d){
		splitThreshold = d;
	}*/
	
	public  void loadForestFromFile_SimpleFormat(String filename){
		try{
			loadForestFromStream_SimpleFormat(new FileInputStream(filename));
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	public  void loadForestFromStream_SimpleFormat(InputStream is){
		try{
			BufferedReader br = new BufferedReader(new InputStreamReader(is));
			String ln  = null;
			ArrayList<String> buff = new ArrayList<String>();
			trees = new ArrayList<RRFTree>();
			
			while((ln = br.readLine()) != null){
				if(ln.indexOf("#trees") > -1){
					break;
				}
			}
			while((ln = br.readLine()) != null){
				if(ln.indexOf("==") == 0){
					if(buff.size() > 0){
						RRFTree r = generateTreeFromString(buff);
						if(r != null){
							trees.add(r);
						}
					}
					buff.clear();
				}else if(ln.indexOf("#")>-1){
					break;
				}else{
					buff.add(ln);
				}
			}
			if(buff.size() > 0){
				RRFTree r = generateTreeFromString(buff);
				if(r != null){
					trees.add(r);
				}
			}
			
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
	}
	
	public RRFTree generateTreeFromString(ArrayList<String> al){
		Hashtable<Integer,RRFNode> map = new Hashtable<Integer,RRFNode>();
		//0	1	2	-0.16	-	1
		
		RRFTree ret = new RRFTree();
		int size = 0;
		for(int ii = 0;ii < al.size();ii++){
			String[] pt = al.get(ii).split("[\\s]+");
			if(pt.length > 4){
				RRFNode rf = new RRFNode();
				map.put(ii,rf);
				size++;
			}
		}
		if(size == 0){
			return null;
		}
		int index = 0;
		ret.nodes = new RRFNode[size];
		
		for(int ii = 0;ii < al.size();ii++){
			String[] pt = al.get(ii).split("[\\s]+");
			if(pt.length > 4){
				RRFNode rf = map.get(ii);
				ret.nodes[index] = rf;
				index++;
				if(pt[2].indexOf("-") < 0){
					rf.setChildren(map.get(Integer.parseInt(pt[1])),map.get(Integer.parseInt(pt[2])));
					rf.setSplitVarNum(Integer.parseInt(pt[3]));
					rf.setSplitThreshold(Double.parseDouble(pt[4]));
					if(Integer.parseInt(pt[6]) == 1){
						rf.setIncludeEquals_Up(true);
					}else{
						rf.setIncludeEquals_Up(false);
					}
					rf.asLeaf(false);
				}else{
					rf.setPredValue(pt[5]);
					rf.asLeaf(true);
				}
			}
		}
		return ret;
	}
	
	public void saveForestToFile_SimpleFormat(String filename){
		
		try{
			BufferedWriter bw = new BufferedWriter(new FileWriter(new File(filename)));
			Enumeration<String> en = exLabel.keys();
			while(en.hasMoreElements()){
				String ex = en.nextElement();
				bw.write("#"+ex+"\n");
				bw.write(exLabel.get(ex));
				bw.write("\n");
			}
			
			
			bw.write("#trees\n");
			for(int ii = 0;ii < trees.size();ii++){
				ArrayList<RRFNode> al = getAllDescendants(trees.get(ii).nodes[0]);
				al.add(0,trees.get(ii).nodes[0]);
				Hashtable<RRFNode,Integer> idmap = new Hashtable<RRFNode,Integer>();
				for(int jj = 0;jj < al.size();jj++){
					idmap.put(al.get(jj),jj);
				}
				for(int jj = 0;jj < al.size();jj++){
					RRFNode n = al.get(jj);
					bw.write(String.valueOf(jj)+"\t");
					if(n.isLeaf()){
						bw.write("-\t-\t");
						bw.write("-\t"+"-\t"+n.predValue+"\t");
						bw.write("0\n");
					}else{
						bw.write(String.valueOf(idmap.get(n.childNodes[0]))+"\t"+String.valueOf(idmap.get(n.childNodes[1]))+"\t");
						bw.write(String.valueOf(n.splitVarNum)+"\t"+String.valueOf(n.splitThreshold)+"\t-\t");
						if(n.includeEquals_Up){
							bw.write("1\n");
						}else{
							bw.write("0\n");
						}
					}
					/*
					RRFNode[] childNodes;
					RRFNode parentNode;
					int splitVarNum;
					boolean leafFlag;
					double splitThreshold;
					String predValue = null;
					boolean includeEquals_Up = true;
					*/
				}
				bw.write("===================\n");
				
			}
			bw.close();
		}catch(Exception exx){
			exx.printStackTrace();
		}
		
	}
	
	public String predict_Majority(double d[]){
		Hashtable<String,Integer> hash = new Hashtable<String,Integer>();
		int max = 0;
		String ret = "";
		
		for(int ii = 0;ii < trees.size();ii++){
			String v = trees.get(ii).predict(d);
			if(hash.get(v) == null){
				hash.put(v,0);
			}
			int current = hash.get(v)+1;
			hash.put(v,current);
			if(current > max){
				max = current;
				ret = v;
			}else if(current == max){
				ret = ret+";"+v;
			}
		}
		return ret;
		
		
	}
	
	
	
	/**
	 * returns HashMap object which stores class and probability against values posted
	 * 
	 * 
	 */
	public HashMap<String,Double> predictProbability(double d[]){
		HashMap<String,Integer> hash = new HashMap<String,Integer>();
		int max = 0;
		String ret = "";
		
		for(int ii = 0;ii < trees.size();ii++){
			String v = trees.get(ii).predict(d);
			if(hash.get(v) == null){
				hash.put(v,0);
			}
			int current = hash.get(v)+1;
			hash.put(v,current);
			if(current > max){
				max = current;
				ret = v;
			}else if(current == max){
				ret = ret+";"+v;
			}
		}
		HashMap<String,Double> res = new HashMap<String,Double>();
		Iterator<String> ite = hash.keySet().iterator();
		while(ite.hasNext()){
			String k = ite.next();
			double v = hash.get(k);
			res.put(k,v/trees.size());
		}
		
		return res;
	}
	
	
	public ArrayList<String> getResultList(HashMap<String,Double> hash,double threshold){
		ArrayList<Double> prob = new ArrayList<Double>();
		ArrayList<String> lab = new ArrayList<String>();
		
		
		Iterator<String> ite = hash.keySet().iterator();
		while(ite.hasNext()){
			String k = ite.next();
			double v = hash.get(k);
			if(v > threshold){
				prob.add(v);
				lab.add(k);
			}
		}
		Collections.sort(prob);
		Collections.reverse(prob);
		ArrayList<String> ret = new ArrayList<String>();
		
		for(Double p:prob){
			Iterator<String> iite = lab.iterator();
			while(iite.hasNext()){
				String ss = iite.next();
				if(hash.get(ss) >= p){
					ret.add(ss);
					iite.remove();
				}
				
			}
			
		}
		return ret;
	}
	public ArrayList<String> getResultList(HashMap<String,Double> hash,int num){
		ArrayList<Double> prob = new ArrayList<Double>();
		ArrayList<String> lab = new ArrayList<String>();
		
		
		Iterator<String> ite = hash.keySet().iterator();
		while(ite.hasNext()){
			String k = ite.next();
			double v = hash.get(k);
			prob.add(v);
			lab.add(k);
			
		}
		Collections.sort(prob);
		Collections.reverse(prob);
		ArrayList<String> ret = new ArrayList<String>();
		
		for(Double p:prob){
			Iterator<String> iite = lab.iterator();
			while(iite.hasNext()){
				String ss = iite.next();
				if(hash.get(ss) >= p){
					ret.add(ss);
					iite.remove();
				}
				
			}
			if(ret.size() >= num){
				break;
			}
		}
		return ret;
	}

	
	
	
	
	public double predict_Mean(double d[]){
		double sum = 0;
		for(int ii = 0;ii < trees.size();ii++){
			String v = trees.get(ii).predict(d);
			sum += Double.parseDouble(v);
		}
		return sum/trees.size();
	}
	
	
	public static HashMap<String,String> parseArg(String[] args){
		HashMap<String,String> ret = new HashMap<>();
		for(int ii = 0;ii < args.length;ii++){
			if(args[ii].indexOf("-") == 0){
				if(args.length > ii+1){
					ret.put(args[ii],args[ii+1]);
				}else{
					ret.put(args[ii],"");
				}
			}
		}
		return ret;
	}
	
	
	public static void normalPrediction(HashMap<String,String> arghash,RandomForestProcess rg){
		boolean header = arghash.containsKey("-hasheader");
		boolean hasnamecol = arghash.containsKey("-hasnamecol");
		boolean hasanswercol = arghash.containsKey("-hasanswercol");
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(arghash.get("-infile"))));
			BufferedWriter fw = new BufferedWriter(new FileWriter(new File(arghash.get("-outfile")),false));
			
			String line = null;
			if(header){
				line = br.readLine();
			}
			int okcount = 0;
			int allcount = 0;
			while((line = br.readLine()) != null){
				if(line.replaceAll("[\\s]","").length() == 0){
					continue;
				}
				
				ArrayList<String> al  = getParsedArray_String(line);
				
				int slen = al.size();
				
				String sname ="";
				String answer = "";
				if(hasnamecol){
					slen--;
					sname = al.remove(0);
				}
				if(hasanswercol){
					slen--;
					answer = al.remove(al.size()-1);
				}
				double[] d = new double[slen];
				for(int ii = 0;ii < al.size();ii++){
					d[ii] = Double.parseDouble(al.get(ii));
				}
				
				String predres = rg.predict_Majority(d);
				
				if(hasanswercol){
					if(predres.equals(answer)){
						okcount++;
					}
					allcount++;
				}
				
				
				if(hasanswercol){
					predres = "true:\t"+answer+"\tpred:\t"+predres;
				}
				if(hasnamecol){
					fw.write(sname+"\t"+predres+"\n");
				}else{
					fw.write(predres+"\n");
				}
			}
			if(hasanswercol){
				System.err.println(okcount+"/"+allcount);
			}
			fw.close();
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
			
		}
	}
	
	public static double average(ArrayList<Double> a){
		double sum = 0;
		for(Double d:a){
			sum += d;
		}
		return sum/a.size();
	}
	public static double var(ArrayList<Double> a){
		double av = average(a);
		double dsum = 0;
		for(Double d:a){
			dsum += (d-av)*(d-av);
		}
		return dsum/a.size();
	}
	public static double cor(ArrayList<Double> d1,ArrayList<Double> d2){
		double s1 = Math.sqrt(var(d1));
		double s2 = Math.sqrt(var(d2));
		double a1 = average(d1);
		double a2 = average(d2);
		double cv = 0;
		for(int ii = 0;ii < d1.size();ii++){
			cv += (d1.get(ii)-a1)*(d2.get(ii)-a2);
		}
		return cv/s1/s2/d1.size();
	}
	
	
	/**
	 * ノードを通るクラスの数を数える。
	 * normalPrediction と同じルーチンに入れたかったが、クラス数とか取らないといけなかったので中止
	 * @param arghash
	 * @param rg 
	 */
	public static void debugCount(HashMap<String,String> arghash,RandomForestProcess rg){
		boolean header = arghash.containsKey("-hasheader");
		boolean hasnamecol = arghash.containsKey("-hasnamecol");
		boolean hasanswercol = arghash.containsKey("-hasanswercol");
		if(!hasanswercol){
			throw new RuntimeException("Debug count function needs answer column.");
		}
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File(arghash.get("-infile"))));
			String line = null;
			if(header){
				line = br.readLine();
			}
			HashSet<String> cl = new HashSet<>();
			while((line = br.readLine()) != null){
				if(line.replaceAll("[\\s]","").length() == 0){
					continue;
				}
				ArrayList<String> al  = getParsedArray_String(line);
				cl.add(al.remove(al.size()-1));
			}
			br.close();
			
			br = new BufferedReader(new FileReader(new File(arghash.get("-testfile"))));
			if(header){
				line = br.readLine();
			}
			while((line = br.readLine()) != null){
				if(line.replaceAll("[\\s]","").length() == 0){
					continue;
				}
				ArrayList<String> al  = getParsedArray_String(line);
				cl.add(al.remove(al.size()-1));
			}
			br.close();

			ArrayList<String> classes = new ArrayList<>(cl);
			Collections.sort(classes);
			HashMap<Integer,String> index_label_map = new HashMap<>();
			HashMap<String,Integer> label_index_map = new HashMap<>();
			for(int ii = 0;ii < classes.size();ii++){
				index_label_map.put(ii,classes.get(ii));
				label_index_map.put(classes.get(ii),ii);
			}
			
			
			
			BufferedWriter fw = new BufferedWriter(new FileWriter(new File(arghash.get("-outfile")),false));
			String[] filenames = {arghash.get("-infile"),arghash.get("-infile")};
			ArrayList<ArrayList<Double>> firstresult = new ArrayList<>();
			for(int gg = 0;gg < 2;gg++){
				String fname = filenames[gg];
				br = new BufferedReader(new FileReader(new File(fname)));
				for(RRFTree t:rg.trees){
					t.prepareDebug(classes.size());
				}
				line = null;
				if(header){
					line = br.readLine();
				}
				while((line = br.readLine()) != null){
					if(line.replaceAll("[\\s]","").length() == 0){
						continue;
					}

					ArrayList<String> al  = getParsedArray_String(line);

					int slen = al.size();

					String sname ="";
					String answer = "";
					if(hasnamecol){
						slen--;
						sname = al.remove(0);
					}
					if(hasanswercol){
						slen--;
						answer = al.remove(al.size()-1);
					}
					double[] d = new double[slen];
					for(int ii = 0;ii < al.size();ii++){
						d[ii] = Double.parseDouble(al.get(ii));
					}

					for(RRFTree t:rg.trees){
						t.predict_debug(d,label_index_map.get(answer));
					}

				}

				for(int ii = 0;ii < rg.trees.size();ii++){
					RRFTree t = rg.trees.get(ii);
					
					t.countToRatio();
					for(RRFNode n:t.nodes){
						ArrayList<Double> dd = new ArrayList<>();
						for(double d:n.classCount){
							dd.add(d);
						}
						if(gg == 0){
							firstresult.add(dd);
						}else{
							fw.write(cor(firstresult.get(ii),dd)+"\n");
						}
					}
				}
				
				br.close();
			}
			//for(RRFTree t:rg.trees){
			//	t.printCountInfo(fw, index_label_map);
			//}
			fw.close();
		}catch(Exception exx){
			exx.printStackTrace();
			
		}
	}
	
	public static void main(String args[]){

		String[] aa = {
			"-treefile_simple","C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\3fold\\tree_simple.0.dat",
			"-infile","C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\input.0.train.dat"
			,"-testfile","C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\input.0.test.dat"
			,"-outfile","C:\\dummy\\vbox_share\\bioo\\database\\for_energyfunction\\14_makej48trees\\3fold\\ttout.dat"
			,"-countclasses",
			"-hasheader",
			"-hasnamecol",
			"-hasanswercol"
		};
		main__(aa);
	}
	public static void main__(String args[]){
		
		HashMap<String,String> arghash = parseArg(args);
		String[] argcheck ={
			"-treefile",
			"-treefile_simple",
			"-output_simple",
			"-infile",
			"-testfile",
			"-outfile",
			"-hasheader",
			"-hasnamecol",
			"-hasanswercol",
			"-countclasses"
		};
		HashSet<String> argchecker = new HashSet<>();
		for(String c:argcheck){
			argchecker.add(c);
		}
		if((!arghash.containsKey("-treefile") && !arghash.containsKey("-treefile_simple")) 
			|| (!arghash.containsKey("-outfile") && !arghash.containsKey("-output_simple")) ){
			System.err.println("usage: RandomForestProcess -treefile <treefile> (or -treefile_simple) -infile <tablefile> -outfile <result file> [-hasheader] [-hasnamecol] [-hasanswercols]");
			System.err.println("usage: RandomForestProcess -treefile <treefile> -output_simple <tree file in simple format>");
			System.exit(1);
		}
		
		for(String str:arghash.keySet()){
			if(!argchecker.contains(str)){
				throw new RuntimeException("Unknown option "+str);
			}
		}
		RandomForestProcess rg = new RandomForestProcess();
		if(arghash.containsKey("-treefile")){
			rg.loadFromFile(arghash.get("-treefile"),true);
		}
		if(arghash.containsKey("-treefile_simple")){
			rg.loadForestFromFile_SimpleFormat(arghash.get("-treefile_simple"));
		}
		if(arghash.containsKey("-output_simple")){
			rg.saveForestToFile_SimpleFormat(arghash.get("-output_simple"));
			System.exit(0);
		}
		if(arghash.containsKey("-countclasses")){
			debugCount(arghash,rg);
		}else{
			normalPrediction(arghash,rg);
		}
		//System.out.println(rg.trees.size()+":::");
		
		
	}
	
	
	public static void main_(String args[]){
		
		RandomForestProcess rg = new RandomForestProcess();
		rg.loadFromFile("testforest2.dat",true);
		rg.saveForestToFile_SimpleFormat("testforest_iris_simple.dat");
		rg.loadForestFromFile_SimpleFormat("testforest_iris_simple.dat");
		System.out.println(rg.trees.size()+":::");
		try{
			BufferedReader br = new BufferedReader(new FileReader(new File("iris_chk.dat")));
			String line = br.readLine();
			while((line = br.readLine()) != null){
				ArrayList<String> al  = getParsedArray_String(line);
				
				double[] d = {Double.parseDouble(al.get(0))
				,Double.parseDouble(al.get(1))
				,Double.parseDouble(al.get(2))
				,Double.parseDouble(al.get(3))};//iris
				System.out.println(rg.predict_Majority(d));
				
			}
			br.close();
		}catch(Exception exx){
			exx.printStackTrace();
			
		}
		
		
		
	}
	
	
	
	
	public static ArrayList<Integer> getParsedArray_Int(String str){
		String[] ar = str.split("[\\s]+");
		ArrayList<Integer> ret = new ArrayList<Integer>();
		for(int ii = 0;ii < ar.length;ii++){
			try{
				if(ar[ii].length() > 0){
					int pp = Integer.parseInt(ar[ii]);
					ret.add(pp);
				
				}
			}catch(Exception exx){
				exx.printStackTrace();
			}
		}
		return ret;
	}
	
	public static ArrayList<String> getParsedArray_String(String str){
		String[] ar = ("#"+str).replaceAll("[\r\n]","").split("[\\s]+");
		ar[0] = ar[0].replaceFirst("^#","");
		ArrayList<String> ret = new ArrayList<String>();
		for(int ii = 0;ii < ar.length;ii++){
			try{
				if(ar[ii].length() > 0){
					ret.add(ar[ii]);
				
				}
			}catch(Exception exx){
				exx.printStackTrace();
			}
		}
		return ret;
	}
	public static ArrayList<Double> getParsedArray_Double(String str){
		String[] ar = ("#"+str).replaceAll("[\r\n]","").split("[\\s]+");
		ar[0] = ar[0].replaceFirst("^#","");
		ArrayList<Double> ret = new ArrayList<Double>();
		for(int ii = 0;ii < ar.length;ii++){
			try{
				if(ar[ii].length() > 0){
					double pp = Double.parseDouble(ar[ii]);
					ret.add(pp);
				}
			}catch(Exception exx){
				exx.printStackTrace();
			}
		}
		return ret;
	}
}

class RRFTree{
	RRFNode[] nodes;
	RRFTree(){
	}
	
	public int getNodeNum(){
		return nodes.length;
	}
	public void setNodeNum(int num){
		nodes = new RRFNode[num];
		for(int ii = 0;ii < num;ii++){
			nodes[ii] = new RRFNode();
		}
	}
	public void setTreeMap(int[] tmap,boolean shift_one){
		for(int ii = 0;ii < tmap.length/2;ii++){
			int childnum_a = tmap[ii];
			int childnum_b = tmap[tmap.length/2+ii];
			if(shift_one){
				childnum_a--;
				childnum_b--;
			}
			if(ii < nodes.length){
				RRFNode[] cnode = new RRFNode[2];
				if(childnum_a > -1){
					cnode[0] = nodes[childnum_a];
				}else{
					cnode[0] = null;
				}
				if(childnum_b > -1){
					cnode[1] = nodes[childnum_b];
				}else{
					cnode[1] = null;
				}
				nodes[ii].setChildren(cnode);
			}
		}
				
	}
	
	public void setSplitVarNum(int sv[],boolean shift_one){
		for(int ii = 0;ii < nodes.length;ii++){
			if(shift_one){
				nodes[ii].setSplitVarNum(sv[ii]-1);
			}else{
				nodes[ii].setSplitVarNum(sv[ii]);
			}
		}
	}
	public void setSplitThreshold(double st[]){
		for(int ii = 0;ii < nodes.length;ii++){
			nodes[ii].setSplitThreshold(st[ii]);
		}
	}
	
	public void setPredValue(String st[]){
		for(int ii = 0;ii <  nodes.length;ii++){
			nodes[ii].setPredValue(st[ii]);
		}
	}
	public String predict(double[] d){
		RRFNode currentnode = nodes[0];
		while(true){
			if(currentnode.isLeaf()){
				return currentnode.getPredValue();
			}else{
				currentnode = currentnode.getNextNode(d);
			}
			if(currentnode == null){
				return null;
			}
		}
	}
	
	public void prepareDebug(int classnum){
		for(RRFNode n:nodes){
			n.prepareCount(classnum);
		}
	}
	public void countToRatio(){
		for(RRFNode n:nodes){
			n.countToRatio();
		}
	}
	public void printCountInfo(BufferedWriter bw,HashMap<Integer,String> code){
		HashMap<String,Integer> dcode = new HashMap<>();
		for(Integer i:code.keySet()){//インデクスとラベルが一対一に対応しているべき
			dcode.put(code.get(i),i);
		}
		ArrayList<String> labels = new ArrayList<>(dcode.keySet());
		Collections.sort(labels);
		try{
			for(RRFNode n:nodes){
				for(int ii = 0;ii < n.classCount.length;ii++){
					bw.write(labels.get(ii)+":"+((float)n.classCount[dcode.get(labels.get(ii))])+"\t");
				}
				bw.write("\n");
			}
		}catch(Exception exx){
			exx.printStackTrace();
		}
	}
	
	
	/**
	 * ノードを通ったクラスをカウントしていく
	 * @param d
	 * @param classcode
	 * @return 
	 */
	public String predict_debug(double[] d,int classcode){
		RRFNode currentnode = nodes[0];
		
		while(true){
			currentnode.classCount[classcode]++;
			if(currentnode.isLeaf()){
				return currentnode.getPredValue();
			}else{
				currentnode = currentnode.getNextNode(d);
			}
			if(currentnode == null){
				return null;
			}
		}
	}
	
	
}

class RRFNode{
	RRFNode[] childNodes;//threshold より小さい時 0 大きい時 1 に進む
	RRFNode parentNode;
	int splitVarNum;
	boolean leafFlag;
	double splitThreshold;
	String predValue = null;
	boolean includeEquals_Up = true;
	int depth = 0;
	double[] classCount;//そのノードに来たクラスの数を数える。デバッグ用
	
	
	RRFNode(){
		
	}
	
	public void prepareCount(int siz){
		classCount = new double[siz];
		for(int ii = 0;ii < siz;ii++){
			classCount[ii] = 0;
		}
	}
	
	public void countToRatio(){
		double sumcount = 0;
		for(int ii = 0;ii < classCount.length;ii++){
			sumcount += classCount[ii];
		}
		if(sumcount > 0){
			for(int ii = 0;ii < classCount.length;ii++){
				classCount[ii] /= sumcount;
			}
		}
	}
	
	
	public String getPredValue(){
		return predValue;
	}
	public boolean isLeaf(){
		return leafFlag;
	}
	public void asLeaf(boolean flag){
		leafFlag = flag;
	}
	public void setIncludeEquals_Up(boolean b){
		includeEquals_Up = b;
	}
	public RRFNode getNextNode(double[] v){
		if(includeEquals_Up){
			if(v[splitVarNum] >= splitThreshold){
				return childNodes[1];
			}else{
				return childNodes[0];	
			}
		}else{
			if(v[splitVarNum] > splitThreshold){
				return childNodes[1];
			}else{
				return childNodes[0];	
			}
		}
	}
	public void setSplitVarNum(int v){
		splitVarNum = v;//when it is leaf, this should be -1
		//System.out.println(v);
		if(v == -1){
			leafFlag = true;
		}else{
			leafFlag = false;
		}
		
	}
	public void setChildren(RRFNode[] c){
		childNodes = new RRFNode[c.length];
		for(int ii = 0;ii < c.length;ii++){
			childNodes[ii] = c[ii];
			childNodes[ii].depth = this.depth+1;
		}
	}
	public void setChildren(RRFNode c1,RRFNode c2){
		childNodes = new RRFNode[2];
		childNodes[0] = c1;
		childNodes[1] = c2;
		
	}
	public void setPredValue(String s){
		if(s == null)return;
		
		predValue = s;
		//if(splitVarNum != -1){
		//	System.out.println("Node must not have pred value."+splitVarNum+" "+s);
		//}
		
		//leafFlag = true;
		
		
	}
	
	public void setSplitThreshold(double d){
		splitThreshold = d;
	}
}



