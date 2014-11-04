package supportClasses;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

public class BlastFile {

	File file;
	boolean xml;
	Hashtable<String,BlastIteration> iterations;
	
	public BlastFile (File file)throws IOException{
		this.file = file;
		BufferedReader in = new BufferedReader(new FileReader(file));
		String inputline = in.readLine();
		if(inputline.trim().startsWith("<?xml")){
			xml=true;
		}else{
			xml = false;
		}
		in.close();
		this.loadCompleteFile();
		
	}
	
	
	
	
	public void writeBlastTSVFile(File outputFile)throws IOException{
		
		Vector<String> v = new Vector<String> ();
		for(Enumeration<String> myenum = iterations.keys(); myenum.hasMoreElements();){
			v.add(myenum.nextElement());
		}
		Collections.sort(v);
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write(BlastIteration.getHeaderForIPKParserTSV());
		out.newLine();
		
		for(Enumeration<String> myenum = v.elements(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			BlastIteration it = iterations.get(key);
			out.write(it.getIPKParserTSV());
		}
		
		
		out.close();
	}
	
	
	
	
	
	public void writeFileSwitchQueryAndHit(File outputFile)throws IOException{
		Hashtable<String,Vector<BlastHit>> h = new Hashtable<String,Vector<BlastHit>>();
		
		for(Enumeration<BlastIteration> myenum1 = iterations.elements(); myenum1.hasMoreElements();){
			BlastIteration it = myenum1.nextElement();
			for(Enumeration<BlastHit> myenum2 = it.getHits().elements(); myenum2.hasMoreElements();){
				BlastHit hit = myenum2.nextElement();
				BlastHit revHit = switchHit(hit);
				
				String key = revHit.getQueryID();
				System.out.println(hit.getQueryID() + " " + key);
				
				Vector<BlastHit> v = new Vector<BlastHit>();
				if(h.containsKey(key)){
					v= h.get(key);
				}
				v.add(revHit);
				h.put(key,v );
			}
			
			
		}
		
		
		
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write(BlastIteration.getHeaderForIPKParserTSV());
		out.newLine();
		for(Enumeration<Vector<BlastHit>> myenum1 = h.elements(); myenum1.hasMoreElements();){
			Vector<BlastHit> v = myenum1.nextElement();
			Collections.sort(v);
			for(Enumeration<BlastHit> myenum2 = v.elements(); myenum2.hasMoreElements();){
				BlastHit b = myenum2.nextElement();
				out.write(b.getIPKParserTSV());
				
			}
			
			
			
		}
		
		
		
		
		out.close();
		
		
		
		
		
	}
	
	
	private static BlastHit switchHit (BlastHit inputHit){
		BlastHit outputHit = new BlastHit(inputHit.getHitID(), inputHit.getQueryID());
		for(Enumeration<BlastHSP> myenum = inputHit.getHSPs().elements(); myenum.hasMoreElements();){
			BlastHSP i = myenum.nextElement();
			BlastHSP outputHSP = new BlastHSP(i.getHitName(), i.getQueryName(), i.getHitDescription(), i.getQueryDescription(),
					i.getHspLength(), i.getHitLength(), i.getQueryLength(),
					i.getNumIdentical(), i.getNumberOfMismatches(), i.getNumberOfGaps(), i.getNumberOfGapOpenings(), i.getNumPositives(),
					i.getHitStart(),i.getHitEnd(),i.getQueryStart(),i.getQueryEnd(),
					i.getRawScore(),i.getBitScore(),i.getEvalue(),i.getHitString(),i.getQueryString(),i.getHitFrame(),i.getQueryFrame(),i.getHitStrand(),i.getQueryStrand());
			
			outputHSP.setAlignmentString(i.getAlignmentString());
			outputHSP.setHspRank(i.getHspRank());
			outputHSP.setNumHSPs(i.getNumHSPs());
			
			
			
			
			/*
			 String queryName, String hitName, String queryDescription, String hitDescription, 
			        int hspLength, int queryLength, int hitLength, 
			        int numIdentical, int numMisMatches, int numGaps, int numGapOpenings, int numPositives,
			        int queryStart, int queryEnd, int hitStart, int hitEnd,
			        double rawScore, double bitScore, String evalue, String queryString, String hitString,
			        int queryFrame, int hitFrame, int queryStrand, int hitStrand 
			  
			 
			 */
			
			outputHit.addHSP(outputHSP);
		}
		
		
		
		
		return outputHit;
	}
	
	
	
	private void loadCompleteFile()throws IOException{
		iterations = new Hashtable<String,BlastIteration> ();
		if(xml){   //xml version
			BlastXMLReader reader = new BlastXMLReader(file);
			BlastIteration it = reader.readIteration();
			while(it != null){
				iterations.put(it.getQueryID(), it);
				it = reader.readIteration();
			}
			reader.close();
		}
		else{    //Tab separated version
			BlastIPKReader reader = new BlastIPKReader(file);
			BlastIteration it = reader.readIteration();
			while(it != null){
				iterations.put(it.getQueryID(), it);				
				it = reader.readIteration();
			}
			reader.close();
		}
		
		
	}
	
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			File f = new File("/Users/steuernb/Documents/projects/test/test.blastn.xml");
			
			BlastFile file = new BlastFile(f);
			file.writeBlastTSVFile(new File(f.getParentFile(), "test.parsed"));
			file.writeFileSwitchQueryAndHit(new File(f.getParentFile(), "test_switched.parsed"));
		} catch (Exception e) {
			e.printStackTrace();
		}
		

	}

	
	
	
	
	
	
	
	
	
}
