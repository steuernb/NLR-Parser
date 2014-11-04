package coreClasses;


import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Collections;
import java.util.Comparator;
import java.util.Date;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.Vector;

import supportClasses.BlastIteration;
import supportClasses.BlastReader;
import supportClasses.BlastXMLReader;
import supportClasses.FastaReader;
import supportClasses.FastaSequence;
import supportClasses.OutOfFrameException;

public class MastFile {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		try {
			MastFile file = new MastFile(new File("/Users/steuernb/Documents/projects/MastParser/mast.xml"));
			file.writeOutput(new File("/Users/steuernb/Documents/projects/MastParser/test.mast.txt"));
			
		} catch (Exception e) {
			e.printStackTrace();
		}

	}

	
	
	File mastFile;
	Hashtable<String, MastEntry> entries;
	Hashtable<String, MastMotifHitList > nblrrs;
	
	double pvalue_threshold;
	String split_pattern;
	
	Hashtable<String,Vector<int[]>> blast;
	
	
	public MastFile(File mastFile)throws IOException{
		this.mastFile         = mastFile;
		this.entries          = new Hashtable<String, MastEntry>();
		this.split_pattern    = "frame_";  //default split pattern
		this.pvalue_threshold = 1e-5;      //default value
		this.readEntries();
	}
	
	public MastFile(File mastFile, String split_pattern, double pvalue_threshold)throws IOException{
		this.split_pattern    = split_pattern;
		this.pvalue_threshold = pvalue_threshold;
		this.mastFile         = mastFile;
		this.entries          = new Hashtable<String, MastEntry>();
		this.readEntries();
	}
	
	private void readEntries()throws IOException{
		MastXMLReader reader = new MastXMLReader(this.mastFile);
		MastEntry entry = reader.readEntry();
		while(entry != null){
			entries.put(entry.getName(), entry);
			entry = reader.readEntry();
		}
		reader.close();
	}
	
	
	/**
	 * Merge found motifs from different sequences according to a common prefix. The prefix is determined 
	 * by the split_pattern field. All sequences that have an equal beginning of the identifier until the first occurence 
	 * of the split_pattern. This is useful if a genomic locus was 6-frame-translated. 
	 * Motifs on different exons are potentially on a different reading frame.  
	 * 
	 * @throws IOException
	 */
	public void mergeEntries()throws IOException{
		nblrrs = new Hashtable<String,MastMotifHitList>();
		if(entries.size()==0){
			this.readEntries();
		}
		
		for(Enumeration<String> myenum = entries.keys(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			String name = key.split(split_pattern)[0];
			Vector<MastMotifHit> v = entries.get(key).getHits();
			MastMotifHitList list = new MastMotifHitList();
			list.setLength(entries.get(key).getLength());
			
			if(nblrrs.containsKey(name)){
				list= nblrrs.get(name);
			}
			for(Enumeration<MastMotifHit> myenum2 = v.elements(); myenum2.hasMoreElements();){
				MastMotifHit mmh = myenum2.nextElement();
				if(mmh.getPvalue()<=this.pvalue_threshold){
					list.addMotif(mmh);
				}
			}
			if(list.getSize()>0){
				nblrrs.put(name, list);
			}
			
		}
		
	}
	

	/**
	 * 
	 * Add a blast file. Currently the file is only used to determine a potential protein start. 
	 * The blast file has to be in xml format (for ncbi blast use -outfmt 5) and the query-ids of the blast file
	 * should match the sequence ids or prefixes (see split pattern) used in the mast file.  
	 * 
	 * @param blastXML
	 * @throws IOException
	 */
	public void setBlast(File blastXML)throws IOException{
		blast = new Hashtable<String,Vector<int[]>>();
		
		BlastReader reader = new BlastXMLReader(blastXML);
		
		for( BlastIteration it = reader.readIteration(); it != null; it = reader.readIteration()){
			
			it.removeHSP(80);
			
			
			int[] a = it.getQueryCoverage();
			int start = -1;
			Vector<int[]> v = new Vector<int[]>();
			for( int i = 0; i< a.length; i++){
				
				if(a[i] > 0){
					if( start == -1){
						start = i+1;
					}
				}
				
				if(a[i] == 0){
					if(start>0){
						int[] aa = {start, i};
						
						v.add(aa);
						start = -1;
					}
				}
			}
			if(start != -1){
				int[] aa = {start, a.length };
				v.add(aa);
			}
			if( v.size()>0){
				blast.put(it.getQueryID(), v);
			}	
		}
		
	}
	
	
	/**
	 * Write a tab separated output file. 
	 * 
	 * The columns are: 
	 * SequenceName	class	complete	start	end	NB-ARC	2-NBLRR-Signal	MotifList
	 * 
	 * SequenceName:	The (merged) name of the sequence
	 * class:			The class of the NB-LRR gene. CNL or TNL
	 * complete:		Either "complete" or "partial" depending if all essential motifs have been found.
	 * start:			A potential start of the protein. If a blast file was set it is used. Othewise it is the beginning of the first motif.
	 * end:				A potential end of the protein. This is the first stop codon or undefined amino acid after the last motif.
	 * NB-ARC			If the motif 1 is present this gives 210 amino acids starting 10 before the start of motif 1. Note that this is always the same reading frame as the motif 1. So this might be misleading.
	 * 2-NBLRR-Signal	If we get a motif1 after and LRR motif (9, 11 oder 19), the sequence might be more than one NBLRR. In that case this will be set to true.
	 * MotifList		The sequence of motifs found for this sequence.
	 * 
	 * @param outputFile
	 * @throws IOException
	 */
	public void writeOutput(File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		out.write("SequenceName\tclass\tcomplete\tstart\tend\tNB-ARC\t2-NBLRR-Signal\tMotifList");
		out.newLine();
		for(Enumeration<String> myenum = this.getNBLRRs().keys(); myenum.hasMoreElements();){
			String sequenceName = myenum.nextElement();
			MastMotifHitList list = this.getNBLRRs().get(sequenceName);
			
			String result = list.checkMotifCombination();
			if(!result.equalsIgnoreCase("FALSE")){
				String nblrr_class = list.getNBLRR_Class();
				boolean is_complete = list.isComplete();
				String complete = "complete";
				if(!is_complete){
					complete = "partial";
				}
				
				String start = list.getStart();
				String end = list.getEnd();
				String aa_nbarc = list.getNBARC_Sequence();
				boolean hasTwo = list.hasPotentiallyTwoNBLRRs();
				String motifs = list.getMotifListString();
				
				if(blast!= null){
					if(blast.containsKey(sequenceName)){
						Vector<int[]> v = blast.get(sequenceName);
						Collections.sort(v, new Comparator<int[]>(){
							public int compare(int[] a1, int[] a2){
								if(a1[0] ==a2[0]){
									return 0;
								}else if( a1[0] <=a2[0]){
									return -1;
								}else {
									return 1;
								}
								
							}
						});
						
						int motivpos = Integer.parseInt(start.split(":")[1]);
						int blastpos = v.get(0)[0];
						if(blastpos > motivpos){
							start = blastpos+"";
						}
						
						
						
					}
					
				}
				
				out.write(sequenceName+"\t"+nblrr_class+"\t"+complete+"\t"+start+"\t"+end+"\t"+aa_nbarc+"\t"+hasTwo+"\t"+motifs);
				out.newLine();
			}
		}
		
		
		
		out.close();
	}
	
	
	public void addSequences(File sequenceFile)throws IOException{
		
		FastaReader fastaReader = new FastaReader(sequenceFile);
		for (FastaSequence seq = fastaReader.readEntry(); seq != null; seq = fastaReader.readEntry()) {
			if( this.entries.containsKey(seq.getIdentifier())){
				MastEntry entry = this.entries.get(seq.getIdentifier());
				for(Enumeration<MastMotifHit> myenum = entry.getHits().elements(); myenum.hasMoreElements();){
					myenum.nextElement().setAaSequence(seq.getSequence());
				}
				
			}

		}
		fastaReader.close();
		
		
		
	}
	
	
	/*
	
	
	public static void shortOutput(Hashtable<String, Vector<MastMotifHit>> h , File outputFile, File blastFile, String splitPattern)throws IOException{
		Hashtable<String,String[]> starts = getStarts( blastFile, splitPattern);
		
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("SequenceName\tclass\tcomplete\tstart\tend\tNB-ARC");
		out.newLine();
		
		Vector<String> vv= new Vector<String>();
		for(Enumeration<String> myenum = h.keys(); myenum.hasMoreElements();){
			vv.add(myenum.nextElement());
		}
		Collections.sort(vv);
		
		for(Enumeration<String> myenum = vv.elements(); myenum.hasMoreElements();){
			String key = myenum.nextElement();
			Vector<MastMotifHit> v = h.get(key);
			Collections.sort(v);
			String firstMotif = MastMotifHit.checkMotifCombination(v);
			if(!firstMotif.equalsIgnoreCase("FALSE")){
				out.write(key + "\t");
				
				out.write(MastMotifHit.getNBLRR_Class(v) + "\t");
				
				boolean complete = MastMotifHit.isComplete(v);
				if(complete){
					out.write("complete");
				}else{
					out.write("partial");
				}
				
				String startPosition = "undetermined";
				String endPosition = MastMotifHit.getEnd(v);
				if(starts.containsKey(v.get(0).getSequence_name())){
					startPosition = starts.get(v.get(0).getSequence_name())[0];
				}
				
				
				
				
				out.write("\t" + startPosition + "\t" + endPosition );
				out.write("\t" + MastMotifHit.getNBARC_Sequence(v));
				out.newLine();
			}
			
		}
		
		
		out.close();
		
	}
	
	*/
	
	/*
	public static Hashtable<String,String[]> getStarts(File blastFile, String splitPattern)throws IOException{
		Hashtable<String,String[]> h = new Hashtable<String,String[]>();
		BlastReader reader = new BlastXMLReader(blastFile);
		for( BlastIteration it = reader.readIteration(); it != null; it = reader.readIteration()){
			String start = "fw:"+it.getHits().get(0).getQueryStart();
			String end = "fw:"+it.getHits().get(0).getQueryEnd();
			if(it.getHits().get(0).getHSPs().get(0).getQueryStrand()==-1){
				start = "rv:"+it.getHits().get(0).getQueryEnd()+"";
				end = "rv:"+it.getHits().get(0).getQueryStart()+"";
			}
			String[] i ={start,end};
			h.put(it.getQueryID(), i);
			h.put(it.getQueryID().split(splitPattern)[0], i);
		}
		
		return h;
		
	}
	
	*/
	
	
	public void writeGFF(File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));

		out.write("##gff-version 2");
		out.newLine();
		out.write("##source-version MastParser " + NLRParser.getVersion());
		out.newLine();
		
		DateFormat dateFormat = new SimpleDateFormat("yyyy-MM-dd HH:mm");
		Date date = new Date();
		out.write("##date "+dateFormat.format(date)); 
		out.newLine();
		out.write("##Type DNA");
		out.newLine();
		out.write("#seqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute");
		out.newLine();
		
		
		for(Enumeration<String> myenum1 = this.getNBLRRs().keys(); myenum1.hasMoreElements();){
			String sequenceName = myenum1.nextElement();
			MastMotifHitList list = this.getNBLRRs().get(sequenceName);
		
			for( Enumeration<MastMotifHit> myenum2 = list.getMotifs().elements(); myenum2.hasMoreElements();){
				MastMotifHit hit = myenum2.nextElement();
				out.write(sequenceName +"\tMastParser\tMastMotif\t"  );
				
				int frame = Integer.parseInt(hit.getSequence_name().split("_frame")[1].substring(1,2) );
				
				boolean isForwardStrand = true;
				if(hit.getSequence_name().split("_frame")[1].substring(0, 1).equalsIgnoreCase("-")){
					isForwardStrand= false;
				}
				int aa_start = hit.getPos();
				int aa_end = hit.getPos() + hit.getMatch().length();
				int nucl_start = -1;
				int nucl_end = -1;
				try{
					nucl_start = MastFile.aminoAcidPositionToNucleotidePosition(aa_start, frame,isForwardStrand, list.getLength(), true);
					nucl_end   = MastFile.aminoAcidPositionToNucleotidePosition(aa_end, frame,isForwardStrand, list.getLength(), false);
					if( !isForwardStrand){
						nucl_start = MastFile.aminoAcidPositionToNucleotidePosition(aa_end, frame,isForwardStrand, list.getLength(), false);
						nucl_end   = MastFile.aminoAcidPositionToNucleotidePosition(aa_start, frame,isForwardStrand, list.getLength(), true);
					}
				}catch (OutOfFrameException e ){
					System.out.println(hit.getSequence_name() +" does not encode a proper frame");
				}
				String strand = "+";
				if(!isForwardStrand){
					strand = "-";
				}
				
				
				out.write(nucl_start + "\t" + nucl_end + "\t" +  hit.getPvalue() + "\t" +strand +"\t" + (Math.abs(frame) ) + "\t" + "name "+hit.getMotif() );
				out.newLine();
				
				
				
			}
			
		}	
		
		
		
		if(this.blast != null){
			for(Enumeration<String> myenum1 = this.blast.keys(); myenum1.hasMoreElements();){
				String sequenceName = myenum1.nextElement();
				for(Enumeration<int[]> myenum2= blast.get(sequenceName).elements(); myenum2.hasMoreElements();){
					int[] a = myenum2.nextElement();
					out.write(sequenceName+"\tBlast\texon\t" + a[0] + "\t" +a[1]+"\t.\t.\t.\t"  );
					out.newLine();
				}
			}
		}
		
		
		
		
		out.close();
		
		
		
		
	}
	
	/**
	 * 
	 * translates an aminoacid position back to the corresponding nucleotide position. It can be defined if the output position should be the first or the last base of the triplet
	 * 
	 * 
	 * 
	 * @param aaPos
	 * 			amino acid position
	 * @param frame
	 * 			the frame of the sequence (i.e. -3, -2, -1, 1, 2 or 3)
	 * @param sequenceLength
	 * 			the length of the amino acid sequence. In the reverse frames the position is dependent on the length.
	 * @param start
	 * 			true will give the position of the first base in the triplet, false the third.
	 * @return
	 * 			the corresponding nucleotide position.
	 * @throws OutOfFrameException
	 * 			thrown if the denoted reading frame is bullshit. 
	 */
	public static int aminoAcidPositionToNucleotidePosition(int aaPos, int frame,boolean isForwardStrand, int sequenceLength, boolean start)throws OutOfFrameException{
		
		
		
		if(  frame < -2 || frame >2){
			throw new OutOfFrameException("You gave this method a frame of " + frame + ". We are talking about reading frames on nucleotide sequences here.");
		}
	
		
		if( start){
		
			if(isForwardStrand){
				return ((aaPos-1) * 3 + frame + 1);
			
			}else{
			
				return (sequenceLength - ( (aaPos-1) *3) - (Math.abs(frame) )  );
			}
			
			
		}else{
			
			if(frame >0){
				return ((aaPos-1) * 3 + frame) +3;
			}
			
			if(frame <0){
				return (sequenceLength - ( (aaPos-1) *3) - (Math.abs(frame) ) -2  );
			}
			
			
		}
		return -1;
	}
	
	public Hashtable<String, MastMotifHitList> getNBLRRs()throws IOException{
		if(nblrrs == null){
			this.mergeEntries();
		}
		return this.nblrrs;
	}
	
	
	public Hashtable<String,MastEntry> getEntries(){
		return this.entries;
	}
	
}
