package coreClasses;

import java.util.Collections;
import java.util.Enumeration;
import java.util.Vector;

public class MastMotifHitList {
	Vector<MastMotifHit> motivlist;
	int length;
	boolean sorted;
	
	public MastMotifHitList(Vector<MastMotifHit> v){
		this.motivlist = v;
		Collections.sort(v);
		sorted = true;
	}
	public MastMotifHitList(){
		this.motivlist = new Vector<MastMotifHit>();
		sorted = true;
	}
	public void addMotif(MastMotifHit hit){
		motivlist.add(hit);
		sorted = false;
	}
	
	public void sort(){
		Collections.sort(motivlist);
		sorted = true;
	}
	
	public void setLength(int length){
		this.length = length;
	}
	public int getLength(){
		return this.length;
	}
	
	public  String checkMotifCombination(){
		if(!sorted){
			sort();
		}
		
		
		Vector<int[]> motivCombinations = new Vector<int[]>();
		int[] b = {17, 16}; motivCombinations.add(b);
		b = new int[3]; b[0] = 1; b[1] = 6; b[2] = 4; motivCombinations.add(b);
		b = new int[3]; b[0] = 6; b[1] = 4; b[2] = 5; motivCombinations.add(b);
		b = new int[3]; b[0] = 4; b[1] = 5; b[2] = 10;motivCombinations.add(b);
		b = new int[3]; b[0] = 5; b[1] = 10; b[2] = 3;motivCombinations.add(b);
		b = new int[3]; b[0] = 10; b[1] = 3; b[2] = 12; motivCombinations.add(b);
		b = new int[3]; b[0] = 3; b[1] = 12; b[2] = 2; motivCombinations.add(b);
		b = new int[3]; b[0] = 12; b[1] = 2; b[2] = 8; motivCombinations.add(b);
		b = new int[3]; b[0] = 2; b[1] = 8; b[2] = 7; motivCombinations.add(b);
		b = new int[3]; b[0] = 8; b[1] = 7; b[2] = 9; motivCombinations.add(b);
		b = new int[3]; b[0] = 7; b[1] = 9; b[2] = 11;motivCombinations.add(b);
		b = new int[2]; b[0] = 9; b[1] = 11;  motivCombinations.add(b);
		b = new int[2]; b[0] = 18; b[1] = 15; motivCombinations.add(b);
		b = new int[2]; b[0] = 15; b[1] = 13; motivCombinations.add(b);
		b = new int[2]; b[0] = 13; b[1] = 1; motivCombinations.add(b);
		b = new int[3]; b[0] = 1; b[1] = 4; b[2]=5 ; motivCombinations.add(b);
		
		
		
		
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		
		
		for( Enumeration<int[]> myenum = motivCombinations.elements(); myenum.hasMoreElements();){
			b = myenum.nextElement();
			for (int i = 0; i<= a.length-b.length; i++ ){
				if( a[i] == b[0]){
					boolean found = true;
					for(int j = 1; j<b.length; j++){
						if(b[j] != a[i+j]){
							found = false;
						}
					}
					if(found){
						
						String s = "";
						for( int j = 0; j< b.length; j++){
							s = s + ","+b[j];
						}
						return s.substring(1);
					}
				}
				
			}
			
			
		}	
		return "FALSE";
	}
	
	
	/**
	 * This checks if the NBLRR is complete based on the existence of certain motifs. Note that this method not specifically checks if it is a NBLRR. 
	
	 * @return
	 */
	public  boolean isComplete(){
		if(!sorted){
			sort();
		}
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		
		boolean hasNTerminal = false;
		boolean hasNBARC = false;
		//boolean hasLRR = false;
		
		for(int i = 0; i< a.length; i++){
			int m = a[i];
			if(m ==16 || m == 18){
				hasNTerminal = true;
			}
			if(m == 3 && hasNTerminal){
				hasNBARC=true;
			}
			if(m == 11 && hasNBARC){
				//hasLRR = true;
				return true;
			}
			
		}
		return false;
		
	}
	
	public  int[] getMotifList(){
		
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		return a;
	}	
	
	public  String getNBLRR_Class(){
		if(!sorted){
			sort();
		}
		
		int[] a = new int[motivlist.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(motivlist.get(i).getMotif().split("_")[1]);
		}
		
		
		
		
		for( int i = 0; i< a.length; i++){
			
			
			try{
				if(a[i] == 2 || (a[i]==1 && a[i+1] ==6  ) ||  (a[i] == 17 && a[i+1] ==16) ){
					return "CNL";
				}
			}catch (ArrayIndexOutOfBoundsException e){}	
			
		}
		
		
		for( int i = 0; i< a.length; i++){
			try{
				if(a[i] ==18 && a[i+1] == 15 && a[i+2] ==13){
					return "TNL";
				}
			}catch (ArrayIndexOutOfBoundsException e){}
			
				
			try{
				if( a[i] ==1 && a[i+1]!=6   && (a[i+1] ==4   || a[i+2] ==4)    ){
					return "TNL";
				}
			}catch (ArrayIndexOutOfBoundsException e){}	
			
		}
		
		
		
		
		return "N/A";
		
	}
	
	/**
	 * gives the name of the sequence and the position the first motif is found. Sequence-name and position are separated with " : "
	 * 
	 * 
	 * 		A vector of MastMotifHits
	 * @return
	 */
	public  String getStart(){
		if(!sorted){
			sort();
		}
		
		MastMotifHit hit = this.motivlist.get(0);
		return hit.getSequence_name() +":" + hit.getPos();
	}
	
	/**
	 * gives the name of the sequence and a putative end position. This is the first stop codon or undefined base after the last motif. 
	 * 
	 * 
	 * 		A vector of MastMotifHits
	 * @return
	 */
	public  String getEnd(){
		if(!sorted){
			sort();
		}
		MastMotifHit hit = motivlist.get(motivlist.size()-1);
		return hit.getSequence_name() + ":" + (hit.getPos() + hit.getMatch().length());
		/*
		char[] sequence = hit.getAaSequence().toCharArray();
		
		for( int i = hit.getPos()+hit.getMatch().length(); i < sequence.length; i++){
			if(sequence[i] =='*' || sequence[i] == 'X'){
				return hit.getSequence_name() + ":" + (i);
			}
		}
		return hit.getSequence_name() + ":" + sequence.length;
		*/
		
		
		
	}
	
	
	/**
	 * gives the sequence starting from 10 aa before the Motif 1 and 200 aa after beginning of Motif 1. 
	 * 
	 * 
	 * @return
	 */
	public  String getNBARC_Sequence(){
		if(!sorted){
			sort();
		}
		for( Enumeration<MastMotifHit> myenum = motivlist.elements(); myenum.hasMoreElements();){
			MastMotifHit hit = myenum.nextElement();
			if(hit.getAaSequence() == null){
				return "N/A";
			}
			if(Integer.parseInt(hit.getMotif().split("_")[1]) == 1 ){
				int start = Math.max(0, hit.getPos()-11);
				int end = Math.min(hit.getPos() + 200, hit.getAaSequence().length()-1);
				//System.out.println(hit.sequence_name + "\t" + start + "\t" + end + "\t" + hit.getAaSequence());
				return hit.getAaSequence().substring(start, end);
			}
		}
		return "N/A";
	}

	public  boolean hasPotentiallyTwoNBLRRs(){
		if(!sorted){
			sort();
		}
		
		
		
		boolean foundLRR = false;
		
		int[] a = getMotifList();
		
		for( int i = 0; i< a.length; i++){
			
			if(a[i] == 9   ||  a[i] == 11 || a[i] ==19){
				foundLRR=true;
			}
			if(a[i] ==1 && foundLRR){
				return true;
			}
		}
		
		
		return false;
	}
	
	public String getMotifListString(){
		
		if(!sorted){
			sort();
		}
		
		int[] a = getMotifList();
		String s= "";
		for( int i = 0; i< a.length; i++){
			s = s + "," + a[i];
		}
		return s.substring(1);
		
	}
	
	public int getSize(){
		return this.motivlist.size();
	}
	public Vector<MastMotifHit> getMotifs(){
		return this.motivlist;
	}
}
