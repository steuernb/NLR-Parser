package coreClasses;

import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MastMotifHit implements Comparable<MastMotifHit>{

	String sequence_name;
	int pos;
	int gap;
	String motif;
	double pvalue;
	String match;
	String aaSequence;
	
	public MastMotifHit(String sequence_name, String hitString){
		Pattern p = Pattern.compile("<hit\\s+pos=\"(\\d+)\"\\s+gap=\"(\\d+)\"\\s+motif=\"(\\w+)\"\\s+pvalue=\"([\\w-\\.]+)\"\\s+match=\"([\\s+]+)\"/>");
		Matcher m = p.matcher(hitString);
		m.find();
		this.pos=Integer.parseInt(m.group(1));
		this.gap=Integer.parseInt(m.group(2));
		this.motif = m.group(3);
		this.pvalue = Double.parseDouble(m.group(4));
		this.match = m.group(5);
		this.sequence_name = sequence_name;
			
	}

	public String getSequence_name() {
		return sequence_name;
	}

	public int getPos() {
		return pos;
	}

	public int getGap() {
		return gap;
	}

	public String getMotif() {
		return motif;
	}

	public double getPvalue() {
		return pvalue;
	}

	public String getMatch() {
		return match;
	}
	
	public String getAaSequence(){
		return this.aaSequence;
	}
	
	
	public void setAaSequence(String aaAequence){
		this.aaSequence = aaAequence;
	}
	
	public int compareTo(MastMotifHit m2){
		if( this.getPos() < m2.getPos()){
			return -1;
		}
		else if( this.getPos() > m2.getPos()){
			return 1;
		}else{
			return 0;
		}
		
		
	}
	/*
	public static String checkMotifCombination(Vector<MastMotifHit> v){
		
		Collections.sort(v);
		
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
		
		
		
		
		
		int[] a = new int[v.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(v.get(i).getMotif().split("_")[1]);
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
	
	
	
	public static boolean isComplete(Vector<MastMotifHit> v){
		Collections.sort(v);
		
		int[] a = new int[v.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(v.get(i).getMotif().split("_")[1]);
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
	
	public static int[] getMotifList(Vector<MastMotifHit> v){
		Collections.sort(v);
		
		int[] a = new int[v.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(v.get(i).getMotif().split("_")[1]);
		}
		return a;
	}	
	
	public static String getNBLRR_Class(Vector<MastMotifHit> v){
		Collections.sort(v);
		
		int[] a = new int[v.size()];
		for( int i = 0; i< a.length; i++){
			a[i] = Integer.parseInt(v.get(i).getMotif().split("_")[1]);
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
	
	
	public static String getStart(Vector<MastMotifHit> v){
		Collections.sort(v);
		MastMotifHit hit = v.get(0);
		return hit.getSequence_name() +":" + hit.getPos();
	}
	
	
	public static String getEnd(Vector<MastMotifHit> v){
		Collections.sort(v);
		MastMotifHit hit = v.get(v.size()-1);
		char[] sequence = hit.getAaSequence().toCharArray();
		
		for( int i = hit.getPos()+hit.getMatch().length(); i < sequence.length; i++){
			if(sequence[i] =='*' || sequence[i] == 'X'){
				return hit.getSequence_name() + ":" + (i);
			}
		}
		return hit.getSequence_name() + ":" + sequence.length;
		
		
		
		
	}
	
	
	
	public static String getNBARC_Sequence(Vector<MastMotifHit> v){
		Collections.sort(v);
		for( Enumeration<MastMotifHit> myenum = v.elements(); myenum.hasMoreElements();){
			MastMotifHit hit = myenum.nextElement();
			if(Integer.parseInt(hit.getMotif().split("_")[1]) == 1 ){
				int start = Math.max(0, hit.getPos()-11);
				int end = Math.min(hit.getPos() + 200, hit.getAaSequence().length()-1);
				return hit.getAaSequence().substring(start, end);
			}
		}
		return "N/A";
	}

	public static boolean hasPotentiallyTwoNBLRRs(Vector<MastMotifHit> v){
		boolean foundLRR = false;
		
		int[] a = getMotifList(v);
		
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
	
	public static String giveMotifList(Vector<MastMotifHit> v){
		
		int[] a = getMotifList(v);
		String s= "";
		for( int i = 0; i< a.length; i++){
			s = s + "," + a[i];
		}
		return s.substring(1);
		
	}
	*/
}
