package coreClasses;

import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

public class MastEntry {

	
	String id;
	String name;
	int length;
	double combined_p_value;
	double p_value;
	double e_value;
	
	
	Vector<MastMotifHit> hits;
	String aaSequence;
	
	public MastEntry(String entryString){
		
		Pattern p = Pattern.compile("<sequence\\s+id=\"(\\w+)\"\\s+db=\"(\\w+)\"\\s+num=\"(\\d+)\"\\s+name=\"(.+?)\"\\s+comment=\"(.*?)\"\\s+length=\"(\\d+)\">");
		Matcher m = p.matcher(entryString);
		m.find();
		
		this.id = m.group(1);
		this.name = m.group(4);
		this.length = Integer.parseInt(m.group(6));
		this.combined_p_value = Double.parseDouble(entryString.split("combined_pvalue=\"")[1].split("\"")[0]);
		this.e_value = Double.parseDouble(entryString.split("evalue=\"")[1].split("\"")[0]);
		
		
		
		
		
		StringBuilder builder = new StringBuilder("");
		int seqstart = 1;
		String sequence = "";
		p = Pattern.compile("<seg start=\"(\\d+)\">\\s*<data>([\\w*]+)\\s*</data>");
		m = p.matcher(entryString);
		while(m.find ()){
			seqstart = Integer.parseInt(m.group(1));
			sequence = m.group(2);
	
			int addX = seqstart - 1 - builder.length();
			for( int i = 0; i< addX; i++){
				builder.append("X");
			}
			builder.append(sequence);
		}
		this.aaSequence = builder.toString();
				
		
		
		
		
		hits = new Vector<MastMotifHit>();
		
		   p = Pattern.compile("<hit\\s+pos=\"(\\d+)\"\\s+gap=\"(\\d+)\"\\s+motif=\"(\\w+)\"\\s+pvalue=\"([\\w-\\.]+)\"\\s+match=\"([\\s+]+)\"/>");
		 m = p.matcher(entryString);
		
		while(m.find()){
			MastMotifHit hit = new MastMotifHit(this.name, m.group(0) );
			//hit.setAaSequence(this.aaSequence);
			//if( hit.getSequence_name().equalsIgnoreCase("seq_638_frame+3")){
				
			//}
			hits.add( hit );
		}
		
		
		
	}


	public String getId() {
		return id;
	}


	public String getName() {
		return name;
	}


	public int getLength() {
		return length;
	}


	public double getP_value() {
		return p_value;
	}


	public double getE_value() {
		return e_value;
	}


	public Vector<MastMotifHit> getHits() {
		return hits;
	}
	
	
	public String getAASequence(){
		return this.aaSequence;
	}
	public void setAASequence(String aaSequence){
		this.aaSequence = aaSequence;
	}
	
}






