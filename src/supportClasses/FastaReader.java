package supportClasses;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;
/**
 * fastaReader is for extracting sequences from a
 * @author steuerna
 *
 */
public class FastaReader {
	
	File file;

	
	BufferedReader in ;
	
//	String header;
//	StringBuilder sequence;
	
	String lastline;
	
	boolean isQuality;
	
	
	/**
	 * 
	 * @param path
	 * 			path to the fastaFile
	 * @throws FileNotFoundException
	 */
	public FastaReader(String path)throws IOException{
		super();
		
		initialize(new File(path));
		
	}
	public FastaReader( File file)throws IOException{
		super();
		initialize(file);
		
	}
	
	

	
	public FastaReader(File file, boolean isQuality)throws IOException{
		super();
		this.isQuality = isQuality;
		initialize(file);
		
	}
	
	
	
	public void initialize(File file)throws IOException{
		this.file = file;
		//check if the fastq is gzipped
				FileInputStream fis = new FileInputStream(file);
				byte[] bytes = new byte[2];
				fis.read(bytes);
				int head = ((int) bytes[0] & 0xff) | ((bytes[1] << 8) & 0xff00);
				boolean gzip = GZIPInputStream.GZIP_MAGIC == head;
				fis.close();
				
				
					
				if(gzip){
					this.in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
				}else{
					this.in = new BufferedReader((new FileReader(file)));
				}
				
				this.lastline = in.readLine();
				while(lastline != null && !lastline.trim().startsWith(">")){
					this.lastline=in.readLine();
				}
				
				
	}
	
	
	
	
	public void close()throws IOException{
		in.close();
	}
	
	
	/**
	 * 
	 * 
	 * 
	 * @return
	 * @throws IOException
	 */
	public FastaSequence readEntry()throws IOException{
		if(lastline == null){
			return null;
		}
		if(!lastline.trim().startsWith(">")){
			throw new IOException("Problem reading the fasta file. Actually, this should not have happened. Review source code in FastaReader.java");
		}
		
		String id =  lastline.trim().split("\\s")[0].substring(1);
		String desc = "";
		try {
			desc = lastline.trim().substring(2+id.length());
		} catch (StringIndexOutOfBoundsException e) {
			
		}
		
		StringBuilder builder = new StringBuilder();
		lastline = in.readLine();
		while(lastline != null && !lastline.trim().startsWith(">")){
			builder.append(lastline.trim());
			lastline = in.readLine();
		}
		
		return new FastaSequence(id, desc, builder.toString());
	}
	
	
	
	/**
	 * moves quickliy to the first sequence with the given identifier. To get the FastaSequence of this identifier call readEntry after this method.
	 * @param identifier
	 * 			the identifier of the sequence
	 * @return
	 * 			true if the sequence was found,
	 * 			false if the fasta file does not contain this sequence.
	 * @throws IOException
	 */
	public boolean gotoEntry(String identifier)throws IOException{
		
		String inputline = lastline;
		if( inputline == null){
			inputline = "";
		}
		
		while(inputline != null && !inputline.startsWith(">"+identifier+" ") && !inputline.equalsIgnoreCase(">"+identifier)){
			
			inputline = in.readLine();
		}
		if(inputline == null){
			return false;
		}else{
			lastline = inputline;
			return true;
		}
		
		
	}
	
	/**
	 * @deprecated I use the StringBuilder now. That should be quick enough.
	 * Read the next entry from the fasta file but only return a substring of that sequence. 
	 * That's much faster than read all first an then clip. 
	 * Recommended for large sequences.
	 * 
	 * @param start
	 * 			left clip
	 * @param stop
	 * 			right clip
	 * @return
	 * @throws IOException
	 */
	public FastaSequence readSubSequenceFromEntry(int start, int stop) throws IOException {
		String header = new String();
		StringBuilder sequence = new StringBuilder();
		boolean endOfFile = false; 
		
		int baseCounta = 0;
		String inputline;
		
		if(        ( lastline != null )    &&     ( lastline.startsWith(">") )        ){ 
			inputline = lastline;
		}else{
			
			inputline = in.readLine();
	       
			while ( inputline != null && !inputline.startsWith( ">" )){
				inputline = in.readLine();
			}
		}
		
		if( inputline != null){
			header = inputline.substring(1);
		}else{
			endOfFile = true;
		}
		
		inputline = in.readLine();
		sequence = new StringBuilder();
		
		while ( ( inputline != null ) && ( !inputline.startsWith(">") ) ){
				if(!inputline.trim().equalsIgnoreCase("")){
					int i = inputline.length();
					if (baseCounta+i<start || baseCounta> stop){
						inputline = "";
					}
					if (baseCounta < start && baseCounta + i >=start){
						inputline = inputline.substring(  start - baseCounta);
					}
					if( baseCounta <= stop && baseCounta + i > stop){
						//System.out.println(baseCounta + " " + i + " " + stop);
						int j = baseCounta + i - stop;
						inputline = inputline.substring(0, inputline.length()-j);
					}
					baseCounta = baseCounta + i;
					sequence = sequence.append(inputline);
					//System.out.println(sequence);
					if(baseCounta > stop ){
						break;
					}
				}	
			
			inputline = in.readLine();
		}
		
		try{
			lastline = inputline;
		} catch(NullPointerException e){}	
		String id ="";
		try {
			id = header.split(" ")[0];
		} catch (Exception e) {
			return null;
		}
		

		String description = new String();
		
		//check whether the header contains a description. (to avoid StringIndexOutOfBoundsException)
		if(header.length()>id.length()+1){
			description = header.substring(id.length()+1);
		}
		
		
		
		FastaSequence f = new FastaSequence(id, description, sequence.toString());
		if( endOfFile ){ 
			return null;
		}else{
			return f;
		}	
		
		
		
	}
	
	
	
	/**
	 * @deprecated
	 * readEntry works as readLine but returns a fastaEntry
	 * 
	 */
	 public FastaSequence readEntry_old()throws IOException{
		String header = new String();
		StringBuilder sequence = new StringBuilder();
		boolean endOfFile = false; 
	 
		String inputline;
		
		if(        ( lastline != null )    &&     ( lastline.startsWith(">") )        ){ 
			inputline = lastline;
		}else{
			
			inputline = in.readLine();
	       
			while ( inputline != null && !inputline.startsWith( ">" )){
				inputline = in.readLine();
			}
		}
		
		if( inputline != null){
			header = inputline.substring(1);
		}else{
			endOfFile = true;
		}
		
		inputline = in.readLine();
		sequence = new StringBuilder();
		
		while ( ( inputline != null ) && ( !inputline.startsWith(">") ) ){
				if(!inputline.trim().equalsIgnoreCase("")){
					sequence = sequence.append(inputline);
				}	
			
			inputline = in.readLine();
		}
		
		try{
			lastline = inputline;
		} catch(NullPointerException e){}	
		String id ="";
		try {
			id = header.split(" ")[0];
		} catch (Exception e) {
			return null;
		}
		

		String description = new String();
		
		//check whether the header contains a description. (to avoid StringIndexOutOfBoundsException)
		if(header.length()>id.length()+1){
			description = header.substring(id.length()+1);
		}
		
		
		
		FastaSequence f = new FastaSequence(id, description, sequence.toString());
		if( endOfFile ){ 
			return null;
		}else{
			return f;
		}	
	}
	
	 

	 
	 
	 
	 public File getFile(){
		 return this.file;
	 }
	
	public String getFileName(){
		return file.getName();
	}
	
	
	public static void main(String[] args){
		try {
			File dir = new File("/Users/steuernb/Documents/projects/test");
			FastaReader reader = new FastaReader(new File(dir, "test.fasta"));
			FastaSequence seq = reader.readEntry();
			
			while (seq != null) {
				seq.setDescription(seq.getLength()+"");
				System.out.println(seq.getFormatedFastaString(20));
				
					seq = reader.readEntry();
			}
			
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
}
