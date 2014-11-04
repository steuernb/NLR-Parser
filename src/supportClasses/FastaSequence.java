package supportClasses;

import java.util.Hashtable;
import java.util.Vector;


/**
 * A sequence in FASTA format has following format:
 *
 * >identifier description the raw sequence in one or more lines.
 *
 * The identifier is required the description is optional.
 *
 * @author steuerna,
 * @author spies
 *
 */
public class FastaSequence implements Comparable<FastaSequence>{
    String identifier;
    String description;
    String sequence;
    
    
    Hashtable<String,String> geneticCode;
    
    
    public FastaSequence(String identifier, String sequence){
    	super();
    	this.identifier = identifier;
    	this.description = "";
    	this.sequence = sequence;
    }
    
    /**
     *
     * @param identifier
     *            identifier the identifier of the sequence (without the leading
     *            ">")
     * @param description
     *            the description of the sequence. It must not contain any
     *            newline character.
     * @param sequence
     *            the raw sequence
     */
    public FastaSequence(String identifier, String description, String sequence) {
        super();
        this.identifier = identifier;
        this.description = description;
        this.sequence = sequence.replaceAll("\n", "");
        
    }

    /**
     *
     * @param fasta
     *            ONE sequence in FASTA fomat. >identifier description <newline>
     *            sequence
     * @throws Exception
     *             If the ingoing String does not begin with ">" or if the
     *             String contains more than one sequence.
     */
    public FastaSequence(String fasta) throws FastaFormatException {
        String myfasta = fasta.trim();
        while (myfasta.startsWith("\n")) {
            myfasta = myfasta.substring(1);
            myfasta = myfasta.trim();
        }
        if (!myfasta.startsWith(">")) {
            throw (new FastaFormatException("this is no fastaformat. Fastaformat should start with >."));
        }

        String[] lines = myfasta.split("\n");

        this.identifier = lines[0].substring(1).split(" ")[0];
        this.description = "";
        try{
        	this.description = lines[0].substring(this.identifier.length()+2);
        }catch(StringIndexOutOfBoundsException e){}

        this.sequence = new String();
        for (int i = 1; (i < lines.length); i++) {
            if (lines[i].startsWith(">")) {
                throw new FastaFormatException("there is more than one sequence in the string " + fasta);
            } else {
                this.sequence = this.sequence + lines[i].trim();
            }

        }
        
    }

    
    
    
    
    
    
    /**
     * Returns a new FastaSequence clipped according to a quality threshold.
     * 
     * the new sequence is the largest region of the old sequence where ALL quality values are above the threshold.
     * 
     * @param phredThreshold
     * 		The threshold for clipping. The new sequence will not contain any base with a quality below it.
     * 			
     * @return
     * 		a FastaSequence representing the clipped sequence.
     */
    
    
    
    public String getDescription() {
        return description;
    }

    public void setDescription(String description) {
        this.description = description;
    }

    public String getIdentifier() {
        return identifier;
    }

    public void setIdentifier(String identifier) {
        this.identifier = identifier;
    }

    public String getSequence() {
        return sequence;
    }

    public String getUnpaddedSequence(){
    	return sequence.replace("-", "");
    }
    
    
    public int getLength(){
    	return this.sequence.length();
    }
    
    public int getLengthWithoutNs(){
    	String s = new String(this.getSequence());
    	
    	return s.replaceAll("[Nn]", "").length();
    }
    
    public int length(){
    	return getLength();
    }
    
    public void setSequence(String sequence) {
        this.sequence = sequence;
    }

    /**
     *
     * @return the sequence in fasta format.
     */
    public String getFastaString() {

        return (">" + identifier + " " + description).trim() + "\n" + sequence + "\n";
    }

    /**
     *
     * @param lettersPerLine
     *            the number of sequence letters per line.
     * @return the String with the formatted fasta sequence: >identifier
     *         description the fromatted sequence.
     */
    public String getFormatedFastaString(int lettersPerLine) {
        String header = ">" + identifier + " " + description;
        String formattedSequence = new String();
        String restSequence = new String(sequence);
        int length = sequence.length();

        while (lettersPerLine <= length) {
            formattedSequence = formattedSequence + "\n"
                    + restSequence.substring(0, lettersPerLine);
            restSequence = restSequence.substring(lettersPerLine, length);
            length = length - lettersPerLine;
        }
        formattedSequence = formattedSequence + "\n" + restSequence;

        if (!formattedSequence.endsWith("\n")){
        	formattedSequence = formattedSequence+"\n";
        }
        return header + formattedSequence ;
    }

    /**
	 * returns the clipped sequence. It does not remove the clipped based but replaces them by an X.
	 * @param preX left clipping border
	 * @param sufX right clipping border.
	 * @return only the sequence itself without a header
	 */
	public String getMaskedSequence(int preX, int sufX){
		String s = new String(sequence);
		String pre=new String("");
		String suf=new String("");
		for (int i = 0;i <preX; i++){
			pre = pre+"X";
		}
		for (int i = 0;i <sufX; i++){
			suf = pre+"X";
		}		
		//System.out.println(identifier+":   "+preX+"  "+sufX);
		return pre+s.substring(preX, sequence.length()-sequence.length()+sufX)+suf;
	}
	/**
	 * Tells weather the sequence contains IUPAC Bases or not.
	 * 
	 * @return
	 */
	public boolean containsIUPACSymbol(){
		boolean b = false;
		if(
			   sequence.toUpperCase().contains("R")
			|| sequence.toUpperCase().contains("Y")
			|| sequence.toUpperCase().contains("M")
			|| sequence.toUpperCase().contains("K")
			|| sequence.toUpperCase().contains("S")
			|| sequence.toUpperCase().contains("W")
			|| sequence.toUpperCase().contains("B")
			|| sequence.toUpperCase().contains("D")
			|| sequence.toUpperCase().contains("H")
			|| sequence.toUpperCase().contains("V")
			|| sequence.toUpperCase().contains("N")
		
		)
		
		{
			b = true;
		}
		
		
		return b;
		
	}
	/**
	 * this gives the positions of IUPAC Bases in a sequence.
	 * @return
	 * 		an int[] with the positions of IUPAC Bases in the sequence.
	 */
	public int[] giveIUPACPositions(){
		int[] array = new int[0];
		char[] charsequence = sequence.toUpperCase().toCharArray();
		for(int i = 0; i < charsequence.length; i++){
			char c = charsequence[i];
			if(c == 'R'||c == 'Y'||c == 'M'||c == 'K'||c == 'S'||c == 'W'||c == 'B'||c == 'D'||c == 'H'||c == 'V'||c == 'N'){
				int[] oldarray = array;
				array= new int[oldarray.length+1];
				for(int j = 0; j< oldarray.length; j++){
					array[j]=oldarray[j];
				}
				array[oldarray.length] = i;
			}
			
			
		}
		
		
		
		return array;
	}
	
	
	
	/**
	 * returns the clipped sequence.
	 * @param preCut left clipping border
	 * @param sufCut right clipping border
	 * @return
	 */
	public String getTrimmedSequence(int preCut, int sufCut){
		String s = new String(sequence);
		try {
			return s.substring(preCut, sufCut);
		} catch (StringIndexOutOfBoundsException e) {
			
			
			try {
				s = s.substring(0,sufCut);
			} catch (StringIndexOutOfBoundsException e2) {
			
			}
			try {
				s = s.substring(preCut);
			} catch (StringIndexOutOfBoundsException e2) {
				
			}
			
			return s;
			
		}
		
	}
	
	
	public String getSequenceReplaceBases(String oldBases, String newBases ){
		String newSequence = sequence.replace(oldBases, newBases);
		
		return newSequence;	
		
	}
	/**
	 * This one returns the longest stretch of Upper Case Letters, thus removes lower case ends
	 * It is used for sequences where quality clips are marked as lowerCase.
	 * @return 
	 * 			the sequence where all bases in lowerCase were clipped.
	 */
	public String getSequenceRemoveLowerCase(){
		String[] split = this.getSequence().split("[a-z]+");
		String outSequence = "";
		for( int i=0; i<split.length;i++){
			if(split[i].length() > outSequence.length()){
				outSequence = split[i];
			}
		}
		
		return outSequence;
		/*
		StringBuffer sb = new StringBuffer(sequence);
		
		StringBuffer newSb = new StringBuffer();
		
		for (int i = 0; i< sequence.length(); i++){
			if(Character.isUpperCase(sb.charAt(i))){
				newSb.append(sb.charAt(i));
			}
		}
		
		return newSb.toString();
		*/
	}
	public FastaSequence getFastaSequenceRemoveLowerCase(){
		FastaSequence seq = new FastaSequence(this.getIdentifier(), this.getDescription(), this.getSequence());
		
		StringBuffer sb = new StringBuffer(seq.getSequence());
		
		StringBuffer newSb = new StringBuffer();
		
		for (int i = 0; i< sequence.length(); i++){
			if(Character.isUpperCase(sb.charAt(i))){
				newSb.append(sb.charAt(i));
			}
		}
		seq.setSequence(newSb.toString());
		
		
		return seq;
	}
	
	public FastaSequence getFastaSequenceRemoveNs(){
		
		
		String s = this.getSequence();
		while(s.substring(0,1).equalsIgnoreCase("N")){
			s = s.substring(1);
		}
		while(s.substring(s.length()-1).equalsIgnoreCase("N")){
			s = s.substring(0, s.length()-1);
		}
		FastaSequence seq = new FastaSequence (this.getIdentifier(), this.getDescription(), s);
		
		
		return seq;
	}
	
	/**
	 * set all "N"s in the sequence to lower case. In Blast lower case can be masked. N on the other hand is considered a perfect match to everything else.
	 * 
	 */
	public void setNsToLowerCase(){
		char[] c = this.getSequence().toCharArray();
		for( int i = 0; i< c.length; i++){
			if(c[i] =='N'){
				c[i] = 'n';
			}
		}
		this.setSequence(new String(c));
	}
	

    public String getReverseSequence() {
        return new StringBuilder(getSequence()).reverse().toString();
    }

    public String getReverseComplementarySequence(){
    	StringBuilder seq = new StringBuilder(getTranscribeSequence());
    	
    	return seq.reverse().toString();
    }
    
    /**
     * 
     * @return
     */
    public String getTranscribeSequence() {
        StringBuilder seq = new StringBuilder(getSequence());
        StringBuilder transcribe = new StringBuilder(seq.length());

        for (int i = 0; i < seq.length(); i++) {
            char t = seq.charAt(i);
            switch (t) {
            case 'A':
                transcribe.append('T');
                break;
            case 'T':
                transcribe.append('A');
                break;
            case 'G':
                transcribe.append('C');
                break;
            case 'C':
                transcribe.append('G');
                break;
            case 'a':
                transcribe.append('t');
                break;
            case 't':
                transcribe.append('a');
                break;
            case 'g':
                transcribe.append('c');
                break;
            case 'c':
                transcribe.append('g');
                break;    
            default:
            	transcribe.append(t);
                break;
            }
        }
        return transcribe.toString();
    }
    
    /**
     * 
     */
	   public int compareTo(FastaSequence arg0) {
			FastaSequence o2 =  arg0;
	    	int compared= 0;
	    	if(this.getLength()< o2.getLength()){
	    		compared = 1;
	    	}else
	    	if(this.getLength()> o2.getLength()){
	    		compared = -1;
	    	}else
	    	if(this.getLength()== o2.getLength()){
	    		compared = 0;
	    	}
	    	
	      return compared;		
	   }
	   
	   
	   
	   
	   
	   
    public char[] toArray(){
    	return this.sequence.toCharArray();
    }
    
    /**
     * calculates the number of Bases written in upper case for the sequence. In case clipped bases are marked by lower case this will give the actual length of the sequence.
     * @return
     * 			number of Bases in Upper Case
     */
    public int getNumUpperCaseBases(){
    	char[] c = this.toArray();
    	int numBases = 0;
    	
    	for(int i=0; i< c.length; i++){
    		if(Character.isUpperCase(c[i])){
    			numBases++;
    		}
    		
    	}
    	return numBases;
    }
   
    
    
    
    private void loadGeneticCodeTable(){
    	this.geneticCode = new Hashtable<String,String>();
    	this.geneticCode.put("ATG", "M");
    	this.geneticCode.put("TTG", "W");
    	this.geneticCode.put("TAT", "Y");
    	this.geneticCode.put("TAC", "Y");
    	this.geneticCode.put("TTT", "F");
    	this.geneticCode.put("TTC", "F");
    	this.geneticCode.put("TGT", "C");
    	this.geneticCode.put("TGC", "C");
    	this.geneticCode.put("AAT", "N");
    	this.geneticCode.put("AAC", "N");
    	this.geneticCode.put("GAT", "D");
    	this.geneticCode.put("GAC", "D");
    	this.geneticCode.put("CAA", "Q");
    	this.geneticCode.put("CAG", "Q");
    	this.geneticCode.put("GAA", "E");
    	this.geneticCode.put("GAG", "E");
    	this.geneticCode.put("CAT", "H");
    	this.geneticCode.put("CAC", "H");
    	this.geneticCode.put("AAA", "K");
    	this.geneticCode.put("AAG", "K");
    	this.geneticCode.put("ATT", "I");
    	this.geneticCode.put("ATC", "I");
    	this.geneticCode.put("ATA", "I");
    	this.geneticCode.put("GGT", "G");
    	this.geneticCode.put("GGC", "G");
    	this.geneticCode.put("GGA", "G");
    	this.geneticCode.put("GAT", "A");
    	this.geneticCode.put("GCC", "A");
    	this.geneticCode.put("GCA", "A");
    	this.geneticCode.put("GTT", "V");
    	this.geneticCode.put("GTC", "V");
    	this.geneticCode.put("GTA", "V");
    	this.geneticCode.put("ACT", "T");
    	this.geneticCode.put("ACC", "T");
    	this.geneticCode.put("ACA", "T");
    	this.geneticCode.put("CCT", "P");
    	this.geneticCode.put("CCC", "P");
    	this.geneticCode.put("CCA", "P");
    	this.geneticCode.put("CTT", "L");
    	this.geneticCode.put("CTC", "L");
    	this.geneticCode.put("CTA", "L");
    	this.geneticCode.put("TCT", "S");
    	this.geneticCode.put("TCC", "S");
    	this.geneticCode.put("TCA", "S");
    	this.geneticCode.put("CGT", "R");
    	this.geneticCode.put("CGC", "R");
    	this.geneticCode.put("CGA", "R");
    	this.geneticCode.put("TAA", "*");
    	this.geneticCode.put("TAG", "*");
    	this.geneticCode.put("TGA", "*");
    	
    	
    	
    	
    }
    
    
    public char translateTriplet(String s)throws StringIndexOutOfBoundsException{
    	char aa = 'X';
    	if(s.length()!= 3){
    		throw new StringIndexOutOfBoundsException ("Length of the inputString has to be 3. Otherwise it is not a triplet and cannot be translated.");
    	}
    	
    	if(this.geneticCode == null){
    		this.loadGeneticCodeTable();
    	}
    	
    	if( this.geneticCode.containsKey(s)){
    		return this.geneticCode.get(s).charAt(0);
    	}
    	
    	
    	return aa;
    }
    
    
    
    
    public static char translateTriplet(char c1, char c2, char c3){
    	
    	
    	char aa = 'X';
    	
/*
M	Met 	AUG        
W	Trp 	UGG        
Y	Tyr 	UAU UAC    
F	Phe 	UUU UUC    
C	Cys 	UGU UGC    
N	Asn 	AAU AAC    
D	Asp 	GAU GAC    
Q	Gln 	CAA CAG    
E	Glu 	GAA GAG    
H	His 	CAU CAC    
K	Lys 	AAA AAG    
I	Ile 	AUU AUC AUA
G	Gly 	GGU GGC GGA
A	Ala 	GCU GCC GCA
V	Val 	GUU GUC GUA
T	Thr 	ACU ACC ACA
P	Pro 	CCU CCC CCA
L	Leu 	CUU CUC CUA
S	Ser 	UCU UCC UCA
R	Arg 	CGU CGC CGA
X	STOP 	UAA UAG UGA
*/
    		String triplet = new String();
    		triplet = (Character.toString(c1)+Character.toString(c2)+Character.toString(c3)).toUpperCase();
    		if( triplet.equalsIgnoreCase(	"ATG")		){	aa='M';	}
    		if( triplet.equalsIgnoreCase(	"TGG")      ){  aa='W';	}
    		if( triplet.equalsIgnoreCase(	"TAT")||triplet.equalsIgnoreCase("TAC") ){  aa='Y';	}
    		if( triplet.equalsIgnoreCase(	"TTT")||triplet.equalsIgnoreCase("TTC") ){  aa='F';	}
    		if( triplet.equalsIgnoreCase(	"TGT")||triplet.equalsIgnoreCase("TGC") ){  aa='C';	}
    		if( triplet.equalsIgnoreCase(	"AAT")||triplet.equalsIgnoreCase("AAC") ){  aa='N';	}
    		if( triplet.equalsIgnoreCase(	"GAT")||triplet.equalsIgnoreCase("GAC") ){  aa='D';	}
    		if( triplet.equalsIgnoreCase(	"CAA")||triplet.equalsIgnoreCase("CAG") ){  aa='Q';	}
    		if( triplet.equalsIgnoreCase(	"GAA")||triplet.equalsIgnoreCase("GAG") ){  aa='E';	}
    		if( triplet.equalsIgnoreCase(	"CAT")||triplet.equalsIgnoreCase("CAC") ){  aa='H';	}
    		if( triplet.equalsIgnoreCase(	"AAA")||triplet.equalsIgnoreCase("AAG") ){  aa='K';	}
    		if( triplet.equalsIgnoreCase(	"ATT")||triplet.equalsIgnoreCase("ATC")||triplet.equalsIgnoreCase("ATA") ){  aa='I';	}
    		if( triplet.equalsIgnoreCase(	"GGT")||triplet.equalsIgnoreCase("GGC")||triplet.equalsIgnoreCase("GGA")||triplet.equalsIgnoreCase("GGG") ){  aa='G';	}
    		if( triplet.equalsIgnoreCase(	"GCT")||triplet.equalsIgnoreCase("GCC")||triplet.equalsIgnoreCase("GCA")||triplet.equalsIgnoreCase("GCG") ){  aa='A';	}
    		if( triplet.equalsIgnoreCase(	"GTT")||triplet.equalsIgnoreCase("GTC")||triplet.equalsIgnoreCase("GTA")||triplet.equalsIgnoreCase("GTG") ){  aa='V';	}
    		if( triplet.equalsIgnoreCase(	"ACT")||triplet.equalsIgnoreCase("ACC")||triplet.equalsIgnoreCase("ACA")||triplet.equalsIgnoreCase("ACG") ){  aa='T';	}
    		if( triplet.equalsIgnoreCase(	"CCT")||triplet.equalsIgnoreCase("CCC")||triplet.equalsIgnoreCase("CCA")||triplet.equalsIgnoreCase("CCG") ){  aa='P';	}
    		if( triplet.equalsIgnoreCase(	"CTT")||triplet.equalsIgnoreCase("CTC")||triplet.equalsIgnoreCase("CTA")||triplet.equalsIgnoreCase("CTG")||triplet.equalsIgnoreCase("TTA")||triplet.equalsIgnoreCase("TTG")){  aa='L';	}
    		if( triplet.equalsIgnoreCase(	"TCT")||triplet.equalsIgnoreCase("TCC")||triplet.equalsIgnoreCase("TCA")||triplet.equalsIgnoreCase("TCG")||triplet.equalsIgnoreCase("AGT")||triplet.equalsIgnoreCase("AGC")){  aa='S';	}
    		if( triplet.equalsIgnoreCase(	"CGT")||triplet.equalsIgnoreCase("CGC")||triplet.equalsIgnoreCase("CGA")||triplet.equalsIgnoreCase("CGG")||triplet.equalsIgnoreCase("AGA")||triplet.equalsIgnoreCase("AGG")){  aa='R';	}
    		if( triplet.equalsIgnoreCase( 	"TAA")||triplet.equalsIgnoreCase("TAG")||triplet.equalsIgnoreCase("TGA")){  aa='*';	;}

	 	
    	return aa;
    }
    
    
    /**
     * translate an nucleotide sequence into an amino acid sequence. 
     * This results in 6 different reading frames: three on the forward strand and three on the reverse complementary strand.
     * 
     * 
     * @return An Array of sex FastaSequences
     */
    public FastaSequence[] translate2Protein(){
    	FastaSequence[] sixFrameTranslation = new FastaSequence[6];
    	
    	StringBuilder seq = new StringBuilder(this.getSequence());
    	
    	int nucleotideLength = seq.length();
    	
    	StringBuilder seq1 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq2 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq3 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq4 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq5 = new StringBuilder(nucleotideLength/3);
    	StringBuilder seq6 = new StringBuilder(nucleotideLength/3);
    	
    	int terminus = seq.length()%3;
    	int end=0;
    	
    	if(terminus ==0){end=-2;}
    	if(terminus ==1){end=-3;}
    	if(terminus ==2){end=-4;}   	    	
    	for (int i = 0; i< seq.length()+end; i=i+3){
    		seq1.append(this.translateTriplet(seq.charAt(i), seq.charAt(i+1), seq.charAt(i+2)));
    	}
    	sixFrameTranslation[0]= new FastaSequence(this.getIdentifier()+"_frame+0", seq1.toString());
    	//System.out.println("translated first");
    	if(terminus ==0){end=-4;}
    	if(terminus ==1){end=-2;}
    	if(terminus ==2){end=-3;}
    	for (int i = 1; i< seq.length()+end; i=i+3){
    		seq2.append(this.translateTriplet(seq.charAt(i), seq.charAt(i+1), seq.charAt(i+2)));
    	}
    	sixFrameTranslation[1]= new FastaSequence(this.getIdentifier()+"_frame+1", seq2.toString());
    	
    	if(terminus ==0){end=-3;}
    	if(terminus ==1){end=-4;}
    	if(terminus ==2){end=-2;}
    	for (int i = 2; i< seq.length()+end; i=i+3){
    		seq3.append(this.translateTriplet(seq.charAt(i), seq.charAt(i+1), seq.charAt(i+2)));
    	}
    	sixFrameTranslation[2]= new FastaSequence(this.getIdentifier()+"_frame+2", seq3.toString());
    	
    	
    	seq = new StringBuilder(this.getReverseComplementarySequence());
    	
    	if(terminus ==0){end=-2;}
    	if(terminus ==1){end=-3;}
    	if(terminus ==2){end=-4;}   	    	
    	for (int i = 0; i< seq.length()+end; i=i+3){
    		seq4.append(this.translateTriplet(seq.charAt(i), seq.charAt(i+1), seq.charAt(i+2)));
    	}
    	sixFrameTranslation[3]= new FastaSequence(this.getIdentifier()+"_frame-0", seq4.toString());
    	
    	if(terminus ==0){end=-4;}
    	if(terminus ==1){end=-2;}
    	if(terminus ==2){end=-3;}
    	for (int i = 1; i< seq.length()+end; i=i+3){
    		seq5.append(this.translateTriplet(seq.charAt(i), seq.charAt(i+1), seq.charAt(i+2)));
    	}
    	sixFrameTranslation[4]= new FastaSequence(this.getIdentifier()+"_frame-1", seq5.toString());
    	
    	if(terminus ==0){end=-3;}
    	if(terminus ==1){end=-4;}
    	if(terminus ==2){end=-2;}
    	for (int i = 2; i< seq.length()+end; i=i+3){
    		seq6.append(this.translateTriplet(seq.charAt(i), seq.charAt(i+1), seq.charAt(i+2)));
    	}
    	sixFrameTranslation[5]= new FastaSequence(this.getIdentifier()+"_frame-2", seq6.toString());
    	
    	
    	
    	
    	return sixFrameTranslation;
    }
    
    /*
    public Hashtable<Integer, FastaSequence> get6FrameTranslation(){
    	Hashtable<Integer, FastaSequence> h = new Hashtable<Integer, FastaSequence>();
    	int frame = 1;
    	char[] c = this.getSequence().toCharArray();
    	String prot = "";
    	for( int i = 0; i<c.length-2; i=i+3){
    		prot = prot + translateTriplet(""+c[i]+c[i+1]+c[i+2]);
    	}
    	h.put(frame, new FastaSequence(this.getIdentifier()+"_"+frame), this.description, prot)
    	
    	
    	
    }
    */
    public FastaSequence getLongestOpenReadingFrame(){
    	
    	int l = this.getLength();
    	
    	FastaSequence[] sixFrame = this.translate2Protein();
    	FastaSequence prot = new FastaSequence(this.identifier, "","");
    	
    	String intermediateProt = "";
    	int startIndex=-1;
    	for( int i = 0; i< sixFrame.length;i++){
    		int strand = 1;
    		int frame = i;
    		if(i>2){
    			strand = -1;
    			frame = i-3;		
    		}
    		
    		
    		FastaSequence seq = sixFrame[i];
    		char[] aa = seq.getSequence().toCharArray();
    		boolean reading = false;
    		for( int j = 0; j< aa.length; j++){
    			if(!reading && aa[j] == 'M'){
    				reading = true;
    				startIndex = j;
    				intermediateProt = intermediateProt + aa[j];
    			}else if(reading) {
    				intermediateProt = intermediateProt + aa[j];
    				if(aa[j] =='*'){
    					reading = false;
    					if(prot.getLength()<intermediateProt.length()){
    						prot.setSequence(intermediateProt );
    						
    						int startPosition = startIndex*3+frame+1;
    						int stopPosition = (j+1) *3 +frame;
    						if(strand ==-1){
    							startPosition = l - ((j+1) *3 +frame);
    							stopPosition = l-( startIndex*3 + frame);
    						}
    						
    						
    						prot.setDescription("strand="+strand+";frame="+frame+";"+startPosition+"-"+stopPosition );
    						//System.out.println(prot.getFastaString());
    					}
    					intermediateProt = "";
    					
    				}
    			}
    		}
    	}
    	return prot;
    	
    }
    
    
    /**
     * 
     * get all open reading frames (ORFs) within a the sequence above a minimum length
     * 
     * @param minLength
     * 				minimum Length of the orf
     * @return
     * 				a Vector of int arrays. Each int array has three entries: start of the orf; end of the orf and the strand.
     */
    public Vector<int[]> getOpenReadingFrames(int minLength){
    	Vector<int[]> orfs = new Vector<int[]>();
    	
    	
    	int l = this.getLength();
    	
    	FastaSequence[] sixFrame = this.translate2Protein();
    	
    	
    	int startIndex=-1;
    	for( int i = 0; i< sixFrame.length;i++){
    		int strand = 1;
    		int frame = i;
    		if(i>2){
    			strand = -1;
    			frame = i-3;		
    		}
    		
    		FastaSequence seq = sixFrame[i];
    		char[] aa = seq.getSequence().toCharArray();
    		boolean inORF = false;
    		for( int j = 0; j< aa.length; j++){
    			if(!inORF && aa[j] == 'M'){
    				inORF = true;
    				startIndex = j;
    				
    			}else if(inORF) {
    				
    				if(aa[j] =='X'){
    					inORF = false;
    					
    					int startPosition = startIndex*3+frame+1;
						int stopPosition = (j+1) *3 +frame;
						if(strand ==-1){
							startPosition = l - ((j+1) *3 +frame) +1;
							stopPosition = l-( startIndex*3 + frame);
						}
    					
    					if( stopPosition+1-startPosition >= minLength){
    						int[] orf = {startPosition,stopPosition,strand };
    						orfs.add(orf);
    					}
    					
    					
    				}
    			}
    		
    		}
    	}
    	
    	return orfs;
    }
    
    
    
    /**
     * get an int[] with the positions of occurance of a subsequence.
     * 
     * @param subString 
     * 				the motiv that will be searched for in the FastaSequence
     * @return
     * 				the int[] with all start positions of the motiv in the sequence.
     */
    public int[] getSubSequencePositions(String motive){
    	
    	int[] positions = new int[0];
    	int fromIndex = 0;
    	int index = 0;
    	while(index != -1){
    		index = sequence.indexOf(motive, fromIndex);
    		fromIndex = index+1;
    		if(index != -1){
    			int[] tmpPositions = new int[positions.length+1];
    			for (int i = 0; i< positions.length; i++){
    				tmpPositions[i]=positions[i];
    			}
    			tmpPositions[positions.length]=index;
    			positions=tmpPositions ;
    		}
    	}
    	return positions;
    }
    
    
  
    
    
    
    
 
   
    
    
 
    
    public void maskToLowerCase(int maskStart, int maskEnd){
    	char[] a = this.getSequence().toCharArray();
    	
    	for( int i = maskStart-1; i< maskEnd; i++){
    		a[i] = Character.toLowerCase(a[i]);
    	}
    	String s = new String(a);
    	
    	this.setSequence(s);
    }
    
    
    
}