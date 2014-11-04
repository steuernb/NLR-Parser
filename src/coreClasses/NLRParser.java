package coreClasses;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

public class NLRParser {
	static final double version = 0.12;
	
	public static double getVersion(){
		return version;
	}
	
	
	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Options options = new Options();
		options.addOption( OptionBuilder.withLongOpt( "input")
										.withDescription("The XML output from a MAST run")
										.isRequired()
										.hasArg()
										.withArgName("MastOutput.xml")
										.create('i'));
		
		options.addOption( OptionBuilder.withLongOpt( "output")
										.withDescription("TSV output File")
										.hasArg()
										.isRequired()
										.withArgName("outputFile")
										.create('o'));
		options.addOption( OptionBuilder.withLongOpt( "splitpattern")
										.withDescription("All different entries that share the prefix before the first occurence of this pattern will be merged for analysis. Default is _frame")
										.hasArg()
										.withArgName("splitPattern")
										.create('s'));	
		options.addOption(OptionBuilder.withLongOpt("pValue")
									   .withDescription("p-value threshold for each individual motif.")
									   .hasArg()
									   .withArgName("pvalue")
									   .create('p'));
		options.addOption(OptionBuilder.withLongOpt("blastFile")
				   					   .withDescription("A blast file of sequences against a reference protein NBLRR set. The format has to be xml (-outfmt 5) and the query has to be the same inputfile as for the mast")
				   					   .hasArg()
				   					   .withArgName("blastFile")
				   					   .create('b'));
		options.addOption(OptionBuilder.withLongOpt("sequenceFile")
				   						.withDescription("A sequence file with the AA sequences for each protein submitted to MAST")
				   						.hasArg()
				   						.withArgName("aasequence.fasta")
				   						.create('a'));
		options.addOption(OptionBuilder.withLongOpt("writeGFF")
										.withDescription("Output is general feature format")
										.create('g'));
		options.addOption(OptionBuilder.withLongOpt("help")
				   						.withDescription("display this help")
				   						.create('h'));
		
		
		CommandLineParser parser = new PosixParser();
		
		
		
		try {
			
			CommandLine line = parser.parse( options, args );
	        
			if(line.hasOption('h')){
				HelpFormatter formatter = new HelpFormatter();
		        formatter.printHelp( "MastParser", options ); 
			}else{
			
			
				File mastFile   =    new File(line.getOptionValue('i'));
		        File outputFile =    new File(line.getOptionValue('o'));
		        
		        
		        
		        double pvalue = 1e-5;
		        if(line.hasOption('p')){
		        	try{
		        		pvalue = Double.parseDouble(line.getOptionValue('p'));
		        	}catch(NumberFormatException e){
		        		throw new ParseException ("Value for option -p has to be a number");
		        	}
		        }
		        
		        
		        
		       
		        String splitPattern = "_frame";
		        if (line.hasOption('s')){
		        	splitPattern = line.getOptionValue('s');
		        }
		       
		        MastFile file = new MastFile(mastFile, splitPattern, pvalue);
		        
		        if(line.hasOption('b')){
		        	File blastFile = new File(line.getOptionValue('b'));
		        	file.setBlast(blastFile);
		        }
		        
		        if(line.hasOption('a')){
		        	file.addSequences(new File(line.getOptionValue('a')));
		        }
		        
		        
		        file.mergeEntries();
		        
		        if(line.hasOption('g')){
		        	file.writeGFF(outputFile);
		        }else{
		        	file.writeOutput(outputFile);
		        }
		        
		        	
			}    	
	        	
	        
		}catch( ParseException exp ) {
	       
	        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
	        HelpFormatter formatter = new HelpFormatter();
	        formatter.printHelp( "java -jar NLR-Parser.jar -i mast.xml -o output.txt", options ); 
		}   
		catch(IOException e){
			e.printStackTrace();
		}
	}

}
