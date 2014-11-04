package addonPrograms;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;

import supportClasses.FastaReader;
import supportClasses.FastaSequence;


/**
 * 
 * @author steuernb
 * @version 0.3
 */
public class Translate6Frame {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Options options = new Options();
		options.addOption( OptionBuilder.withLongOpt("input")
										.withDescription("Input File in Fasta format")
										.hasArg()
										.isRequired()
										.withArgName("input")
										.create('i'));
		options.addOption(OptionBuilder.withLongOpt("output")
										.withDescription("outputFile containting filtered sequence set")
										.hasArg()
										.create('o'));
		
		options.addOption("h" ,"help", false, "Print this help");
		
		
		
		CommandLineParser parser = new PosixParser();
		
		 try {
		        // parse the command line arguments
		        CommandLine line = parser.parse( options, args );
		        
		      
		        String input = line.getOptionValue('i');
		        
		        String output = input + ".6frame.fasta";
		        if(line.hasOption('o')){
		        	output = line.getOptionValue('o');
	        	}
		        
		       
		        
		        File inputFile = new File(input);
		        File outputFile = new File(output);
		      
		       try {
		    	   translate( inputFile,  outputFile);
			} catch (Exception e) {
				e.printStackTrace();
				throw new ParseException ("");
			} 
		       
		        
		        
		    }
		    catch( ParseException exp ) {
		       
		        System.err.println( "Parsing failed.  Reason: " + exp.getMessage() );
		        HelpFormatter formatter = new HelpFormatter();
		        formatter.printHelp( "Translate nucleotide sequence to 6 reading frames", options );
		    }
		

	}

	public static void translate(File inputFile, File outputFile)throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputFile));
		FastaReader reader = new FastaReader(inputFile);
		FastaSequence seq = reader.readEntry();
		while (seq != null) {
			FastaSequence[] a = seq.translate2Protein();
			for( int i = 0; i< a.length; i++){
				out.write(a[i].getFormatedFastaString(100));
			}
			seq = reader.readEntry();
		}
		reader.close();
		out.close();
		
		
	}
	
}
