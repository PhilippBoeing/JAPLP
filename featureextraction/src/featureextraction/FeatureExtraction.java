package featureextraction;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.LinkedList;
import org.biojava3.aaproperties.PeptideProperties;

public class FeatureExtraction {
	
	  
	public class FastaElement {
		public String id;
		public String sequence;
		
		public FastaElement(String id, String sequence){
			this.id = id;
			this.sequence = sequence;
		}
	}
	
	private int featureIndex = 0;
	
	// expected input: path to a directory with .fasta
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		if(args.length != 1){
			System.out.println("Error - exactly one argument needed.\nProvide path to folder containing .fasta input data files.");
			}
		else{
			@SuppressWarnings("unused")
			FeatureExtraction fe = new FeatureExtraction(args[0]);
		}

	}
	
	public FeatureExtraction(String inputFolder){
		File directory = new File(inputFolder);
		if(!(directory.exists() && directory.isDirectory())){
			System.out.println("Error - argument must be a folder");
			return;
		}
		
		File[] files = directory.listFiles();
		BufferedWriter writer = null;
		try {
			writer = new BufferedWriter(new FileWriter(directory.getAbsoluteFile() + "/features"));
		} catch (IOException e1) {
			System.out.println("Error writing to file.");

		}

		HashMap<String,Integer> labelling = new HashMap<String,Integer>();
		int label = 1;
		for(File file : files){
			if(file.getName().endsWith(".fasta")){
				System.out.println("Extracting features from "+file.getName());
				try {
					extractFeatures(file,writer, label);
					labelling.put(file.getName(), label);
					label++;
				} catch (IOException e) {
					System.out.println("Error reading the file.");
				}
			}
		}
		
		System.out.println("Labelling:");
		for(String key : labelling.keySet()){
			System.out.println(key + ": " + Integer.toString(labelling.get(key)));
		}
		try {
			writer.close();
		} catch (IOException e) {
			System.out.println("Error writing to file.");
		}
		System.out.println("Done.");
	}
	
	private String getFeatureIndex(){
		return ' '+Integer.toString(++featureIndex)+':';
	}
	
	private void resetFeatureIndex(){
		featureIndex = 0;
	}

	private void extractFeatures(File file,BufferedWriter writer, int label) throws IOException {
		LinkedList<FastaElement> items = getFastaElements(file);
		
		String outputLine = "";
		int counter = 0;
		
		for(FastaElement item : items){
			++counter;
			
			// this can be used to limit the number of samples to create a reduced training set. 
			if (counter > 150)
				break;
			
			writer.write(Integer.toString(label));
			resetFeatureIndex();
			
			
			// each of these features can be turned off by commenting one of these lines
			
			// Sequence Length
			writer.write(getFeatureIndex()+Integer.toString(item.sequence.length()));
			
			// Global Amino Acid Composition
			writer.write(getAminoComposition(item.sequence));
			
			// Local (first 50) Amino Acid Composition
			writer.write(getAminoComposition(item.sequence.substring(0,Math.min(50, item.sequence.length()))));
			
			// Local (last 50) Amino Acid Composition
			writer.write(getAminoComposition(item.sequence.substring(Math.max(0, item.sequence.length()-50))));			
			
			// Molcular Weight
			writer.write(getFeatureIndex()+Double.toString(PeptideProperties.getMolecularWeight(item.sequence)));
			
			// Isoelectric Point
			writer.write(getFeatureIndex()+Double.toString(PeptideProperties.getIsoelectricPoint(item.sequence)));
			
			// Apliphatic Index
			writer.write(getFeatureIndex()+Double.toString(PeptideProperties.getApliphaticIndex(item.sequence)));
			
			// Average Hydropathy
			writer.write(getFeatureIndex()+Double.toString(PeptideProperties.getAvgHydropathy(item.sequence)));
			
			// Instability Index
			writer.write(getFeatureIndex()+Double.toString(PeptideProperties.getInstabilityIndex(item.sequence)));
			
			writer.write('\n');
			System.out.print(item.id + ", ");
		}
		
		System.out.print('\n');
	}

	private String getAminoComposition(String sequence) {
		Character[] aminoAcids = {'G','P','A','V','L','I','M','C','F','Y','W','H','K','R','Q','N','E','D','S','T'};
		// we will ignore the additonal U,B,Z,X of FASTA, since only two inputs use them.

		HashMap<Character,Integer> map = new HashMap<Character,Integer>();
		

		
		for(int i = 0; i<sequence.length();i++){
			Integer count = map.get(sequence.charAt(i));
			if(count == null){
				map.put(new Character(sequence.charAt(i)), 1);
			}
			else{
				map.put(new Character(sequence.charAt(i)), count+1);
			}
		}
		
		int totalCount = 0;
		
		for(Character aminoAcid : aminoAcids){
			Integer count =  map.get(aminoAcid);
			if(count != null){
				totalCount  += count;
			}
		}
		
		String result = "";
		for(Character aminoAcid : aminoAcids){
			Integer count =  map.get(aminoAcid);
			if(count == null){
				result += getFeatureIndex()+"0";
			}
			else{
				result += getFeatureIndex()+Double.toString((double) count / totalCount);
			}
			
			
		}
		
		return result;
	}

	private LinkedList<FastaElement> getFastaElements(File file) throws IOException {
		BufferedReader reader = new BufferedReader(new FileReader(file));
		LinkedList<FastaElement> resultList = new LinkedList<FastaElement>();
		
		String line = reader.readLine();
		String id = "";
		String seq = "";
		while(true){
			if(line == null || line.startsWith(">")){
				if(id != ""){
					resultList.add(new FastaElement(id,seq));
				}
				if(line == null){
					break;
				}
				seq = "";
				id = line.substring(1); 
			}
			else{
				seq = seq + line;
			}
			line = reader.readLine();
		}
		
		reader.close();
		return resultList;
	}

}
