package iHMMuneAlign;

import java.io.*;
import java.util.ArrayList;
import org.biojava.bio.Annotation;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.*;
import org.biojava.bio.seq.DNATools;
//import org.biojava.bio.seq.Sequence;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojavax.bio.seq.RichSequence;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.dp.SimpleMarkovModel;

import java.util.Iterator;
import java.util.Random;

// Referenced classes of package iHMMuneAlign:
//            GlobalDefines, ExponentialDecay, MutationSpectrum, FastaReader, 
//            ProbabilityHolder, StateInfo, FakeAnnotation, MutabilityScore

public class VnDnJCnoC
    implements GlobalDefines
{

    final int GENE_TYPE_V = 1;
    final int GENE_TYPE_D = 2;
    final int GENE_TYPE_J = 3;
    boolean DEBUG;
    double MIN_PROB_LIMIT;
    MutationSpectrum G_mutationSpectrum;
    MutabilityScore G_mutability_score;
    ExponentialDecay G_exponential_decay;
    PrintWriter G_fostream;
    double toMagicalStateProb;
    boolean DEBUGGING = false;

	//array lists for vdj states
	ArrayList vStates = new ArrayList();
	ArrayList dStates = new ArrayList();
	ArrayList jStates = new ArrayList();
	//array lists for vdj distributions
	ArrayList vStateDistributions = new ArrayList();
	ArrayList dStateDistributions = new ArrayList();
	ArrayList jStateDistributions = new ArrayList();
	DotState X1a = null; // end of V when no exo
	DotState X1b = null; // end of V when exo
	DotState X2 = null; // start of N1
	DotState X3 = null; // end of N1
	DotState X8 = null;
	DotState X9 = null; 
	DotState X5a_states[] = null;
	DotState X5b_states[] = null;
	DotState X7a_states[] = null; 
	DotState X7b_states[] = null;
	DotState X11a_states[] = null;
	DotState X11b_states[] = null;
	EmissionState VD_N_states[] = null;
	EmissionState DJ_N_states[] = null;
	EmissionState Kappa_N_states[] = null;
	EmissionState Lambda_N_states[] = null;
	SimpleMarkovModel vdj = null;
	int chainType;
	ProbabilityHolder probHolder;
	boolean removed_C_Region;
    private static PrintStream v_out;

    public VnDnJCnoC(MutabilityScore mutabilityScore, String mutationSpecFileName)
    {
        DEBUG = false;
        MIN_PROB_LIMIT = 9.9999999999999995E-007D;
        G_mutability_score = mutabilityScore;
        G_exponential_decay = new ExponentialDecay();
        G_mutationSpectrum = new MutationSpectrum(mutationSpecFileName);//new MutationSpectrum("Mutation_spectrum.txt");

        try
        {
            Random rnd = new Random();
            String G_fostreamFileName = "gene_nucleotide_mutation_probabilities" + rnd.nextInt() + ".txt";
            G_fostream = new PrintWriter(new FileWriter(G_fostreamFileName));
            
            File G_fostreamFile = new File(G_fostreamFileName);
            G_fostreamFile.deleteOnExit();
        }
        catch(IOException ioex)
        {
            throw new Error("output stream could not be opened");
        }
    }

    public MarkovModel createModel(RichSequence vSequence, int VGene_start_offset, int completeVGeneLength, int from_alignment_start_UMS_length, ArrayList dSequences, ArrayList jSequences, ProbabilityHolder probHolder, boolean removed_C_Region, 
            double A_probability, int alignmentModelType, int chainType)
    {
        FiniteAlphabet dna = DNATools.getDNA();
        //FastaReader fr = new FastaReader();
        int advance[] = {1};
		this.probHolder = probHolder;
		this.chainType = chainType;
		this.removed_C_Region = removed_C_Region;
        
        //set the to magical state transition based on the alignmentModelType
        toMagicalStateProb = 0.0D;
        if (DEBUGGING) {
        	System.out.println("toMagicalState prob set to " + toMagicalStateProb + " alignmentModelType is " + alignmentModelType);
        }
        if (alignmentModelType == 2) {
           if (from_alignment_start_UMS_length > 0) {
        	  double from_alignment_start_VGene_length = completeVGeneLength - VGene_start_offset - 1;
        	  toMagicalStateProb = (from_alignment_start_VGene_length/from_alignment_start_UMS_length - 0.75) * 4 * 0.5; //0.02D;
              // ensure toMagicalStateProb is not negative
        	  if (toMagicalStateProb < 0) {
            	  toMagicalStateProb = 0.0D;
              }
        	  
              if (!removed_C_Region) {
                 //toMagicalStateProb = 0.05D;
                 if (DEBUGGING) {
                 	System.out.println("toMagicalState prob set to " + toMagicalStateProb);
                 }
              }
           }
        }
		if (DEBUGGING) {
			System.out.println("toMagicalState prob set to " + toMagicalStateProb);
		}
        // do the following long list of stuff
        // anything to do with D states, put if (chainType == 1) around them
        
		
        /*try
        {
            if (chainType == 1) {// heavy chain
        	    dSequences = fr.readFile(D_GENES_FILE);
				System.out.println("number of D germline sequences read from file = " + dSequences.size());
            }
            jSequences = fr.readFile(J_GENES_FILE);
			System.out.println("number of J germline sequences read from file = " + jSequences.size());
            System.out.println("IGHV germline sequence found from alignment has length = " + vSequence.length());
            
        }
        catch(Exception e)
        {
            System.out.println("fr.readFile(): " + e.getMessage() + D_GENES_FILE + " not a fasta file ???");
            System.exit(0);
        }*/
        //create states to represent the start and end of different regions
        //these are non-emitting states - just place holders
        X1a = createDotState("X1a"); // end of V when no exo
        X1b = createDotState("X1b"); // end of V when exo
        X2 = createDotState("X2"); // start of N1
        X3 = createDotState("X3"); // end of N1
        if (chainType == 1) {// heavy chain
        	X8 = createDotState("X8"); // start of N2
            X9 = createDotState("X9"); // end of N2
	        X5a_states = createDotStateArray("X5a", dSequences.size()); //start of D when no exo - for each germline D
	        X5b_states = createDotStateArray("X5b", dSequences.size()); // start of D when exo
	        X7a_states = createDotStateArray("X7a", dSequences.size()); // end of D with no exo
	        X7b_states = createDotStateArray("X7b", dSequences.size()); // end of D with exo
        }
	    
		X11a_states = createDotStateArray("X11a", jSequences.size()); //start of J with no exo - for each germline J
	    X11b_states = createDotStateArray("X11b", jSequences.size()); //start of J with exo
        
        // create states for each nucleotide in the germline gene found to be the best match to the input sequence
        createVStates(advance, dna, vSequence, vStates, vStateDistributions, VGene_start_offset, completeVGeneLength, probHolder, A_probability);
        //creates states for each of the germline D genes
        if (chainType == 1) {// heavy chain
    	    createDStates(advance, dna, dSequences, dStates, dStateDistributions, probHolder, completeVGeneLength, A_probability);
        }
        //creates states for each of the germline J genes
        createJStates(advance, dna, jSequences, jStates, jStateDistributions, probHolder, completeVGeneLength, A_probability);

        //create the states for the N1 and N2 regions
        if (chainType == 1) {// heavy chain
        	VD_N_states = createNStates(advance, dna, "V-D N-Region", true, probHolder.VD_N_Fields.length);
            DJ_N_states = createNStates(advance, dna, "D-J N-Region", false, probHolder.DJ_N_Fields.length);
        } else if (chainType == 2) { // kappa chain
        	Kappa_N_states = createNStates(advance, dna, "Kappa N-Region", true, probHolder.Kappa_N_Fields.length);
        } else if (chainType == 3) { // lambda chain
        	Lambda_N_states = createNStates(advance, dna, "Lambda N-Region", true, probHolder.Lambda_N_Fields.length);
        }

        // create a SMM object to hold the model
		vdj = new SimpleMarkovModel(1, dna, "VDJGene");
        
        //add the states that define the start and end of the V (with and without exo), the N1 and N2
        if (chainType == 1) {//heavy chain
        	addDotStates(vdj, X1a, X1b, X2, X3, X8, X9);
        } else {
        	addDotStates(vdj, X1a, X1b, X2, X3);
        }
        
        if (chainType == 1) {//heavy chain
        	//add the states from the arrays for start of the IGHD without and with exonuclease
	        addDotStateArray(vdj, X5b_states);
	        addDotStateArray(vdj, X5a_states);
	        //add the states from the arrays for the end of the IGHD with and without exonuclease
	        addDotStateArray(vdj, X7a_states);
	        addDotStateArray(vdj, X7b_states);
        }
        
        //add the states from the arrays for the start of the IGHJ without and with exonuclease
        addDotStateArray(vdj, X11b_states);
        addDotStateArray(vdj, X11a_states);

        //add the N-regions states into the model
        if (chainType == 1) {// heavy chain
        	addNStates(VD_N_states, vdj);
        	addNStates(DJ_N_states, vdj);
        } else if (chainType == 2) { //kappa chain
         	addNStates(Kappa_N_states, vdj);
        } else if (chainType == 3) { //lambda chain
        	addNStates(Lambda_N_states, vdj);
        }

        //add the V states (just one sequences with per nucleotide emissions) - so each nt gets a prob
        //germline nt, then each of the other three based on potential of mutation to that nt
        addVStates(vdj, vStates);

        Utils.writeTo("v_probs.txt", "Start_Model");
        printStates(vSequence, vStates);
        if (chainType == 1) {// heavy chain
        	//add the D states (for all specified germline sequences from repertoire file) - each nt in a given seq gets its own state with emm dist set
        	addDStates(vdj, dStates);
            for (Object sequence : dStates) {
                printStates(vSequence, (ArrayList) sequence);
            }
        }
        //add the J states (for all specified germline sequences from repertoire files - each nt in a given J gets own state based on emm data
        addJStates(vdj, jStates);
        for (Object sequence : jStates) {
            printStates(vSequence, (ArrayList) sequence);
        }

        /*
        // the V end exo counts - assuming this means that there were 145 sequence with 0, 82 with 1 nt lost, 32 with 2 and so on
        // not sure why these are here rather than in the probability holder like everything else
        int V_end_exo_count[] = {
            145, 82, 32, 21, 7, 12, 6, 1, 1, 1
        };*/

        //create the transitions between the states for the V gene
        createVTransitions(vdj, vStates, X1a, X1b, probHolder, chainType);

        if (chainType == 1) {// heavy chain
	        //create transitions for the end of the V and the N regions
	        create_VD_N_part_1_Transitions(vdj, VD_N_states, X1a, X1b, X2, X3);
	        //create transitiong for the end of the N and the start of the D (with and without exo)
	        create_VD_N_part_2_Transitions(vdj, X3, X5a_states, X5b_states);
	        
	        //create the transits for the D
	        createDTransitions_multi(vdj, X5b_states, X5a_states, dStates, X7a_states, X7b_states, probHolder, X8);
	
	        //create the D-N2-J regions similar to N1 function, but coming from multiple D rather than single IGHV
	        //this does the end of the D and the N
	        create_VD_N_part_1_Transitions_multi(vdj, DJ_N_states, X7a_states, X7b_states, X8, X9);
	        //this does the end of the N and the start of the J (with and without exo)
	        create_VD_N_part_2_Transitions(vdj, X9, X11a_states, X11b_states);
        } else if (chainType == 2) {// kappa chain
        	create_VD_N_part_1_Transitions(vdj, Kappa_N_states, X1a, X1b, X2, X3);
        	create_VD_N_part_2_Transitions(vdj, X3, X11a_states, X11b_states);
	    } else if (chainType == 3) {// lambda chain
	    	create_VD_N_part_1_Transitions(vdj, Lambda_N_states, X1a, X1b, X2, X3);
	    	create_VD_N_part_2_Transitions(vdj, X3, X11a_states, X11b_states);
	    }
        
        //creates the J transitions from the with and without exo j start to the end of the markov model
        createJTransitions(vdj, X11a_states ,X11b_states, jStates, probHolder, chainType, removed_C_Region);

        /*
        //probs for 0, 1, 2 removals from the IGHV gene end
        //once again unsure why these are included here as should come from the prob holder
        double V_end_exo_probs[] = {
            0.10000000000000001D, 0.20000000000000001D, 0.40000000000000002D
        };*/
        
        //set the V transitions
        setVTransitions(vdj, vStates, X1a, X1b, probHolder, chainType);

        if (chainType == 1) {// heavy chain {
	        //set the trans probs for the first part of the V-N1-D
	        set_VD_N_part_1_Transitions(vdj, VD_N_states, probHolder.VD_N_Fields, X1a, X1b, X2, X3);
	        //set the probs for the rest of the V-N-D, the N to D part
	        set_VD_N_part_2_Transitions_multi(vdj, X3, X5b_states);
	        
	        //set the D trans probs
	        setDTransitions_multi(vdj, X5b_states, X5a_states, dStates, X7a_states, X7b_states, probHolder, X8);
	
	        //now repeat the same as for the V-N1-D, but this time using the D-N2-J
	        set_VD_N_part_1_Transitions_multi(vdj, DJ_N_states, probHolder.DJ_N_Fields, X7a_states, X7b_states, X8, X9);
	        set_VD_N_part_2_Transitions_multi(vdj, X9, X11b_states);
        } else if (chainType ==2) { //kappa
        	set_VD_N_part_1_Transitions(vdj, Kappa_N_states, probHolder.Kappa_N_Fields, X1a, X1b, X2, X3);
	        set_VD_N_part_2_Transitions_multi(vdj, X3, X11b_states);
	    } else if (chainType == 3) { //lambda
	    	set_VD_N_part_1_Transitions(vdj, Lambda_N_states, probHolder.Lambda_N_Fields, X1a, X1b, X2, X3);
	        set_VD_N_part_2_Transitions_multi(vdj, X3, X11b_states);
	    }
        
        //finally set the J trans probs and probs for ending the model
        setJTransitions(vdj, X11b_states, X11a_states, jStates, probHolder, chainType, removed_C_Region);

        G_fostream.close();
		
		//return the model
        return vdj;
    }

    private static void printStates(RichSequence vSequence, ArrayList vStates) {
        if (v_out == null) {
            v_out = Utils.getWriter("v_probs.txt");
        }
        Utils.writeTo("v_names.txt", vSequence.getName());
        for (Object o : vStates) {
            SimpleEmissionState s = (SimpleEmissionState) o;
            FiniteAlphabet a = (FiniteAlphabet) s.getDistribution().getAlphabet();
            for (Iterator i = a.iterator(); i.hasNext();) {
                Symbol sym = (Symbol) i.next();
                try {
                    v_out.print(s.getDistribution().getWeight(sym) + " ");
                } catch (IllegalSymbolException e) {
                    e.printStackTrace();
                }
            }
            v_out.println();
        }
    }

    private void addVStates(SimpleMarkovModel vdj, ArrayList vStates)
    {
        try
        {
            for(int i = 0; i < vStates.size(); i++)
            {
                vdj.addState((EmissionState)vStates.get(i));
            }

        }
        catch(Exception e)
        {
            throw new Error("Can't add v states to model");
        }
    }

    private void addEmissionState(SimpleMarkovModel vdj, EmissionState emission_state)
    {
        try
        {
            vdj.addState(emission_state);
        }
        catch(Exception e)
        {
            throw new Error("addEmissionState():  " + e.getMessage());
        }
    }
	
	private	void removeEmissionState(SimpleMarkovModel vdj, EmissionState emission_state) {
		try {
			vdj.removeState(emission_state);
		}
		catch (Exception e) {
			throw new Error("removeEmissionState(): " + e.getMessage());
		}
	
	}

    private void addDStates(SimpleMarkovModel vdj, ArrayList dStates)
    {
        try
        {
            for(int i = 0; i < dStates.size(); i++)
            {
                ArrayList arraylist = (ArrayList)dStates.get(i);
                for(int j = 0; j < arraylist.size(); j++)
                {
                    vdj.addState((EmissionState)arraylist.get(j));
                }

            }

        }
        catch(Exception e)
        {
            throw new Error("Can't add D states to model");
        }
    }

    private void addJStates(SimpleMarkovModel vdj, ArrayList jStates)
    {
        try
        {
            for(int i = 0; i < jStates.size(); i++)
            {
                ArrayList arraylist = (ArrayList)jStates.get(i);
                for(int j = 0; j < arraylist.size(); j++)
                {
                    vdj.addState((EmissionState)arraylist.get(j));
                }

            }

        }
        catch(Exception e)
        {
            throw new Error("Can't add J states to model");
        }
    }

    //sets the emission probs for the N1 region
    private void set_VD_Emission(Distribution nDist)
    {
        try
        {
            nDist.setWeight(DNATools.a(), 0.23599999999999999D);
            nDist.setWeight(DNATools.c(), 0.2387D);
            nDist.setWeight(DNATools.g(), 0.37130000000000002D);
            nDist.setWeight(DNATools.t(), 0.15379999999999999D);
        }
        catch(Exception e)
        {
            throw new Error("set_VD_Emission: Can't set State Emission");
        }
    }

    //sets the emission probs for the N2 region
    private void set_DJ_Emission(Distribution nDist)
    {
        try
        {
            nDist.setWeight(DNATools.a(), 0.2079D);
            nDist.setWeight(DNATools.c(), 0.3135D);
            nDist.setWeight(DNATools.g(), 0.33000000000000002D);
            nDist.setWeight(DNATools.t(), 0.14849999999999999D);
        }
        catch(Exception e)
        {
            throw new Error("set_DJ_Emission: Can't set State Emission");
        }
    }



    private DotState[] createDotStateArray(String name, int size)
    {
        DotState dotStateArray[] = new DotState[size];
        for(int i = 0; i < size; i++)
        {
            StateInfo stateInfo = new StateInfo("Dot State " + name + i, 'X');
            FakeAnnotation annotation = new FakeAnnotation(stateInfo);
            dotStateArray[i] = new SimpleDotState(name, annotation);
        }

        return dotStateArray;
    }

    private DotState createDotState(String name)
    {
        StateInfo stateInfo = new StateInfo("Dot State " + name, 'X');
        FakeAnnotation annotation = new FakeAnnotation(stateInfo);
        return new SimpleDotState(name, annotation);
    }

    private void addDotStateArray(MarkovModel vdj, DotState dotStateArray[])
    {
        try
        {
            for(int i = 0; i < dotStateArray.length; i++)
            {
                vdj.addState(dotStateArray[i]);
            }

        }
        catch(Exception e)
        {
            throw new Error("addDotStateArray(): " + e.getMessage());
        }
    }

    private void addDotStates(MarkovModel vdj, DotState X1a, DotState X1b, DotState X2, DotState X3, DotState X8, DotState X9)
    {
        try
        {
            vdj.addState(X1a);
            vdj.addState(X1b);
            vdj.addState(X2);
            vdj.addState(X3);
            vdj.addState(X8);
            vdj.addState(X9);
        }
        catch(Exception e)
        {
            throw new Error("addDotStates: " + e.getMessage());
        }
    }

    private void addDotStates(MarkovModel vdj, DotState X1a, DotState X1b, DotState X2, DotState X3)
    {
        try
        {
            vdj.addState(X1a);
            vdj.addState(X1b);
            vdj.addState(X2);
            vdj.addState(X3);
        }
        catch(Exception e)
        {
            throw new Error("addDotStates: " + e.getMessage());
        }
    }
	
	private void removeDotStates(MarkovModel vdj, DotState X1a, DotState X1b, DotState X2, DotState X3) {
		try {
			vdj.removeState(X1a);
			vdj.removeState(X1b);
			vdj.removeState(X2);
			vdj.removeState(X3);
		} catch (Exception e) {
			throw new Error("removeDotStates: " + e.getMessage());
		}
	}


    private EmissionState[] createNStates(int advance[], FiniteAlphabet dna, String geneName, boolean VD, int size)
    {
        try
        {
            EmissionState nStates[] = new EmissionState[size];
            for(int i = 0; i < nStates.length; i++)
            {
                Distribution tempDist = DistributionFactory.DEFAULT.createDistribution(dna);
                char token;
                if(VD)
                {
                    set_VD_Emission(tempDist);
                    token = '1';
                } else
                {
                    set_DJ_Emission(tempDist);
                    token = '2';
                }
                StateInfo stateInfo = new StateInfo(geneName, token);
                FakeAnnotation annotation = new FakeAnnotation(stateInfo);
                nStates[i] = new SimpleEmissionState(geneName, annotation, advance, tempDist);
            }

            return nStates;
        }
        catch(Exception e)
        {
            throw new Error("createNState(): Can't create N distributions");
        }
    }

    private void addNStates(EmissionState n_states[], MarkovModel vdj)
    {
        for(int i = 0; i < n_states.length; i++)
        {
            try
            {
                vdj.addState(n_states[i]);
            }
            catch(Exception ex)
            {
                throw new Error("addNStates():  " + ex.getMessage());
            }
        }

    }

    private void createVStates(int advance[], FiniteAlphabet dna, RichSequence vSequence, ArrayList vStates, ArrayList vStateDistributions, int VGene_start_offset, int completeVGeneLength,
            ProbabilityHolder probHolder, double A_probability)
    {
        try
        {
            Symbol sym = null;
            double totalMutationProb = probHolder.Gene_Mutation_Field;
            
			//get the germline IGHV sequence and name
			String seq_string = vSequence.seqString();
			String geneName = vSequence.getName();
			//for each nt from the portion of germline IGHV that matches to query/input sequence build the state information
			for(int i = 1; i <= vSequence.length(); i++) {
                //create a DNA distribution for the current position
				Distribution tempDist = DistributionFactory.DEFAULT.createDistribution(dna);
                vStateDistributions.add(tempDist);
                //set an emission probability for current position
				setGeneEmissionProb(vSequence, seq_string, i, tempDist, dna, 1, completeVGeneLength, VGene_start_offset, A_probability);
                //set the state info for current position
				StateInfo stateInfo = new StateInfo(geneName, 'V', -1, i + VGene_start_offset, getChar(vSequence.symbolAt(i)), completeVGeneLength);
                FakeAnnotation annotation = new FakeAnnotation(stateInfo);
                vStates.add(new SimpleEmissionState(geneName, annotation, advance, tempDist));
				

				//clean-up before exiting
				/*tempDist = null;
				stateInfo = null;
				annotation = null;*/
            }
			//clean-up
			//seq_string = null;
			//geneName = null;
        }
        catch(Exception e)
        {
            e.printStackTrace();
            throw new Error("createVStates: Can't create distributions(): " + e.getMessage());
        }
    }

    private void createDStates(int advance[], FiniteAlphabet dna, ArrayList dSequences, ArrayList dStates, ArrayList dStateDistributions, ProbabilityHolder probHolder, int completeVGeneLength,
            double A_probability)
    {
        double totalMutationProb = probHolder.Gene_Mutation_Field;
        try
        {
            Symbol sym = null;
			//for each germline IGHD gene in the repertoire
            for(int i = 0; i < dSequences.size(); i++)
            {
                //get sequence from repertoire
				RichSequence seq = (RichSequence)dSequences.get(i);
				//get sequence name & string
                String seq_string = seq.seqString();
				String geneName = seq.getName();
		
				//create arraylists to store the state and distribution details for current IGHD sequence
                ArrayList stateArrayList = new ArrayList();
                ArrayList distArrayList = new ArrayList();
				
				//for each nt in the current germline IGHD sequence
				for(int j = 1; j <= seq.length(); j++)
                {
                    //create a temp distribution to be associated with current nt
					Distribution tempDist = DistributionFactory.DEFAULT.createDistribution(dna);
                    distArrayList.add(tempDist);
					//set the emission probs for the current nt
                    setGeneEmissionProb(seq, seq_string, j, tempDist, dna, 2, completeVGeneLength, 0, A_probability);
                    //save the state info/emission probs for the current nt
					StateInfo stateInfo = new StateInfo(geneName, 'D', i, j, getChar(seq.symbolAt(j)), seq.length());
                    FakeAnnotation annotation = new FakeAnnotation(stateInfo);
                    stateArrayList.add(new SimpleEmissionState(geneName, annotation, advance, tempDist));
					
					//clean-up before exiting
					/*tempDist = null;
					geneName = null;
					stateInfo = null;
					annotation = null;*/
                }
				//for the current IGHD gene add nt states to array list
                dStates.add(stateArrayList);
				//for current IGHD add nt distributions to array list
                dStateDistributions.add(distArrayList);
				
				//clean-up before moving to next IGHD gene in repertoire
				/*stateArrayList = null;
				distArrayList = null;
				seq_string = null;
				geneName = null;*/
            }
			
        }
        catch(Exception e)
        {
            throw new Error("createDStates: Can't create distributions");
        }
    }

    private char getChar(Symbol symbol)
    {
        return symbol.getName().charAt(0);
    }

    private void createJStates(int advance[], FiniteAlphabet dna, ArrayList jSequences, ArrayList jStates, ArrayList jStateDistributions, ProbabilityHolder probHolder, int completeVGeneLength,
            double A_probability)
    {
        double totalMutationProb = probHolder.Gene_Mutation_Field;
        try
        {
            Symbol sym = null;
			//for each of the IGHJ genes in the germline repertoire
            for(int i = 0; i < jSequences.size(); i++)
            {
                //get the sequence string and name from the RichSequence object
				RichSequence seq = (RichSequence)jSequences.get(i);
                String seq_string = seq.seqString();
				String geneName = seq.getName();
				
				//create the state and dist arraylists to store the results for nts for current IGHJ in
                ArrayList stateArrayList = new ArrayList();
                ArrayList distArrayList = new ArrayList();
				//work through each nt in current sequence and create dists and emission probs
                for(int j = 1; j <= seq.length(); j++)
                {
                    //create a temp distribution for nt and save
					Distribution tempDist = DistributionFactory.DEFAULT.createDistribution(dna);
                    distArrayList.add(tempDist);
					//set the emission probs for the current nt
                    setGeneEmissionProb(seq, seq_string, j, tempDist, dna, 3, completeVGeneLength, 0, A_probability);
                    StateInfo stateInfo = new StateInfo(geneName, 'J', i, j, getChar(seq.symbolAt(j)), seq.length());
                    FakeAnnotation annotation = new FakeAnnotation(stateInfo);
                    //save emission probs
					stateArrayList.add(new SimpleEmissionState(geneName, annotation, advance, tempDist));
					
					//clean-up
					tempDist = null;
					stateInfo = null;
					annotation = null;
                }
				
				//add results for current sequence's nts to jStates
                jStates.add(stateArrayList);
                jStateDistributions.add(distArrayList);
				
				//clean-up
				stateArrayList = null;
				distArrayList = null;
				seq_string = null;
				geneName = null;
            }

        }
        catch(Exception e)
        {
            throw new Error("createJStates: Can't create distributions");
        }
    }


   private double getRelativePprob(double probabilities[], int offset)
    {
        double result = 1.0D;
        int reversed_offset = probabilities.length - offset;
        int i;
        for(i = 0; i < reversed_offset; i++)
        {
            result *= probabilities[i];
        }

        if(i < probabilities.length)
        {
            double notNextStateProb = 1.0D - probabilities[i];
            result *= notNextStateProb;
        }
        return result;
    }

    private double getCompleteRelativeProb(double probabilities[])
    {
        double sum = 0.0D;
        for(int i = 0; i < probabilities.length; i++)
        {
            sum += getRelativePprob(probabilities, i);
        }

        return sum;
    }

    //private void create_VD_N_part_1_Transitions_multi_P(SimpleMarkovModel vdj, EmissionState VD_N_53_states[], DotState X7a_states[], DotState X7b_states[], DotState X8, DotState X9, ArrayList P_D_end_states_array_list)
    private void create_VD_N_part_1_Transitions_multi(SimpleMarkovModel vdj, EmissionState VD_N_53_states[], DotState X7a_states[], DotState X7b_states[], DotState X8, DotState X9)
    {
        try
        {
            EmissionState currN = null;
            EmissionState nextN = null;
            //for(int i = 0; i < P_D_end_states_array_list.size(); i++)
            //{
            //    EmissionState P_D_end_state_array[] = (EmissionState[])P_D_end_states_array_list.get(i);
            //    vdj.createTransition(X7a_states[i], P_D_end_state_array[0]);
            //still need this part linking the X7a to the state that defines the start of the IGHJ
            //these need to be done per germline d gene -> move to being part of the create D
            //vdj.createTransition(X7a_states[i], X8);
            //vdj.createTransition(X7b_states[i], X8);
            //}

            vdj.createTransition(X8, X9);
            vdj.createTransition(X8, VD_N_53_states[0]);
            
            //the to magical states
            vdj.createTransition(X8, vdj.magicalState());
            vdj.createTransition(X9, vdj.magicalState());

            for(int i = 0; i < VD_N_53_states.length - 1; i++)
            {
                currN = VD_N_53_states[i];
                nextN = VD_N_53_states[i + 1];
                vdj.createTransition(currN, nextN);
                vdj.createTransition(currN, X9);

                //trans to magical states
                vdj.createTransition(currN, vdj.magicalState());
            }

            vdj.createTransition(nextN, X9);
            //trans to magical state for final
            vdj.createTransition(nextN, vdj.magicalState());
        }
        catch(Exception e)
        {
            throw new Error("create_VD_N_part_1_Transitions_multi():  " + e.getMessage());
        }
    }

    //private void create_VD_N_part_1_Transitions(SimpleMarkovModel vdj, EmissionState VD_N_53_states[], DotState X1a, DotState X1b, DotState X2, DotState X3, EmissionState P_V_end_states[])
    private void create_VD_N_part_1_Transitions(SimpleMarkovModel vdj, EmissionState VD_N_53_states[], DotState X1a, DotState X1b, DotState X2, DotState X3)
    {
        try
        {
            //create temp holders for neighbouring emission states
            EmissionState currN = null;
            EmissionState nextN = null;
            //not using the P states so don't need to create links to them
            //vdj.createTransition(X1a, P_V_end_states[0]);
            //link the v end states to the start of the N states
            vdj.createTransition(X1a, X2);
            vdj.createTransition(X1b, X2);
            //link the start of the N to the end of the N so that we can have rearrangements that lack N nts
            vdj.createTransition(X2, X3);
            //link the start of the N to the first N1 nucleotide
            vdj.createTransition(X2, VD_N_53_states[0]);
            
            //for the partial sequences will want to have a prob for ending the sequence
            vdj.createTransition(X1a, vdj.magicalState());
            vdj.createTransition(X1b, vdj.magicalState());
            vdj.createTransition(X2, vdj.magicalState());
            vdj.createTransition(X3, vdj.magicalState());
            //vdj.createTransition(VD_N_53_states[0], vdj.magicalState());  //this is done in the for loop below

            //for each of the N nucleotides create transition states between a position (i) and it downstream neighbour (i+1)
            for(int i = 0; i < VD_N_53_states.length - 1; i++)
            {
                //get the current N state
                currN = VD_N_53_states[i];
                //get the neighbouring N state
                nextN = VD_N_53_states[i + 1];
                //create transition between the current and it neighbour
                vdj.createTransition(currN, nextN);
                //create transition between the current and the end of the N state (so allow current to be the final N nt)
                vdj.createTransition(currN, X3);
                //create transition between the current and the magicalState
                vdj.createTransition(currN, vdj.magicalState());
            }
            //once all N states are exhausted create a link between the final N state and the end of the N region
            vdj.createTransition(nextN, X3);
            //final N to magicalState
            vdj.createTransition(nextN, vdj.magicalState());
        }
        catch(Exception e)
        {
            throw new Error("create_VD_N_part_1_Transitions():  " + e.getMessage());
        }
    }

    private void create_VD_N_part_2_Transitions(SimpleMarkovModel vdj, DotState X3, DotState X5a_states[], DotState X5b_states[])
    {
        try
        {
            //X5b states contains the dot state for start of ighd (with exo) for each germline possibility
            for(int i = 0; i < X5b_states.length; i++)
            {
                //create trans from end of N to start of the current germline ighd
                vdj.createTransition(X3, X5b_states[i]);
                //create trans from non-exo state for the current germline to the exo state (which will then transition to the dstates
                vdj.createTransition(X5b_states[i], X5a_states[i]);
                
                //create the trans to the magical states
                vdj.createTransition(X5b_states[i], vdj.magicalState());
                vdj.createTransition(X5a_states[i], vdj.magicalState());
            }

        }
        catch(Exception e)
        {
            throw new Error("create_VD_N_part_2_Transitions_multi_P(): Can't create transitions");
        }
    }



    //private void set_VD_N_part_1_Transitions_multi_P(SimpleMarkovModel vdj, EmissionState N_53_states[], double N_53_Probs[], DotState X7a_states[], DotState X7b_states[], DotState X8, DotState X9, ArrayList P_end_states_list, double P_probs[])
    private void set_VD_N_part_1_Transitions_multi(SimpleMarkovModel vdj, EmissionState N_53_states[], double N_53_Probs[], DotState X7a_states[], DotState X7b_states[], DotState X8, DotState X9)
    {
        try
        {
            EmissionState currN = null;
            EmissionState nextN = null;
            //double X7aiToP_endProb = P_probs[0];
            //double X7aiToX8Prob = 1.0D - X7aiToP_endProb;

            double X7biToX8Prob = 1.0D - toMagicalStateProb;
            Distribution dist;

            double X8ToN1Prob = N_53_Probs[0];
            //double X8ToN1Prob = N_53_Probs[0] * (1 - toMagicalStateProb);
            double X8ToX9Prob = 1.0D - X8ToN1Prob - toMagicalStateProb;
            dist = vdj.getWeights(X8);
            dist.setWeight(N_53_states[0], X8ToN1Prob);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);
            if (DEBUGGING) {
            	System.out.println("X8 to N2(1) set to " + X8ToN1Prob);
            }
            dist.setWeight(X9, X8ToX9Prob);
            if (DEBUGGING) {
            	System.out.println("X8 to X9 set to " + X8ToX9Prob);
            }
            for(int i = 0; i < N_53_states.length - 1; i++)
            {
                currN = N_53_states[i];
                nextN = N_53_states[i + 1];
                double NToNProb = N_53_Probs[i + 1];
                double NToX9Prob = 1.0D - NToNProb - toMagicalStateProb;
                dist = vdj.getWeights(currN);
                dist.setWeight(nextN, NToNProb);
                if (DEBUGGING) {
                	System.out.println("currN to nextN set to " + NToNProb);
                }
                dist.setWeight(X9, NToX9Prob);
                if (DEBUGGING) {
                	System.out.println("currN to X9 set to " + NToX9Prob);
               	}
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
            }

            double finalNtoX9 = 1.0D - toMagicalStateProb;
            dist = vdj.getWeights(nextN);
            dist.setWeight(X9, finalNtoX9);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);
            if (DEBUGGING) {
            	System.out.println("final N to X9 set to 1.0D");
            }
        }
        catch(Exception e)
        {
            throw new Error("set_VD_N_part_1_Transitions_multi():  " + e.getMessage());
        }
    }

    private void set_VD_N_part_1_Transitions(SimpleMarkovModel vdj, EmissionState N_53_states[], double N_53_Probs[], DotState X1a, DotState X1b, DotState X2, DotState X3)
    {
        try
        {
            // vars to hold the current and the next nt/state in the region
            EmissionState currN = null;
            EmissionState nextN = null;
			
            double X1aToX2Prob = 1.0D - toMagicalStateProb;
            //get weights for X1a state
            Distribution dist = vdj.getWeights(X1a);
            //set prob for moving from X1a to X2
            dist.setWeight(X2, X1aToX2Prob);
            if (DEBUGGING) {
            	System.out.println("X1a to X2 set to " +  X1aToX2Prob);
            }
			
            //prob for moving from exo v end to start of N region
            double X1bToX2Prob = 1.0D - toMagicalStateProb;
            //get weights for the X1b (exo rem gene end)
            dist = vdj.getWeights(X1b);
            //set weight for moving to start of N
            dist.setWeight(X2, X1bToX2Prob);
            if (DEBUGGING) {
            	System.out.println("X1b to X2 set to " +  X1bToX2Prob);
            }
			
            //get prob for moving to first N nt - from the N prob array
            double X2ToN1Prob = N_53_Probs[0];
            //get prob for moving straight to end of the N region - so prob of NOT moving to first N
            double X2ToX3Prob = 1.0D - X2ToN1Prob - toMagicalStateProb;
            //get weights for X2 (start N region)
            dist = vdj.getWeights(X2);
            //set the trans weights for moving to either first N or skipping the N region
            dist.setWeight(N_53_states[0], X2ToN1Prob);
            dist.setWeight(X3, X2ToX3Prob);
            if (DEBUGGING) {
            	System.out.println("X2 to N1(1) set to " +  X2ToN1Prob);
            	System.out.println("X2 to X3 set to " +  X2ToX3Prob);
            }
			
            //set the magical state trans prob in accordance with the alignment model selected
            dist = vdj.getWeights(X1a);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);
            dist = vdj.getWeights(X1b);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);
            dist = vdj.getWeights(N_53_states[0]);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);
            dist = vdj.getWeights(X2);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);
			
			
            //for each of the N states
            for(int i = 0; i < N_53_states.length - 1; i++)
            {
                //get the current state and it's downstream neighbour
                currN = N_53_states[i];
                nextN = N_53_states[i + 1];
				
                //get prob from prob array for adding (i)th N nucleotide
                double NToNProb = N_53_Probs[i + 1];
                //get prob of not adding another nt, hence prob adding no further N
                double NToX3Prob = 1.0D - NToNProb - toMagicalStateProb;
                //get weights for current N state/nt
                dist = vdj.getWeights(currN);
                //set weight of adding another N or moving to the end of the N
                dist.setWeight(nextN, NToNProb);
                dist.setWeight(X3, NToX3Prob);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                if (DEBUGGING) {
                	System.out.println("currN to nextN set to " +  NToNProb);
                	System.out.println("currN to X3 set to " +  NToX3Prob);
                }
            }
			
            //if get to the final N state, the set prob of moving to end of N state to 1
            double finalNtoX3 = 1.0D - toMagicalStateProb;
            dist = vdj.getWeights(nextN);
            dist.setWeight(X3, finalNtoX3);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);
            if (DEBUGGING) {
            	System.out.println("finalN to X3 set to 1.0D");
            }
        }
        catch(Exception e)
        {
            throw new Error("set_VD_N_part_1_Transitions():  " + e.getMessage());
        }
    }

    //altered to get rid of the P reference
    //private void set_VD_N_part_2_Transitions_multi_P(SimpleMarkovModel vdj, DotState X3, DotState X5b_states[])
    private void set_VD_N_part_2_Transitions_multi(SimpleMarkovModel vdj, DotState X3, DotState X5b_states[])
    {
        try
        {
            //holders for the current and neighbouring state
            EmissionState currN = null;
            EmissionState nextN = null;

            //get the number of possible germline IGHD genes
            double NUMBER_OF_DGENES = X5b_states.length;


            //prob of moving from X3 (state between N and gene)
            //equal probs given to each of the germline genes
            double X3ToX5bProb = ((1.0D - toMagicalStateProb) / (double)X5b_states.length);
            //get the weights for the X3 states
            Distribution dist = vdj.getWeights(X3);
            dist.setWeight(vdj.magicalState(), toMagicalStateProb);

            //for each of the germline genes
            for(int i = 0; i < X5b_states.length; i++)
            {
                //changed this from 1.0 to X3ToX5bProb
                dist.setWeight(X5b_states[i],X3ToX5bProb);
                if (DEBUGGING) {
                	System.out.println("X3 to X5b[" + i +"] set to " +  X3ToX5bProb);
                }
            }

        }
        catch(Exception e)
        {
            throw new Error("set_VD_N_part_2_Transitions_multi_P():  " + e.getMessage());
        }
    }


    private void createVTransitions(SimpleMarkovModel vdj, ArrayList vStates, DotState X1a, DotState X1b, ProbabilityHolder probHolder, int chainType)
    {
        try
        {
            //create some temp holders for neighbouring states
            EmissionState currVState = null;
            EmissionState nextVState = null;

            //create state for start of V-D-J rearrangement from the magical start
            //ie commence the markov model
            vdj.createTransition(vdj.magicalState(), (EmissionState)vStates.get(0));

            double V_end_exo_probs[] = null;
            if (chainType == 1) {// heavy chain
            	//get the first vstate = to first nt in the ighv
	            EmissionState vState = (EmissionState)vStates.get(0);
	            //get the annotation for the first nt (contains token and emission dist)
	            Annotation anno = vState.getAnnotation();
	            StateInfo si = (StateInfo)anno.getProperty(null);
	            //get the name of the IGHV gene as stored in teh state info - set in the alignment thread with finding the ighv match
	            //use this to determine the IGHV family number
	            int familyIndex = findFamilyIndex(si.geneName);
	            //get the mean exo probability for the IGHV's gene family
	            double mean = probHolder.V_end_exo_mean_Fields[familyIndex - 1];
	            //get the std dev for the IGHV's gene familiy
	            double stdDev = probHolder.V_end_exo_stdDev_Fields[familyIndex - 1];
	            //calculates the exo probabilities array using a normal distribution with the std dev and mean provided
	            V_end_exo_probs = probHolder.getExoProbArray(mean, stdDev, MIN_PROB_LIMIT, true);
            } else if (chainType ==2) { // kappa chain
            	V_end_exo_probs = probHolder.Kappa_V_end_exo_Fields;
            } else if (chainType == 3) { //lambda chain
            	V_end_exo_probs = probHolder.Lambda_V_end_exo_Fields;
            }

            //for each 'state' in the IGHV (ie each nt in the ighv gene) create transition between
            //it and it's neighbouring state (nt)
            for(int i = 0; i < vStates.size() - 1; i++)
            {
                currVState = (EmissionState)vStates.get(i);
                nextVState = (EmissionState)vStates.get(i + 1);
                vdj.createTransition(currVState, nextVState);
                if(vStates.size() - i <= V_end_exo_probs.length)
                {
                    //for nucleotides that are in teh region that can potentially contain exo removal
                    //based on the number of values in the V_end_exo_probs array then create transition
                    //to the exo removals state
                    vdj.createTransition(currVState, X1b);
                }
            }

            // if the vstate if the final nt of the V then can move to the no exo removal state
            vdj.createTransition(nextVState, X1a);

        }
        catch(Exception e)
        {
            throw new Error("createVTransitions(): " + e.getMessage());
        }
    }

    private int findFamilyIndex(String geneName)
    {
        //var to hold the return value -> return val is the family for the input gene name
        int result = 0;
        
        if(geneName.indexOf("DIR") != -1) {
            result = 5;
        }
        
        //first try to find IGH
        int startIGH = geneName.indexOf("IGH");
        int familyNumberIndex = 0;
        
        //failed to find IGH, so check if can use "-"
        if(startIGH == -1) {
            startIGH = geneName.indexOf("-");
        
            if (startIGH == -1) {
            	
            	startIGH = geneName.indexOf("JH");
            	
            	if (startIGH == -1) {
            		//failed to find -, IGH or JH so throw an error
           			throw new Error("Invalid gene sequence name, neither IGH nor DIR nor JH: " + geneName);
           		} else {
           			familyNumberIndex = startIGH + 2;
           		}
           	} else {
           		//set the location of the family number in the gene name string
           		familyNumberIndex = startIGH - 1;		
           	}
           	
        } else {
        	//set the location of the family number in the gene name string
            familyNumberIndex = startIGH + 4;
        }
        
        //get the char at the position in the gene name string
        char number = geneName.charAt(familyNumberIndex);
        //set the result to the integer parsed from the gene name string
        result = Integer.parseInt("" + number);
        
        //return
        return result;
    }

    private void setVTransitions(SimpleMarkovModel vdj, ArrayList vStates, DotState X1a, DotState X1b, ProbabilityHolder probHolder, int chainType)
    {
        try
        {
        	double V_end_exo_probs[] = null;
        	if (chainType == 1) {// heavy chain
	        	EmissionState vState = (EmissionState)vStates.get(0);
	            Annotation anno = vState.getAnnotation();
	            StateInfo si = (StateInfo)anno.getProperty(null);
	            int familyIndex = findFamilyIndex(si.geneName);
	            double mean = probHolder.V_end_exo_mean_Fields[familyIndex - 1];
	            double stdDev = probHolder.V_end_exo_stdDev_Fields[familyIndex - 1];
	            V_end_exo_probs = probHolder.getExoProbArray(mean, stdDev, MIN_PROB_LIMIT, false);
            } else if (chainType == 2) { //kappa chain
            	V_end_exo_probs = probHolder.Kappa_V_end_exo_Fields;
            } else if (chainType == 3) { //lambda chain
            	V_end_exo_probs = probHolder.Lambda_V_end_exo_Fields;
            }	
            
	        EmissionState currVState = null;
            EmissionState nextVState = null;
            double MagicalToVProb = 1.0D;
            Distribution dist = vdj.getWeights(vdj.magicalState());
            dist.setWeight((EmissionState)vStates.get(0), MagicalToVProb);
            if (DEBUGGING) {
            	System.out.println("magical to V1 set to " +  MagicalToVProb);
            }
            int v_size = vStates.size();
            int probArrayIndex = 1;
            for(int i = 0; i < v_size - 1; i++)
            {
                currVState = (EmissionState)vStates.get(i);
                nextVState = (EmissionState)vStates.get(i + 1);
                dist = vdj.getWeights(currVState);
                if(v_size - i <= V_end_exo_probs.length)
                {
                    double temp = V_end_exo_probs[v_size - i - 1];
                    double VToVProb = 1.0D - temp;
                    dist.setWeight(nextVState, VToVProb);
                    dist.setWeight(X1b, temp);
                    if (DEBUGGING) {
                    	System.out.println("currV to nextV set to " + VToVProb + " and currV to X1b set to " + temp);
                    }
                } else
                {
                    dist.setWeight(nextVState, 1.0D);
                    //System.out.println("currV to nextV set to 1.0");
                }
            }

            dist = vdj.getWeights(nextVState);
            dist.setWeight(X1a, 1.0D);
            if (DEBUGGING) {
            	System.out.println("final V to X1a set to 1.0");
            }
        }
        catch(Exception e)
        {
            throw new Error("SetVTransitions(): " + e.getMessage());
        }
    }

    //changed def to elim p states
    //private void createDTransitions_multi_P(SimpleMarkovModel vdj, DotState X5b_states[], DotState X6_states[], ArrayList dStates, DotState X7a_states[], DotState X7b_states[], ProbabilityHolder probHolder)
    private void createDTransitions_multi(SimpleMarkovModel vdj, DotState X5b_states[], DotState X5a_states[], ArrayList dStates, DotState X7a_states[], DotState X7b_states[], ProbabilityHolder probHolder, DotState X8)
    {
        try
        {
            //arrays to hold the exo nuclease probs for the start and end of the IGHD gene
            double D_start_exo_probs[] = (double[])null;
            double D_end_exo_probs[] = (double[])null;

            int dSize = 0;
            //for each germline gene in the dstates
            for(int i = 0; i < dStates.size(); i++)
            {
                //get the nucleotide states for the given germline gene
                //it is in the form of an array list with entry for each nt
                ArrayList dSequence = (ArrayList)dStates.get(i);

                //create the trans from end exo to start of the N2
                vdj.createTransition(X7a_states[i],X8);
                vdj.createTransition(X7b_states[i],X8);

                //get the length of the current germline gene which is equal to the number of states in the dSequence
                dSize = dSequence.size();
                //get the first nucleotide/state
                EmissionState dState = (EmissionState)dSequence.get(0);
                //get the annotation associated with the first nucleotide
                Annotation anno = dState.getAnnotation();
                StateInfo si = (StateInfo)anno.getProperty(null);
                //use the genename from the annotation to get the gene family number
                //gene family number is needed to get exo nuclease data as this is family based
                int familyIndex = findFamilyIndex(si.geneName);
                //get the mean exo for the given gene family for the 5` end (start) of the ighd
                double mean_start = probHolder.D_start_exo_mean_Fields[familyIndex - 1];
                //get the standard deviation for the given gene family for the 5` end (start) of the ighd
                double stdDev_start = probHolder.D_start_exo_stdDev_Fields[familyIndex - 1];
                //create the exo prob array for the 5` exo removals from the IGHD gene, based on normal dist using
                //the supplied mean and std dev
                D_start_exo_probs = probHolder.getExoProbArray(mean_start, stdDev_start, MIN_PROB_LIMIT, false);
                //get the mean for the 3` exo removals (end) of the ighd gene family
                double mean_end = probHolder.D_end_exo_mean_Fields[familyIndex - 1];
                //get teh std dev for the 3` end of the ighd gene family
                double stdDev_end = probHolder.D_end_exo_stdDev_Fields[familyIndex - 1];
                //create the array of probs for 3` (end) exo removal based on normal distribution using mean and std dev provided

                D_end_exo_probs = probHolder.getExoProbArray(mean_end, stdDev_end, MIN_PROB_LIMIT, true);


                //create for move from non-exo D start to D1
                vdj.createTransition(X5a_states[i], dState);
                vdj.createTransition(X5b_states[i], dState);

                //create D1 to magical state
                //this is covered by the for loop later
                //vdj.createTransition(dState, vdj.magicalState());

                //for each state/nt in the current ighd sequence create the transitions between each d state
                // and the 5` and 3` exonuclease states (this is only done for portions of the sequence that fall
                // within the bounds of the 5` and 3` exo removal - based on the lengths given in the prob arrays)
                for(int j = 1; j < dSize - 1; j++)
                {
                    //get the current dstate
                    dState = (EmissionState)dSequence.get(j);

                    // if it is within the region subject to 5` exo removal
                    if(j <= D_start_exo_probs.length)
                    {
                        //create the trans state
                        vdj.createTransition(X5b_states[i], dState);
                    }

                    //if it is within the region subject to 3` exo removal
                    if(dSize - j <= D_end_exo_probs.length)
                    {
                        //create the trans state
                        vdj.createTransition(dState, X7b_states[i]);
                    }
                }

            }

            //create a holder for the neighbouring downstream dstate/nt
            EmissionState nextD = null;

            //for each of the germline sequences
            for(int i = 0; i < dStates.size(); i++)
            {
                //get the states for the current sequence
                ArrayList dSequence = (ArrayList)dStates.get(i);
                //for each of the dStates/nts create the trans states between it and its neighbour
                for(int j = 0; j < dSequence.size() - 1; j++)
                {
                    //get the current state
                    EmissionState currD = (EmissionState)dSequence.get(j);
                    //get the neighbouring state
                    nextD = (EmissionState)dSequence.get(j + 1);
                    //create the transitsion
                    vdj.createTransition(currD, nextD);

                    //create the transition
                    vdj.createTransition(currD, vdj.magicalState());
                }

                //for the final nt of any germline sequence create transition to the non-exo D end state
                vdj.createTransition(nextD, X7a_states[i]);
                vdj.createTransition(X5b_states[i], nextD);

                //final to magical state, plus the 3' exo to magical state
                vdj.createTransition(nextD, vdj.magicalState());
                vdj.createTransition(X7a_states[i], vdj.magicalState());
                vdj.createTransition(X7b_states[i], vdj.magicalState());
            }

        }
        catch(Exception e)
        {
           throw new Error("createDTransitions_multi(): " + e.getMessage());
        }
    }

    private void setDTransitions_multi(SimpleMarkovModel vdj, DotState X5b_states[], DotState X5a_states[], ArrayList dStates, DotState X7a_states[], DotState X7b_states[], ProbabilityHolder probHolder, DotState X8)
    {
        try
        {
            //arrays to hold the start and end exonuclease probs
            double D_start_exo_probs[] = (double[])null;
            double D_end_exo_probs[] = (double[])null;
            //hold the current and neighbouring states
            EmissionState currD = null;
            EmissionState nextD = null;
            //hold the exo prob and the remainder
            double exo_prob = 0.0D;
            double remainder = 0.0D;
            //start of exo prob index
            int D_end_exo_prob_index = 1;
            //number of germline ighd genes
            int dGeneCount = dStates.size();

            //for each ighd from the germline repertoire
            for(int i = 0; i < dStates.size(); i++)
            {
                //get the nts/states for the current germline gene
                ArrayList dSequence = (ArrayList)dStates.get(i);
                //get the length of the current germline sequence which is equal to number of states
                int dSize = dSequence.size();
                //get the first state for the current germline gene
                EmissionState dState = (EmissionState)dSequence.get(0);
                //get the annotation associated with the first state
                Annotation anno = dState.getAnnotation();
                //use the annotation to determine the germline gene family number
                StateInfo si = (StateInfo)anno.getProperty(null);
                int familyIndex = findFamilyIndex(si.geneName);
                //get the mean exonuclease for the 5` end of the ighd gene for the associated ighd family
                double mean_start = probHolder.D_start_exo_mean_Fields[familyIndex - 1];
                //get the std dev for exonuclease the 51 end of the ighd gene family
                double stdDev_start = probHolder.D_start_exo_stdDev_Fields[familyIndex - 1];
                //create the exo probs array using normal dist and the provided mean and stdev for the current ighd family
                D_start_exo_probs = probHolder.getExoProbArray(mean_start, stdDev_start, MIN_PROB_LIMIT, false);
                //repeat for the 3` (end) end of the ighd gene
                double mean_end = probHolder.D_end_exo_mean_Fields[familyIndex - 1];
                double stdDev_end = probHolder.D_end_exo_stdDev_Fields[familyIndex - 1];
                D_end_exo_probs = probHolder.getExoProbArray(mean_end, stdDev_end, MIN_PROB_LIMIT, false);

                //get prob for zero exo rems from the start of the ighd gene
                double no_exo_start_prob = D_start_exo_probs[0];
                //prob for moving from X5b to the non-exo state (X5a) equal to prob of there being no exo removals
                double X5biToX5aiProb = no_exo_start_prob;
                Distribution dist = vdj.getWeights(X5b_states[i]);
                //set prob for moving to the non exo state equal to prob of zero exo rems
                dist.setWeight(X5a_states[i], X5biToX5aiProb);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                dist.setWeight(dState,0.0D);
                if (DEBUGGING) {
                	System.out.println("X5b[" + i + "] to X5a[" + i + "] set to " + X5biToX5aiProb);
                	System.out.println("X5b[" + i + "] to D1 set to 0.0");
                	System.out.println("X5b[" + i + "] to MagicalState set to " + toMagicalStateProb);
                }
                //no link to D1 from the non-exo state that i can see at the moment
                vdj.getWeights(X5a_states[i]);
                double X5atoD1 = 1.0D - toMagicalStateProb;
                dist.setWeight(dState, X5atoD1);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                if (DEBUGGING) {                
                	System.out.println("X5a[" + i + "] to D1 set to 1.0");
                }
                //sets the link between first and second D nts to 1, as we only start linking the dStates from the second nt below
                double D1toD2 = 1.0D - toMagicalStateProb;
                dist = vdj.getWeights(dState);
                dist.setWeight((EmissionState)dSequence.get(1), D1toD2);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                D_end_exo_prob_index = 1;

                //for each ighd state - ie each d nucleotide - starts from the second nt
                for(int j = 1; j < dSequence.size() - 1; j++)
                {
                    //get the current and the neighbouring state
                    currD = (EmissionState)dSequence.get(j);
                    nextD = (EmissionState)dSequence.get(j + 1);

                    remainder = 1.0D - toMagicalStateProb;
                    //if the current state is within possible 5` exo removals
                    if(j < D_start_exo_probs.length)
                    {
                        //get prob for associated level of exo removals
                        exo_prob = D_start_exo_probs[j];
                        //get the weights for the 5` d exo
                        dist = vdj.getWeights(X5b_states[i]);
                        //set the weight for moving from X5 to currD equal to exo prob
                        dist.setWeight(currD, exo_prob);
                        dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                        //am assuming that the exo probs get lower as they are probs of observing given levels of exo
                        //but as we get more nts into the ighd it should be more likely to move out of this state
                        //as more chance the sequence is part of the D rather than not
                        //see what happens for (1-exo_prob)
                        //double tempProb = 1.0D - exo_prob;
                        //dist.setWeight(currD, tempProb);
                        if (DEBUGGING) {
                        	System.out.println("X5b[" + i + "] to currD[" + j + "] set to " + exo_prob);
                        }
                    }

                    if(j==D_start_exo_probs.length) {
                        //if the 5' exo was completed - ie still in the exo state but up to the final dstate
                        double finalExoToD = 1.0D - toMagicalStateProb;
                        dist = vdj.getWeights(X5b_states[i]);
                        dist.setWeight(currD, finalExoToD);
                        dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                        if (DEBUGGING) {
                        	System.out.println("X5b (" + i + ") to currD(" + j + ") set to 1.0");
                        }
                    }

                    //if current state is within possible 3` exo removals
                    if(dSize - j <= D_end_exo_probs.length)
                    {
                        //get prob for associated level of exo rem
                        exo_prob = D_end_exo_probs[dSize - j - 1];
                        //get the prob of not being exo
                        //1 - prob exo
                        remainder -= exo_prob;
                        //get the weights for the currD
                        dist = vdj.getWeights(currD);
                        //set the weight of moving into the exo state
                        dist.setWeight(X7b_states[i], exo_prob);
                        if (DEBUGGING) {
                        	System.out.println("currD[" + j + "] to X7b[" + i + "] to  set to " + exo_prob);
                        }
                    }

                    //set the prob of moving to the next D state equal to the prob that there isn't exo removal
                    dist = vdj.getWeights(currD);
                    dist.setWeight(nextD, remainder);
                    dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                    if (DEBUGGING) {
                    	System.out.println("currD to nextD (" + j + ") set to " + remainder);
                    }
                }

                //for the final N state set the prob of moving to the no exo 3` state to 1
                double finalDtoNoExo = 1.0D - toMagicalStateProb;
                dist = vdj.getWeights(nextD);
                dist.setWeight(X7a_states[i], finalDtoNoExo);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                if (DEBUGGING) {
                	System.out.println("finalD to X7a (" + i + ") set to 1.0");
                }

                //for moving from the end states to the start of the N2
                double DendToX8 = 1.0D - toMagicalStateProb;
                dist = vdj.getWeights(X7a_states[i]);
                dist.setWeight(X8, DendToX8);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                dist = vdj.getWeights(X7b_states[i]);
                dist.setWeight(X8, DendToX8);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);

                if (DEBUGGING) {
                	System.out.println("X7a and X7b [" + i + " to X8 set to 1.0");
                }
            }

        }
        catch(Exception e)
        {
            throw new Error("setDTransitions_multi():  " + e.getMessage());
        }
    }

    private void createJTransitions(SimpleMarkovModel vdj, DotState X11a_states[], DotState X11b_states[], ArrayList jStates, ProbabilityHolder probHolder, int chainType, boolean removed_C_Region)
    {
        try
        {
            EmissionState currJ = null;
            EmissionState nextJ = null;
            double J_start_exo_probs[] = (double[])null;
            for(int i = 0; i < jStates.size(); i++)
            {
                ArrayList jSequence = (ArrayList)jStates.get(i);
                EmissionState jState = (EmissionState)jSequence.get(0);
                if (chainType == 1) {//heavy chain
                	Annotation anno = jState.getAnnotation();
	                StateInfo si = (StateInfo)anno.getProperty(null);
	                int familyIndex = findFamilyIndex(si.geneName);
	                double mean_start = probHolder.J_start_exo_mean_Fields[familyIndex - 1];
	                double stdDev_start = probHolder.J_start_exo_stdDev_Fields[familyIndex - 1];
	                J_start_exo_probs = probHolder.getExoProbArray(mean_start, stdDev_start, MIN_PROB_LIMIT, false);
                } else if (chainType == 2) { //kappa
                	J_start_exo_probs = probHolder.Kappa_J_start_exo_Fields;
                } else if (chainType == 3) { //lambda
                	J_start_exo_probs = probHolder.Lambda_J_start_exo_Fields;
                }
                
	            //need to create the non exo to J1 transition
                vdj.createTransition(X11a_states[i], jState);
                vdj.createTransition(X11b_states[i], jState);

                //create the trans to the magical for use when searching for partial seqs
                //vdj.createTransition(X11a_states[i], vdj.magicalState());
                //vdj.createTransition(X11b_states[i], vdj.magicalState());

                for(int j = 1; j <= J_start_exo_probs.length; j++)
                {
                    jState = (EmissionState)jSequence.get(j);
                    vdj.createTransition(X11b_states[i], jState);
                }


            }

            for(int i = 0; i < jStates.size(); i++)
            {
                ArrayList jSequence = (ArrayList)jStates.get(i);
                for(int j = 0; j < jSequence.size() - 1; j++)
                {
                    currJ = (EmissionState)jSequence.get(j);
                    nextJ = (EmissionState)jSequence.get(j + 1);
                    vdj.createTransition(currJ, nextJ);
                    //if(!removed_C_Region)
                    //if(!removed_C_Region)
                    //{
                        vdj.createTransition(currJ, vdj.magicalState());
                    //}
                }

                //if(removed_C_Region)
                //{
                //    vdj.createTransition(currJ, vdj.magicalState());
                //}
                vdj.createTransition(nextJ, vdj.magicalState());
            }

        }
        catch(Exception e)
        {
            throw new Error("createJTransitions(): " + e.getMessage());
        }
    }

    //altered to remove the P states
    //private void setJTransitions(SimpleMarkovModel vdj, DotState X11b_states[], DotState X11a_states[], DotState X12_states[], ArrayList jStates, ProbabilityHolder probHolder, boolean removed_C_Region)
    private void setJTransitions(SimpleMarkovModel vdj, DotState X11b_states[], DotState X11a_states[], ArrayList jStates, ProbabilityHolder probHolder, int chainType, boolean removed_C_Region)
    {
        try
        {
            //array for the j exo probs
            double J_start_exo_probs[] = (double[])null;
            //states to hold curr, neighbouring and final j
            EmissionState currJ = null;
            EmissionState nextJ = null;
            EmissionState lastJ = null;

            //number of germline j gene sequences
            int jGeneCount = jStates.size();

            //for each of the germline j gene sequences
            for(int i = 0; i < jStates.size(); i++)
            {
                //get the states/nts for the current germline sequence
                ArrayList jSequence = (ArrayList)jStates.get(i);
                //get the length of the current germline sequence
                int jSize = jSequence.size();
                //get the first ighj state/nt
                EmissionState jState = (EmissionState)jSequence.get(0);
                if (chainType == 1) {// heavy chain
	                //get the annotation associated with the first state/nt
	                Annotation anno = jState.getAnnotation();
	                //use the annotation to determine the ighj family number
	                StateInfo si = (StateInfo)anno.getProperty(null);
	                int familyIndex = findFamilyIndex(si.geneName);
	                //get the mean and std dev corresponding to the current ighj family
	                double mean_start = probHolder.J_start_exo_mean_Fields[familyIndex - 1];
	                double stdDev_start = probHolder.J_start_exo_stdDev_Fields[familyIndex - 1];
	                //use the mean and std dev to create exo probs using a normal distribution
	                J_start_exo_probs = probHolder.getExoProbArray(mean_start, stdDev_start, MIN_PROB_LIMIT, false);
                } else if (chainType ==2 ){ //kappa
                	J_start_exo_probs = probHolder.Kappa_J_start_exo_Fields;
                } else if (chainType == 3) { //lambda
                	J_start_exo_probs = probHolder.Lambda_J_start_exo_Fields;
                }
                //get the prob of no exo rems from the start of the ighj gene
                double J_start_no_exo_prob = J_start_exo_probs[0];
                //prob for moving from exo to non exo state for j start equal to prob that there are zero exo rems
                double X11biToX11aiProb = J_start_no_exo_prob - toMagicalStateProb;
                //set probs for moving from exo to non exo states
                Distribution dist = vdj.getWeights(X11b_states[i]);
                dist.setWeight(X11a_states[i], X11biToX11aiProb);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                dist.setWeight(jState,0.0D);
                if (DEBUGGING) {
                	System.out.println("X11b_states[" + i + "] to X11a_states [" + i + "] set to " + X11biToX11aiProb);
                	System.out.println("X11b_states[" + i + "] to J1 set to 0.0");
				}
                //need to set probs for X11a to J1
                double X11aToJ1 = 1.0D - toMagicalStateProb;
                dist = vdj.getWeights(X11a_states[i]);
                dist.setWeight(jState, X11aToJ1);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                if (DEBUGGING) {
                	System.out.println("X11a_states[" + i + "] to J1 set to 1.0");
                }
                //set prob for first to second j to one
                double J1toJ2 = 1.0D - toMagicalStateProb;
                dist = vdj.getWeights(jState);
                dist.setWeight((EmissionState)jSequence.get(1), J1toJ2);
                dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                if (DEBUGGING) {
                	System.out.println("J1 to J2 set to 1.0D");
                }

                //for each of the jstates/ighj nts up to the second last
                for(int j = 1; j < jSequence.size() - 2; j++)
                {
                    //get the current and the neighbouring nt
                    currJ = (EmissionState)jSequence.get(j);
                    nextJ = (EmissionState)jSequence.get(j + 1);
                    //also get the nt two upstream as need to be check for the end of the V-D-J rearrangement
                    lastJ = (EmissionState)jSequence.get(j + 2);
                    //if the current nt/state falls into the range of possible j start exo rems
                    if(j < J_start_exo_probs.length)
                    {
                        //get the prob for the associated with the current number of exo rems
                        double exo_prob = J_start_exo_probs[j];
                        //set the prob for moving from the exo rem state into the J gene equal to the exo prob
                        //check this...... should it be 1-prob of exo because prob of exo is for staying in exo state
                        //where as we are setting prob for moving from the exo state into the j state... CHECK!!!!
                        dist = vdj.getWeights(X11b_states[i]);
                        dist.setWeight(currJ, exo_prob);
                        dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                        if (DEBUGGING) {
                        	System.out.println("X11b_states[" + i + "] to currJ(" + j + ") set to " + exo_prob);
                        }
                    }

                    if(j == J_start_exo_probs.length) {
                      double JExoFinalToJ = 1.0D - toMagicalStateProb;
                      dist = vdj.getWeights(X11b_states[i]);
                      dist.setWeight(currJ, JExoFinalToJ);
                      dist.setWeight(vdj.magicalState(), toMagicalStateProb);
                      if (DEBUGGING) {
                      	System.out.println("X11b_states[" + i + "] to currJ(" + j + ") set to 1.0D");
                      }
                    }


                    //if the is no downstream sequence that has been removed then need to check if we have reached the end
                    double remainder = 1.0D;
					double JtoMagicalProb = toMagicalStateProb;
                    //if(removed_C_Region) //if downstream J nts have been removed by the JTraillingVectorFinder then need to allow jump to MagicalState prior to final nts of the J
                    //{
                        if (DEBUGGING) {
                        	System.out.println("for removed downstream J nts");                
                        }
                        JtoMagicalProb = 0.01D;  //set a prob of moving to magical state is downstream nts have been removed
                   // }
                    //set the prob of moving onto the next jstate equal to the prob that we are not at the end
					//use magical state prob to determine prob of not moving to magical state, ie. moving to next J nt
                    remainder = 1.0D - JtoMagicalProb;
                    dist = vdj.getWeights(currJ);
                    dist.setWeight(nextJ, remainder);  //set weight for moving to next J
                    dist.setWeight(vdj.magicalState(), JtoMagicalProb);  //set weight for moving the magical state if end of seq has been reached
                    if (DEBUGGING) {
                    	System.out.println("currJ(" + j + ") to nextJ set to " + remainder);
                    	System.out.println("currJ(" + j + ") to magical state " + JtoMagicalProb);
                    }
                }


                //yet to deal with the second last and final nucleotide
                double secondLastJtoMagicalProb;
                double secondLastJtoLastJProb;
                //if the c region has been removed
                if(removed_C_Region)
                {
                    secondLastJtoMagicalProb = 0.5D;
                    secondLastJtoLastJProb = 1.0D - secondLastJtoMagicalProb;
                } else //if the C region hasn't been removed
                {
                    secondLastJtoMagicalProb = 0.02D;
                    secondLastJtoLastJProb = 1.0D - secondLastJtoMagicalProb;
                }
                //set weights for the neighbouring nt moving into the magical state
                // or for moving onto the last J
                dist = vdj.getWeights(nextJ);
                dist.setWeight(vdj.magicalState(), secondLastJtoMagicalProb);
                dist.setWeight(lastJ, secondLastJtoLastJProb);
                if (DEBUGGING) {
                	System.out.println("second last J to magical set to " + secondLastJtoMagicalProb);
                	System.out.println("second last J to last J set to " + secondLastJtoLastJProb);
                }
                double JXlastToMagicalProb = 1.0D;
                //set weights for the last J moving to the magical state
                dist = vdj.getWeights(lastJ);
                dist.setWeight(vdj.magicalState(), JXlastToMagicalProb);
                if (DEBUGGING) {
                	System.out.println("last J to magical set to " +  JXlastToMagicalProb);
                }
            }

        }
        catch(Exception e)
        {
            throw new Error("setJTransitions:  " + e.getMessage());
        }
    }


    private double getSum(double array[])
    {
        double sum = 0.0D;
        for(int i = 0; i < array.length; i++)
        {
            sum += array[i];
        }

        return sum;
    }

    private int getSum(int array[])
    {
        int sum = 0;
        for(int i = 0; i < array.length; i++)
        {
            sum += array[i];
        }

        return sum;
    }

     private void setGeneEmissionProb(RichSequence seq, String seq_string, int nucl_pos, Distribution dist, FiniteAlphabet dna, int gene_type, int completeVGeneLength,
            int VGene_start_offset, double A_probability)
    {
        double exp_mutation_prob;
        switch(gene_type)
        {
        case 1: // '\001'
            exp_mutation_prob = ExponentialDecay.exponentialDecayVGene(seq_string, nucl_pos, VGene_start_offset);
            break;

        case 2: // '\002'
            // becuase of antigen selection rules, the mutability score of the DGene (actually the DGene and its VD and DJ junctions as well, but that's irrelevant in this matter)
            // has to be multiplied by 1.5
            exp_mutation_prob = 1.5D * ExponentialDecay.exponentialDecayDGene(seq_string, nucl_pos, completeVGeneLength);
            break;

        case 3: // '\003'
            exp_mutation_prob = ExponentialDecay.exponentialDecayJGene(seq_string, nucl_pos, completeVGeneLength);
            break;

        default:
            throw new Error("Gene type not V,D or J");
        }
        //get the 5mer with the current nt at the central position
        String pentaNucleotide = getPentaNucleotide(seq_string, nucl_pos);

        //use the 5mer to determine with the tri-nt or hotspot based mutability scores
        double mutability_score = G_mutability_score.pentaNucleotideScore(pentaNucleotide);
        
        if (DEBUGGING) {
        	System.out.println("mutability_score: " + mutability_score + ", A_prob: " + A_probability + ", exp: " + exp_mutation_prob);
        }
        //calculate the probability of mutation occuring at the current position
        //based on the A_probability calculated from the V-gene, the exponential decay based on the position within the sequence and the
        //mutability score
        double probability_of_mutation = A_probability * exp_mutation_prob * mutability_score;

        if (DEBUGGING) {
        	System.out.println("prob_of_mutation: " + probability_of_mutation);
        }
        //output the details to the gene_nucleotide_mutations_probabilities.txt file
        G_fostream.println("gene name: " + seq.getName() + "  nucl. position: " + nucl_pos + "  Probability of mutation  = " + probability_of_mutation);
        
        //set the probability of observing no mutation at the current position
        //this will be the emission prob for the germline nt
        double noMutationProb = 1.0D - probability_of_mutation;

        //get the tri-nucleotide that the current nt is part of --> used for calculation mutated emission probs
        String trinucleotide = getTriNucleotide(seq_string, nucl_pos);

        //get the nucleotide at the current position
        char mutateFrom = seq_string.charAt(nucl_pos - 1);

        //determine 6% and 2% of the total porbability of mutation for the current position
        double six_percent_of_mutation_prob = probability_of_mutation * 0.06D;   //6%
        double two_percent_of_probability = probability_of_mutation * 0.02D;     //2%

        //the total probability of mutation for a position is reduced by 6% to allow 2% to be added
        //to each of the mutated nt emission states to prevent probabilities of zero from occuring
        double reduced_probability_of_mutation = probability_of_mutation - six_percent_of_mutation_prob;
        
        //look at each possible mutation direction for the current nucleotide
        for(int i = 0; i < G_mutationSpectrum.UNIQUE_NUCLEOTIDE_ARRAY.length; i++)
        {
            //get the mutation direction
            char nucl = G_mutationSpectrum.UNIQUE_NUCLEOTIDE_ARRAY[i];
            
            //if the mutation direction is the same as the current nt then this is the same as no mutation occuring and the germline
            //nucleotide being present at the position
            if(nucl == mutateFrom)
            {
                //set the emission prob to the prob of no mutation occuring
                setNucleotideGeneEmission(nucl, noMutationProb, dist);
            } else //emission probs for mutated nts
            {
                //use the mutation_spectrum.txt file to determine the proportion of mutations for the current trint that
                //move in the direction of the current mutated nt
                double currMutationFraction = G_mutationSpectrum.getTNProbability(trinucleotide, nucl);
                G_mutationSpectrum.getClass();
                
                //if there is no entry for the tri-nt and mutation direction, print an error
                if(currMutationFraction == -1D)
                {
                    throw new Error("setGeneEmissionProb(): trinucleotide probability match lookup not found: " + trinucleotide + " : " + nucl);
                }
                
                if (DEBUGGING) {
                	System.out.println("reduced mut prob:" + reduced_probability_of_mutation + " curr mut frac:" + currMutationFraction);
               	}
                double relativeMutationProb = reduced_probability_of_mutation * currMutationFraction + two_percent_of_probability;
                
                //check that the relative probability isn't equal to 0 - throw an error if it is
                //it shouldn't be owing to the addition of the two_percent_of_probability
                if(relativeMutationProb == 0.0D)
                {
                    throw new Error("ZERO probability: " + relativeMutationProb);
                }

                //print a different error message when the probability for mutation direction in the given tri-nt
                //is zero ... print some ***** to highlight zero probs
                if(relativeMutationProb == 0.0D)
                {
                    G_fostream.println("************* mutation to nucl : " + nucl + " = " + relativeMutationProb);
                } else
                {
                    G_fostream.println("mutation to nucl : " + nucl + " = " + relativeMutationProb);
                }

                //set the emission probs for the mutation direction
                setNucleotideGeneEmission(nucl, relativeMutationProb, dist);
            }
        }

    }

    private String getTriNucleotide(String seq_string, int nucl_pos)
    {
        String trinucleotide = null;
        int string_nucl_pos = nucl_pos - 1;
        if(string_nucl_pos == 0)
        {
            G_mutationSpectrum.getClass();
            trinucleotide = 'n' + seq_string.substring(string_nucl_pos, string_nucl_pos + 2);
        } else
        if(string_nucl_pos == seq_string.length() - 1)
        {
            G_mutationSpectrum.getClass();
            trinucleotide = seq_string.substring(string_nucl_pos - 1, string_nucl_pos + 1) + 'n';
        } else
        {
            trinucleotide = seq_string.substring(string_nucl_pos - 1, string_nucl_pos + 2);
        }
        return trinucleotide;
    }

    private String getPentaNucleotide(String seq_string, int nucl_pos)
    {
        int string_nucl_pos = nucl_pos - 1;
        int max_string_pos = seq_string.length() - 1;
        String pentaNucleotide = null;
        if(string_nucl_pos == 0)
        {
            pentaNucleotide = "uu" + seq_string.substring(0, 3);
        } else
        if(string_nucl_pos == 1)
        {
            pentaNucleotide = "u" + seq_string.substring(0, 4);
        } else
        if(string_nucl_pos == max_string_pos)
        {
            pentaNucleotide = seq_string.substring(max_string_pos - 2, max_string_pos + 1) + 'u' + 'u';
        } else
        if(string_nucl_pos == max_string_pos - 1)
        {
            pentaNucleotide = seq_string.substring(max_string_pos - 3, max_string_pos + 1) + 'u';
        } else
        {
            pentaNucleotide = seq_string.substring(string_nucl_pos - 2, string_nucl_pos + 3);
        }
        return pentaNucleotide;
    }

    private void setNucleotideGeneEmission(char nucl, double probability, Distribution dist)
    {
        try
        {
            //to avoid errors caused by negative probabilities
            if (probability < 0) {
              //prob of less than zero (ie negative) is replaced by very small number
              //makes event very improbable without setting the prob to 0
              probability = Double.MIN_VALUE;
              //probability = 0;
            }
            switch(nucl)
            {
            case 97: // 'a'
                dist.setWeight(DNATools.a(), probability);
                break;

            case 99: // 'c'
                dist.setWeight(DNATools.c(), probability);
                break;

            case 103: // 'g'
                dist.setWeight(DNATools.g(), probability);
                break;

            case 116: // 't'
                dist.setWeight(DNATools.t(), probability);
                break;
            case 110: // 'n'
                dist.setWeight(DNATools.n(), 0.0D);
                break;
            }
        }
        catch(Exception ille)
        {
            throw new Error("setNucleotideGeneEmission():  " + ille.getMessage());
        }
    }

}
