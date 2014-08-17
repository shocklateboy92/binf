package iHMMuneAlign;

import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.DP;
import org.biojava.bio.dp.MarkovModel;
import org.biojava.bio.dp.SimpleMarkovModel;
import org.biojava.bio.dp.State;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import java.io.IOException;
import java.io.ObjectOutputStream;
import java.util.Iterator;

/**
 * Created by fernie on 9/08/14.
 */
public class Utils {
    public static void printModel(MarkovModel markov_model) {
        SimpleMarkovModel model = (SimpleMarkovModel) markov_model;
        FiniteAlphabet stateAlphabet = model.stateAlphabet();

        System.out.println("======= BEGIN DEBUG OUTPUT ======");

        for (Iterator<Symbol> it = stateAlphabet.iterator(); it.hasNext();) {
            State fromState = (State) it.next();
            System.out.println(fromState.getName() + ":" + fromState.getClass().getName());
            System.out.println(fromState);
            try {
                Distribution transitions = model.getWeights(fromState);
                for (Iterator<Symbol> st = model.transitionsFrom(fromState).iterator(); st.hasNext();) {
                    State to = (State) st.next();
                    System.out.println(fromState.getName() + " -> " + to.getName() + " = " + transitions.getWeight(to));
                }
            } catch (IllegalSymbolException e) {
                e.printStackTrace();
            }
        }
        System.out.println(DNATools.a().getName());
        System.out.println("======== END DEBUG OUTPUT =======");
    }

    public static void printDP(DP dp) {
        try {
            for (State s : dp.stateList(dp.getModel())) {
                System.out.println(s);
            }
            System.out.println(dp.getForwardTransitions());
        } catch (BioException e) {
            e.printStackTrace();
        }
    }
}
