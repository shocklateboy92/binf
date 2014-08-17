package iHMMuneAlign;

import org.biojava.bio.BioException;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.*;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Symbol;

import java.io.FileOutputStream;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;

/**
 * Created by fernie on 17/08/14.
 */
public class DPWriter {
    DP dp;
    private List<Symbol> observations;
    MarkovModel model;
    PrintStream out;
    List<Symbol> trackSymbols;

    public DPWriter(DP dp, List<Symbol> observations) {
        this.dp = dp;
        this.observations = observations;
        this.model = dp.getModel();
//        this.out = new PrintStream(System.out);

        trackSymbols = new ArrayList<Symbol>(4);
        for (Iterator<Symbol> it = ((FiniteAlphabet) model.emissionAlphabet()).iterator(); it.hasNext();) {
            trackSymbols.add(it.next());
        }
    }

    public void write(OutputStream outputStream) throws BioException {
        this.out = new PrintStream(outputStream);
        State[] states = dp.stateList(model);

        out.println(states.length);
        out.println(trackSymbols.size());
        out.println(observations.size());

        out.println();

        writeTransitions(states);

        out.println();

        writeEmissions(states);

        out.println();

        writeObservations();
    }

    private void writeObservations() {
        for (Symbol s : observations) {
            int symbolIndex = trackSymbols.indexOf(s);
            out.print(symbolIndex);
            out.print(' ');
        }

        out.println();
    }

    private void writeEmissions(State[] states) throws IllegalSymbolException {
        for (State s: states) {
            if (s instanceof EmissionState) {
                Distribution distribution = ((EmissionState) s).getDistribution();
                for (Symbol sym : trackSymbols) {
                    double weight = distribution.getWeight(sym);
                    out.print(weight);
                    out.print(" ");
                }
            } else {
                for (Symbol sym : trackSymbols) {
                    out.print(0d);
                    out.print(' ');
                }
            }
            out.println();
        }
    }

    private void writeTransitions(State[] states) throws IllegalSymbolException {
        for (State i : states) {
            Distribution weights = model.getWeights(i);
            for (State j : states) {
                if (model.containsTransition(i, j)) {
                    double weight = weights.getWeight(j);
                    if (Double.isNaN(weight)) {
                        weight = 0;
                    }
                    out.print(weight + " ");
                } else {
                    out.print(0d + " ");
                }
            }
            out.println();
        }
    }
}
