package iHMMuneAlign;

import com.google.common.collect.Lists;
import org.biojava.bio.dist.AbstractDistribution;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dp.*;
import org.biojava.bio.symbol.*;
import org.biojava.utils.ChangeType;
import org.biojava.utils.ChangeVetoException;

import java.io.OutputStream;
import java.io.PrintStream;
import java.util.*;

/**
 * Created by fernie on 11/08/14.
 */
public class ModelWriter {

    private final MagicalState initState;
    private final MarkovModel model;
    private final PrintStream out;
    private final List<Symbol> trackSymbols;

    public ModelWriter(MarkovModel model, OutputStream out) {
        this.model = model;
        this.out = new PrintStream(out);

        this.initState = model.magicalState();
        trackSymbols = new ArrayList<Symbol>(4);
        for (Iterator<Symbol> it = ((FiniteAlphabet) model.emissionAlphabet()).iterator(); it.hasNext();) {
            trackSymbols.add(it.next());
        }
    }

    public void write() {
        try {
            writePreamble();
            writeLine();
            writeTrackSymbols();
            writeLine();
            writeStates();
            writeLine();
            writeMagicalState();
            writeLine();
        } catch (IllegalSymbolException e) {
            e.printStackTrace();
        }
    }

    private void writeMagicalState() throws IllegalSymbolException {
        out.println("STATE: ");
        out.println("\tNAME: INIT");
        out.println("\tPATH_NAME: S");
        writeTransitions(initState);
        writeLine();

        out.println("STATE: ");
        out.println("\tNAME: END");
        out.println("\tPATH_NAME: F");
        writeTransitions(initState);
        writeEmissions(initState.getDistribution());
        writeLine();
    }

    private void writeLine() {
        out.println();
    }

    private void writeStates() throws IllegalSymbolException {
        out.println("<STATE DEFINITIONS>");

        FiniteAlphabet states = model.stateAlphabet();
        for (Iterator<Symbol> si = states.iterator(); si.hasNext();) {
            State state = (State) si.next();

            // We'll deal with start/end separately
            if (state == initState) {
                continue;
            }

            out.println("STATE:");
            out.println("\tNAME: " + nameOf(state));
            out.println("\tPATH_LABEL: " + nameOf(state).charAt(0));

            writeTransitions(state);

            if (state instanceof SimpleEmissionState) {
                Distribution distribution = ((SimpleEmissionState) state).getDistribution();
                writeEmissions(distribution);
            } else if (state instanceof SimpleDotState) {
                writeEmissions(new AbstractDistribution() {
                    @Override
                    protected void setWeightImpl(AtomicSymbol sym, double weight) throws IllegalSymbolException, ChangeVetoException {

                    }

                    @Override
                    protected void setNullModelImpl(Distribution nullModel) throws IllegalAlphabetException, ChangeVetoException {

                    }

                    @Override
                    protected double getWeightImpl(AtomicSymbol sym) throws IllegalSymbolException {
                        return 0;
                    }

                    @Override
                    public Alphabet getAlphabet() {
                        return null;
                    }

                    @Override
                    public Distribution getNullModel() {
                        return null;
                    }
                });
            }
            writeLine();
        }
    }

    private void writeEmissions(Distribution distribution) throws IllegalSymbolException {
        out.println("EMISSION:\tTRACK: P(X)");
        out.println("\tORDER: 0");
        for (Symbol s : trackSymbols) {
            double weight = distribution.getWeight(s);
            out.print(weight);
            out.print(" ");
        }
        out.println();
    }

    private void writeTransitions(State state) throws IllegalSymbolException {
        out.println("TRANSITION:\tSTANDARD: P(X)");

        Distribution weights = model.getWeights(state);
        for (Iterator<Symbol> sj = ((FiniteAlphabet)
                weights.getAlphabet()).iterator(); sj.hasNext();) {
            Symbol toState = sj.next();

            double weight = weights.getWeight(toState);
            if (Double.isNaN(weight)) {
                weight = 0;
            }

            out.println("\t" + nameOf(toState) + ": " + weight);
        }
    }

    private String nameOf(Symbol toState) {
        return toState == initState ? "END" : toState.getName().replace('/', '_').replace('*', '_').replace(' ', '_') + "@" + System.identityHashCode(toState);
    }

    public void writeTrackSymbols() {
        out.println("<TRACK SYMBOL DEFINITIONS>");
        out.print("TRACK: ");

        for (Iterator<Symbol> it = ((FiniteAlphabet) model.emissionAlphabet()).iterator(); it.hasNext();) {
            Symbol s = it.next();
            out.print(nameOf(s).charAt(0));

            if (it.hasNext()) {
                out.print(", ");
            }
        }
        out.println();
    }

    private void writePreamble() {
        out.println("#STOCHHMM MODEL FILE");
    }
}
