package iHMMuneAlign;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Map;

/**
 * Created by fernie on 29/03/15.
 */
public class Utils {
    static Map<String, PrintStream> files = new HashMap<String, PrintStream>();

    static void writeTo(String fileName, String data) {
        getWriter(fileName).println(data);
    }

    static PrintStream getWriter(String fileName) {
        if (!files.containsKey(fileName)) {
            try {
                files.put(fileName, new PrintStream(new File(fileName)));
            } catch (FileNotFoundException e) {
                e.printStackTrace();
            }
        }

        return files.get(fileName);
    }
}
