
public class Main {
    private static final String infile = "irisflowers.csv";

    //reads training samples from input file and do their corresponding operations
    public static void main(String[] args) {
        Solution s = new Solution();
        s.load_data(infile);
        s.scalefeatures(s.getX());
    }
}
