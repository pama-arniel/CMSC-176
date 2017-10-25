import java.util.ArrayList;

public class Main {
    private static final String infile = "HousePricingRelationship.in";
    //private static final String infile = "sample";
    private static final double alpha = .00000001;
    private static final int iters = 100;

    //reads training samples from input file and do their corresponding operations
    public static void main(String[] args) {
        Solution s = new Solution();
        s.load_data(infile);

        ArrayList<Double> costhistory = s.gradientDescent(s.getX(), s.getY(), alpha , iters);
        s.graph(costhistory);
    }
}
