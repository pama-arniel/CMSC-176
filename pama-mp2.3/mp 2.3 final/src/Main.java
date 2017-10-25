
public class Main {
    private static final String infile = "irisflowers.csv";
    private static final double lambda = 0.001;
    private static final int degree = 1;
    private static final int thetavalset = 1;
    private static final double alpha = 1;
    private static final int iters = 200;
    private static final int restarts = 100;

    //reads training samples from input file and performs operations after loading
    public static void main(String[] args) {
        Solution s = new Solution(lambda, degree, thetavalset, alpha, iters, restarts);
        s.load_data(infile);
    }
}
