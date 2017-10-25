import Jama.Matrix;

public class Main {
    private static final String infile = "irisflowers.csv";
    private static final double lambda = 0.001;
    private static final int degree = 1;
    private static final int thetavalset = 1;

    //reads training samples from input file and do their corresponding operations
    public static void main(String[] args) {
        Solution s = new Solution();
        s.load_data(infile);
        Matrix scaledX = s.scalefeatures(s.getX());
        s.engineerPolynomials(scaledX, degree);
        s.generateTheta(thetavalset);
        s.regularizedCost(lambda);
        System.out.println("CONSTANT VALUES USED:");
        System.out.println("Lambda: " + lambda);
        System.out.println("Degree: " + degree);
        System.out.println("Thetas: " + thetavalset);
    }
}
