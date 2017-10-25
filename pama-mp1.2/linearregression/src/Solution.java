import Jama.Matrix;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.data.xy.XYSeries;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.ui.ApplicationFrame;
import org.jfree.ui.RefineryUtilities;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;

/**
 * Created by Jose Arniel Pama on 8/26/2017.
 * Link to Java Matrices Library, JAMA: http://math.nist.gov/javanumerics/jama/
 * How to graph XY Plot: https://www.tutorialspoint.com/jfreechart/jfreechart_line_chart.htm
 */

public class Solution {
    private Matrix X; //the whole input data
    private Matrix Y; //the whole output data
    private Matrix theta; //theta values representing one hypothesis
    private int M; //number of training examples; row dimension
    private int N; //number of Xs in one training example; col dimension


    public Solution() {
        M = 0;
    }


    public void load_data(String filename) {
        BufferedReader br = null;

        try {
            String currline;
            br = new BufferedReader(new FileReader(filename));

            ArrayList<ArrayList<Double>> X = new ArrayList<>();
            ArrayList<Double> Y = new ArrayList<>();

            while ((currline = br.readLine()) != null) {
                ArrayList currentsample = new ArrayList<>(Arrays.asList(currline.split(",")));
                ArrayList<Double> tobeaddedtoX = new ArrayList<>();

                //parse as doubles, since a JAMA Matrix object is represented as double
                double n;
                int len = currentsample.size();
                for (int i = 0; i < (len - 1); i++) {
                    n = Double.parseDouble((String) currentsample.get(i));
                    tobeaddedtoX.add(n);
                }

                n = Double.parseDouble((String) currentsample.get(len - 1));
                Y.add(n); //the last item, index len-1, is the Y of every training example

                tobeaddedtoX.add(0, 1.0); //for the X0, which is 1; necessary for matrix operation later
                X.add(tobeaddedtoX);
            }

            M = X.size();
            double[][] x2D = alalTo2DArray(X);
            this.X = Matrix.constructWithCopy(x2D);
            N = this.X.getColumnDimension();
            this.theta = Matrix.constructWithCopy(generateTheta(x2D));
            this.Y = Matrix.constructWithCopy(YTo2DArray(Y));

            System.out.println("NUMBER OF TRAINING EXAMPLES: " + this.M);
            System.out.print("THE X MATRIX:");
            this.X.print(this.X.getColumnDimension(), 2);
            System.out.print("THE Y MATRIX:");
            this.Y.print(this.Y.getColumnDimension(), 2);
            System.out.print("INITIAL theta MATRIX: ");
            this.theta.print(this.theta.getColumnDimension(), 2);

        } catch (IOException e) {
            e.printStackTrace();
        } finally {
            try {
                if (br != null) br.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    public Matrix getX() {
        return this.X;
    }

    public Matrix getY() {
        return this.Y;
    }

    public int getN() {
        return this.N;
    }


    /////////////// START OF GRADIENT DESCENT FUNCTIONS  ///////////////
    public double[][] generateTheta(double[] given) {
        int xrows = given.length;
        double[][] res = new double[xrows][];
        for (int i = 0; i < xrows; i++) {
            double[] newrow = new double[1]; //Matrix theta has 1 column always

            //generate random number from 0 to M
            newrow[0] = given[i];
            res[i] = newrow;
        }

        return res;
    }

    //return the derivative of cost J given Matrix X, Y, and theta
    public double derivativeJ(Matrix X, Matrix Y, Matrix theta, int j){
        Matrix H = solveforH(X, theta);
        Matrix hy = solveforHMinusY(H, Y);
        Matrix Xj = X.getMatrix(0, X.getRowDimension() - 1, j, j);
        Matrix derJmatrix = Xj.transpose().times(hy);

        double res = derJmatrix.get(0, 0); //because derJmatrix is scalar
        return (res / M);
    }

    //return the cost history
    public ArrayList<Double> gradientDescent(Matrix X, Matrix Y, double alpha, int iters) {
        ArrayList<Double> costhistory = new ArrayList<>();
        int N = this.getN(); //number of X or number of cols
        double[] thetaArr = new double[N]; //there are N number of theta values

        Matrix theta;
        System.out.println("PROGRESSION OF THETA VALUES OVER TIME:");
        for (int i = 1; i <= iters; i++) {
            theta = Matrix.constructWithCopy(generateTheta(thetaArr));
            costhistory.add( cost(this.X, this.Y, theta) );
            for (int j = 0; j < N; j++) {
                thetaArr[j] = thetaArr[j] - (alpha * derivativeJ(X, Y, theta, j));
                //System.out.print( String.format("%.2f", thetaArr[j]) + " ");
                System.out.print( thetaArr[j] + " ");
            }

            System.out.println();
        }

        System.out.println();
        return costhistory;
    }

    public void graph(ArrayList<Double> costhistory){
        new Graph(costhistory);
    }

    public class Graph extends ApplicationFrame {
        public Graph(ArrayList<Double> costhistory) {
            super("Gradient Descent");

            final XYSeries series = new XYSeries("Cost");
            int iters = costhistory.size();
            for(int i = 0; i < iters; i++){
                series.add(i, costhistory.get(i));
            }
            System.out.println("LAST CALCULATED COST: " + costhistory.get(iters-1));

            final XYSeriesCollection data = new XYSeriesCollection(series);
            final JFreeChart chart = ChartFactory.createXYLineChart(
                    "Cost over iteration",
                    "Iteration",
                    "Cost",
                    data,
                    PlotOrientation.VERTICAL,
                    true,
                    true,
                    false
            );

            final ChartPanel chartPanel = new ChartPanel(chart);
            chartPanel.setPreferredSize(new java.awt.Dimension(700, 400));
            setContentPane(chartPanel);
            pack();
            RefineryUtilities.centerFrameOnScreen(this);
            setVisible(true);
        }

    }
    /////////////// END OF GRADIENT DESCENT FUNCTIONS  ///////////////

    /////////////// START OF HELPER FUNCTIONS  ///////////////
    public void print2D(ArrayList<ArrayList> alal) {
        for (ArrayList al : alal) {
            for (Object num : al) {
                System.out.print(num.toString() + " ");
            }
            System.out.println();
        }
    }

    public void print2D(double[][] alal) {
        for (int i = 0; i < alal.length; i++) {
            double[] curr = alal[i];
            for (int j = 0; j < curr.length; j++) {
                System.out.print(curr[j] + " ");
            }
            System.out.println();
        }
    }

    //transforms a 2D arraylist to a 2D array
    public double[][] alalTo2DArray(ArrayList<ArrayList<Double>> alal) {
        int alallen = alal.size();
        double[][] res = new double[alallen][];
        for (int i = 0; i < alallen; i++) {
            ArrayList<Double> row = alal.get(i);

            //perform equivalent 'toArray' operation
            int rowlen = row.size();
            double[] copy = new double[rowlen];

            //manually loop and set individually; simply putting "Integer[] copy = (Integer[]) row.toArray()" does not work
            for (int j = 0; j < rowlen; j++) {
                copy[j] = row.get(j);
            }

            res[i] = copy;
        }

        return res;
    }

    //transforms Y to a 2D array
    public double[][] YTo2DArray(ArrayList<Double> Y) {
        int ylen = Y.size();
        double[][] res = new double[ylen][];
        for (int i = 0; i < ylen; i++) {
            double[] newrow = new double[1]; //Matrix Y has 1 column always
            newrow[0] = Y.get(i);
            res[i] = newrow;
        }

        return res;
    }

    //this is supposed to generate the theta values randomly, but sets each to 2 at the moment
    public double[][] generateTheta(double[][] X) {
        int xrows = X[0].length;
        double[][] res = new double[xrows][];
        for (int i = 0; i < xrows; i++) {
            double[] newrow = new double[1]; //Matrix theta has 1 column always
            newrow[0] = 0; //set to zero lang sa
            res[i] = newrow;
        }

        return res;
    }

    /////////////// END OF HELPER FUNCTIONS  ///////////////


    /////////////// START OF COMPUTATION FUNCTIONS  ///////////////
    public Matrix solveforH(Matrix X, Matrix theta) {
        return X.times(theta);
    }

    public Matrix solveforHMinusY(Matrix H, Matrix Y) {
        return H.minus(Y);
    }

    public double summationOfHMinusYSquared(Matrix hminusy) {
        double res = 0;
        double curr;
        int row = hminusy.getRowDimension();
        int col = 1; //it is always 1; could be: "int col = hminusy.getColumnDimension()" but unnecessary
        for (int i = 0; i < row; i++) {
            for (int j = 0; j < col; j++) {
                curr = hminusy.get(i, j);
                res += curr * curr;
            }
        }
        return res;
    }

    public double cost(Matrix X, Matrix Y, Matrix theta) {
        Matrix H = solveforH(X, theta);
        Matrix hy = solveforHMinusY(H, Y);
        double sumofhysquared = summationOfHMinusYSquared(hy);
        double res = sumofhysquared / (2 * M); //according to the cost function formula: (( h(x) - y )^2 ) / (2m)
        return res;
    }

    public double cost() {
        return cost(this.X, this.Y, this.theta);
    }

    /////////////// END OF COMPUTATION FUNCTIONS  ///////////////
}