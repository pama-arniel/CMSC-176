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
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by Jose Arniel Pama on 10/16/2017.
 * Link to Java Matrices Library, JAMA: http://math.nist.gov/javanumerics/jama/
 */

public class Solution {
    private Matrix engineeredX; //the polynomially engineered values of X
    private Matrix multiclassY; //multiclass representation of the response output
    private Matrix multithetas; //theta values representing hypotheses
    ArrayList<String> YVals; //e.g. "Iris-setosa", "Iris-versicolor, "Iris-virginica"
    private double lambda;
    private int degree;
    private int thetavalset;
    private double alpha;
    private int iters;
    private int restarts;

    public Solution(double lambda, int degree, int thetavalset, double alpha, int iters, int restarts) {
        this.lambda = lambda;
        this.degree = degree;
        this.thetavalset = thetavalset;
        this.alpha = alpha;
        this.iters = iters;
        this.restarts = restarts;
    }

    public void load_data(String filename) {
        BufferedReader br = null;

        try {
            String currline;
            br = new BufferedReader(new FileReader(filename));

            ArrayList<ArrayList<Double>> X = new ArrayList<>();
            ArrayList<String> Y = new ArrayList<>();

            currline = br.readLine(); //so that the column labels will be skipped
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

                Y.add( (String) currentsample.get(len - 1)); //the last item, index len-1, is the Y of every training example

                tobeaddedtoX.add(0, 1.0); //for the X0, which is 1; necessary for matrix operation later
                X.add(tobeaddedtoX);
            }

            int rows = X.size();
            int cols = X.get(0).size();
            //store the contents of double arraylist X to a matrix
            Matrix xmatrix = new Matrix(rows, cols);
            for(int i = 0; i < rows; i++){
                ArrayList<Double> al = X.get(i);
                for(int j = 0; j < cols; j++){
                    xmatrix.set(i, j, al.get(j));
                }
            }

            System.out.println("NUMBER OF TRAINING EXAMPLES: " + xmatrix.getRowDimension());
            System.out.print("\nTHE X MATRIX:");
            xmatrix.print(xmatrix.getColumnDimension(), 2);

            System.out.println("THE OUTPUTS (Y):");
            for(String st: Y){
                System.out.println(st);
            }


            this.multiclassY = generateMultiClassY( Y );
            Matrix scaledX = scalefeatures(xmatrix);
            this.engineeredX = engineerPolynomials(scaledX, this.degree);
            this.multithetas = generateInitialThetas(this.thetavalset);
            regularizedCost(this.lambda);
            ArrayList<ArrayList<Double>> costhistories = gradientDescent( this.engineeredX, this.multiclassY, this.lambda, this.alpha, this.iters, this.restarts );
            graph(costhistories);

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

    public Matrix generateInitialThetas(int setval){
        int thetarows = engineeredX.getColumnDimension();
        int thetacols = multiclassY.getColumnDimension();

        Matrix res = new Matrix(thetarows, thetacols);

        for(int i=0; i < thetarows; i++){
            for(int j=0; j < thetacols; j++){
                res.set(i, j, setval);
            }
        }

        System.out.print("\nMODIFIED thetas MATRIX: ");
        res.print(res.getColumnDimension(), 2);
        return res;
    }

    public double[] generateRandomizedThetas(){
        int thetacols = engineeredX.getColumnDimension();
        double[] res = new double[thetacols];
        double currvalue = 0.0;

        for(int i=0; i < thetacols; i++){
            currvalue = (int)(-10 + Math.random()*10 ); //generate random values from -10 to 10
            res[i] = currvalue;

        }

        return res;
    }

    public Matrix generateMultiClassY(ArrayList<String> Y){
        //find first all the possible Y values in Y (e.g. "Iris-setosa", "Iris-versicolor, "Iris-virginica") and store them on a container
        ArrayList<String> possibleYVals = new ArrayList<>();
        int ylen = Y.size();
        for(String curr: Y){
            if(!possibleYVals.contains(curr)){
                possibleYVals.add(curr);
            }
        }

        this.YVals = possibleYVals;

        //start generating the multiclass matrix already
        int ycols = possibleYVals.size(); //number of possible Y values would dictate the number of Y columns
        Matrix res = new Matrix(ylen, ycols);
        for(int col = 0; col < ycols; col++){
            String currpossibleYVal = possibleYVals.get(col);

            for(int row = 0; row < ylen; row++){
                String currval = Y.get(row);

                if( currpossibleYVal.equals(currval) ){
                    res.set(row, col, 1);
                }else{
                    res.set(row, col, 0);
                }
            }
        }

        System.out.print("\nMultiple Class Y:");
        res.print(ycols, 2);
        return res;
    }

    public Matrix scalefeatures(Matrix X){
        int M = X.getRowDimension();
        int N = X.getColumnDimension();
        Matrix resmatrix = new Matrix(M, N);

        for(int col = 0; col < N; col++){
            Matrix Xj = X.getMatrix(0, M - 1, col, col);

            //calculate for the mean
            double currval;
            double total = 0;
            for(int i = 0; i < M; i++){
                currval = Xj.get(i, 0);
                total += currval;
            }

            double currmean = total / M;

            //calculate for the standard deviation
            double xmeansquared = 0; //summation of the (x - mean) squared
            double num;
            for(int i = 0; i < M; i++){
                currval = Xj.get(i, 0);
                num = currval - currmean;
                xmeansquared += (num * num);
            }
            double currstdev = Math.sqrt(xmeansquared / M);

            //scale values according to Z score
            double currzcore;
            for(int i = 0; i < M; i++){
                currval = Xj.get(i, 0);
                currzcore = (currval - currmean) / currstdev;

                if(col==0){ //meaning, it is at the bias column
                    resmatrix.set(i, col, 1.0);
                }else{
                    resmatrix.set(i, col, currzcore);
                }
            }
        }

        System.out.print("\nTHE SCALED X MATRIX:");
        resmatrix.print(N, 2);

        return resmatrix;
    }


    //returns the polynomially engineered values of X
    public Matrix engineerPolynomials(Matrix Xscaled, int degree){
        int someconst = degree + 1; //to be used in the calculations, e.g. calculating the row dimension of the resulting matrix
        int totalcols = Xscaled.getColumnDimension() - 1; //minus 1 for the bias
        int totalrows = (int) Math.pow(someconst, totalcols);
        Matrix exponentmatrix = new Matrix(totalrows, totalcols);

        int exponent = totalcols;
        int interval = 0;

        //generate first the table for exponents to be used in the final engineered polynomials generation
        for(int col = 0; col < totalcols; col++){
            exponent -= 1; //used to know the intervals of each to-be-inputted values which range from 0 to degree
            interval = (int) Math.pow( someconst , exponent ); //refers to how frequent a value (from 0 to degree) alternately appears in a certain column

            int input = 0;
            int num = 0;
            for(int row=0; row < totalrows; row++){
                if(num == interval){
                    if(input < degree){
                        input++;
                    }else{
                        input = 0;
                    }
                    num = 1;
                }else{
                    num++;
                }

                exponentmatrix.set(row, col, input);
            }
        }

        int scaledXRows = Xscaled.getRowDimension();
        int expRows = exponentmatrix.getRowDimension();

        //generate the final polynomially engineered X; use the exponentmatrix for computations
        Matrix resmatrix = new Matrix(scaledXRows, expRows);
        for(int scaledXrow = 0; scaledXrow < scaledXRows; scaledXrow++){
            for(int expRow = 0; expRow < expRows; expRow++){
                double res = 1.0;
                for(int col = 0; col < totalcols; col++){
                    res *= Math.pow( Xscaled.get(scaledXrow, col+1), exponentmatrix.get(expRow, col));
                }
                resmatrix.set(scaledXrow, expRow, res);
            }
        }

        System.out.print("\nPolynomially Engineered Values of X: ");
        resmatrix.print(resmatrix.getColumnDimension(), 2);
        return resmatrix;
    }

    public void regularizedCost(double lambda){

        int cols = this.multiclassY.getColumnDimension();
        int yrows = this.multiclassY.getRowDimension();
        int thetarows = this.multithetas.getRowDimension();

        System.out.println("REGULARIZED COSTS: ");
        for(int j = 0; j < cols; j++){
            Matrix currtheta = this.multithetas.getMatrix(0, thetarows-1, j, j);
            Matrix currentY = this.multiclassY.getMatrix(0, yrows-1, j, j);
            double res = regularizedCost(this.engineeredX, currentY, currtheta, lambda);
            System.out.println("COST " + (j + 1) + " : " + res);
        }

        System.out.println();
    }

    public double regularizedCost(Matrix engineeredX, Matrix Y, Matrix currtheta, double lambda){
        int M = engineeredX.getRowDimension();
        double constant = calcAppendedConst(M, currtheta, lambda, "cost");

        int xrows = engineeredX.getRowDimension();
        int xcols = engineeredX.getColumnDimension();

        //calculate now for the regularized cost for Y
        double summation = 0.0;
        for(int xrow=0; xrow < xrows; xrow++){
            double yval = Y.get(xrow, 0);
            double hOfX = calcForH(currtheta, engineeredX.getMatrix(xrow, xrow, 0, xcols-1));
            double gOfX = calcForG(hOfX);

            double a = yval * Math.log(gOfX);

            double b;
            double b1;

            b1 = Math.log(1 - gOfX);
            if (Double.isInfinite(b1) || Double.isNaN(b1)){
                b1 = 0.0;
            }

            b = (1 - yval) * b1;
            summation += ( a + b );
        }

        double res;

        try {
            //BigDecimal div = new BigDecimal((-summation/M) );
            //div = div.setScale(2, BigDecimal.ROUND_HALF_UP);
            BigDecimal a = new BigDecimal(summation).negate();
            BigDecimal b = new BigDecimal(M);
            BigDecimal div = a.divide(b,2, BigDecimal.ROUND_HALF_UP);
            res = div.doubleValue() + constant;
        }catch (NumberFormatException e){
            res = Double.NaN;
        }

        return res;
    }

    //calculate first for the "appended constant"; depends on type supplied
    public double calcAppendedConst(int M, Matrix currtheta, double lambda, String type){
        double res = 0.0;
        double currval;
        int thetarows = currtheta.getRowDimension();

        if(type.equals("cost")){
            for(int thetarow = 0; thetarow < thetarows; thetarow++){
                currval = currtheta.get(thetarow,0);
                res += (currval*currval);
            }
            res = res * (lambda / (2 * M));
        }else if (type.equals("gradient descent")){
            for(int thetarow = 0; thetarow < thetarows; thetarow++){
                currval = currtheta.get(thetarow,0);
                res += currval;
            }
            res = res * (lambda /  M);
        }else{
            System.out.println("Invalid type provided");
        }

        return res;
    }

    public double calcForH(Matrix theta, Matrix givenX){
        Matrix resmatrix = givenX.times(theta);
        return resmatrix.get(0,0);
    }

    public double calcForG(double hOfXvalue){
        return (1 / (1 + Math.pow(Math.E, -hOfXvalue)));
    }


    //return the cost histories of each Y column
    public ArrayList<ArrayList<Double>> gradientDescent(Matrix engineeredX, Matrix multiclassY, double lambda, double alpha, int iters, int restarts) {
        ArrayList<ArrayList<Double>> res = new ArrayList<>();
        int ycols = multiclassY.getColumnDimension();
        int yrows = multiclassY.getRowDimension();
        int xcols = engineeredX.getColumnDimension();

        int M = engineeredX.getRowDimension();
        int N = engineeredX.getColumnDimension(); //number of Xs in one training example

        System.out.println("GRADIENT DESCENT FINAL COSTS: ");
        for(int ycol = 0; ycol < ycols; ycol++){
            ArrayList<Double> mincosthistory = new ArrayList<>();
            double mincost = Double.MAX_VALUE;
            double[] thetaArr = new double[N];

            for(int restart = 1; restart <= restarts; restart++){
                ArrayList<Double> costhistory = new ArrayList<>();
                Matrix theta;
                thetaArr = generateRandomizedThetas();
                double constant;

                for (int i = 1; i <= iters; i++) {
                    theta = Matrix.constructWithCopy(generateTheta(thetaArr));
                    constant = calcAppendedConst(M, theta, lambda, "gradient descent");
                    Matrix currY = multiclassY.getMatrix(0, yrows - 1, ycol, ycol);
                    double currcost = regularizedCost(engineeredX, currY, theta, lambda);
                    costhistory.add( currcost );

                    for(int xcol=0; xcol< xcols; xcol++){
                        thetaArr[xcol] = thetaArr[xcol] - (alpha * derivativeJ(engineeredX, currY, theta, xcol)) + constant;
                    }
                }

                double lastcost = costhistory.get(iters-1);
                if( lastcost < mincost){
                    mincosthistory = costhistory;
                    mincost = lastcost;
                }
            }

            System.out.println("COST " + (ycol + 1) + " : " + mincost);
            res.add(mincosthistory);
        }

        System.out.println();
        return res;
    }

    public void graph(ArrayList<ArrayList<Double>> costhistories){
        int index = 0;
        for (ArrayList<Double> costhistory : costhistories){
            new Graph(costhistory, YVals.get(index++));
        }
    }

    public class Graph extends ApplicationFrame {
        public Graph(ArrayList<Double> costhistory, String title) {
            super("Gradient Descent");

            final XYSeries series = new XYSeries("Cost");
            int iters = costhistory.size();
            double currnum = 0.0;
            for(int i = 0; i < iters; i++){
                currnum = costhistory.get(i);
                series.add(i, currnum);
            }

            //System.out.println("LAST CALCULATED COST: " + costhistory.get(iters-1));
            final XYSeriesCollection data = new XYSeriesCollection(series);
            final JFreeChart chart = ChartFactory.createXYLineChart(
                    title + " cost over iteration",
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

    public double[][] generateTheta(double[] given) {
        int xrows = given.length;
        double[][] res = new double[xrows][];
        for (int i = 0; i < xrows; i++) {
            double[] newrow = new double[1]; //provided theta has 1 column always
            newrow[0] = given[i];
            res[i] = newrow;
        }

        return res;
    }

    //return the derivative of cost J given Matrix X, Y, and theta
    public double derivativeJ(Matrix engineeredX, Matrix currY, Matrix currtheta, int j){
        Matrix hOfX = engineeredX.times(currtheta);
        int rows = hOfX.getRowDimension();
        Matrix gOfX = new Matrix(rows, 1);
        for(int i = 0; i < rows; i++){
            gOfX.set(i, 0, calcForG(hOfX.get(i, 0)));
        }

        Matrix hy = gOfX.minus(currY);
        Matrix Xj = engineeredX.getMatrix(0, engineeredX.getRowDimension() - 1, j, j);
        Matrix derJmatrix = Xj.transpose().times(hy);

        double res = derJmatrix.get(0, 0); //because derJmatrix is scalar
        return (res / engineeredX.getRowDimension());
    }
}
