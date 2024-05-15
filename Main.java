import java.util.Arrays;
import java.util.Objects;
import java.util.Scanner;

import static java.lang.System.exit;


public class Main {

    public static double [][] MatrixTranspon (double [][] A) {//Транспонирование матрицы
        double [][] TransponA = new double[A[0].length][A.length];
        for (int i = 0; i < A.length; i++) {
            for (int j = 0; j < A[0].length; j++) {
                TransponA[j][i] = A[i][j];
            }
        }
        return TransponA;
    }

    public static double [][] MatrixMultiply (double [][] A, double[][] B) {//Умножение матриц
        double [][] newMatrix = new double[A.length][B[0].length];
        for (var i = 0; i < newMatrix.length; i++) {
            for (var j = 0; j < newMatrix[0].length; j++) {
                newMatrix[i][j] = 0;
                for (var k = 0; k < A[0].length; k++) {
                    newMatrix[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        return newMatrix;
    }
    public static double [][] MatrixInversion(double [][]MatrixA, int N) {
        double temp;
        double[][] A = new double[N][N];
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                A[i][j] = MatrixA[i][j];
            }
        }
        double[][] E = new double[N][N];


        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++) {
                E[i][j] = 0f;

                if (i == j)
                    E[i][j] = 1f;
            }

        for (int k = 0; k < N; k++) {
            temp = A[k][k];

            for (int j = 0; j < N; j++) {
                A[k][j] /= temp;
                E[k][j] /= temp;
            }

            for (int i = k + 1; i < N; i++) {
                temp = A[i][k];

                for (int j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }

        for (int k = N - 1; k > 0; k--) {
            for (int i = k - 1; i >= 0; i--) {
                temp = A[i][k];

                for (int j = 0; j < N; j++) {
                    A[i][j] -= A[k][j] * temp;
                    E[i][j] -= E[k][j] * temp;
                }
            }
        }

        for (int i = 0; i < N; i++)
            for (int j = 0; j < N; j++)
                A[i][j] = E[i][j];
        return A;
    }
    public static double [][] IdentityMatrix (int n) {//Едининичная матрица
        double [][] IdentityMatrix = new double[n][n];
        for (int i = 0; i < n; i++){
            IdentityMatrix[i][i] = 1;
        }
        return IdentityMatrix;
    }
    public static double  Norma (double [][] A) {//Транспонирование матрицы
        double  Norma = 0;
        for (int i = 0; i < A.length; i++) {
            Norma+= Math.pow(A[i][0],2);
        }
        Norma=Math.sqrt(Norma);
        return Norma;
    }
    public static double Function (double [] X) { //Минимизируемая функция
        return 3*Math.pow(X[0],2) + 4*X[0]*X[1] +  5*Math.pow(X[1],2);
    }
    public static double [][] Limitations (double [] X, int m) { //Вектор значений ограничений в точке
        double[][] LinitationsValue = new double[1][m];
        LinitationsValue [0][0] = -X[0];
        LinitationsValue [0][1] = -X[1];
        LinitationsValue [0][0] = -X[0] - X[1] + 4;
        return LinitationsValue;
    }
    public static double Tk1 (double [] X, double[] DeltaXk) { //Вычесление Tk
        double Tk1 = X[0]/-DeltaXk[0];
        return Tk1;
    }
    public static double Tk2 (double [] X, double[] DeltaXk) { //Вычесление Tk
        double Tk2 = X[1]/-DeltaXk[1];
        return Tk2;
    }
    public static double Tk3 (double [] X, double[] DeltaXk) { //Вычесление Tk
        double Tk3 = -(-1+X[0]+X[1])/(DeltaXk[0]+DeltaXk[1]);
        return Tk3;
    }
    public static double [][] A (double [] X, int m) {//Матрица первых производных ограничений
        double [][] A = new double[m][X.length];
        A[0][0]=-1;
        A[0][1]=0;
        A[1][0]=0;
        A[1][1]=-1;
        A[2][0]=-1;
        A[2][1]=-1;
        return A;
    }

    public static double [][] Gradient (double [] X, int n) { //Градиенты
        double[][] GradientValue = new double[n][1];
        GradientValue [0][0] = 6*X[0] + 4*X[1];
        GradientValue [1][0] = 4*X[0] + 10*X[1];
        return GradientValue;
    }

    public static double [] Lagranzh (double [] X, double [] Lambda) { //Невязки множителей Лагранжа
        double[] LagranzhValue = new double[X.length + Lambda.length];
        LagranzhValue [0] = 2*(X[0]-4) + Lambda[0];
        LagranzhValue [1] = 2*(X[1]-5) + 2*Lambda[0]*X[1];
        LagranzhValue [2] = X[0] + Math.pow(X[1],2) - 4;
        return LagranzhValue;
    }
    public static double ValueFunction(double[] X, double[] S, double x) {
        double F;
        F = 3*Math.pow(X[0] + x*S[0], 2) +4*(X[0] + x*S[0])*(X[1]+x*S[1]) + 5*Math.pow(X[1]+x*S[1], 2);
        return F;
    }
    public static double Firstderivative(double[] X, double[] S, double x) {
        double P = 0.000001;
        return (ValueFunction(X, S, x + P) - ValueFunction(X, S, x)) / P;
    }

    public static double Secondderivative(double[] X, double[] S, double x) {
        double P = 0.000001;
        return (Firstderivative(X, S, x + P) - Firstderivative(X, S, x)) / P;

    }
    public static void Newton(double E, double[] X, double[] S) {
        boolean z = false;
        double x = 5;
        double x_0 = x;
        do {
            x = x_0 - Firstderivative(X, S, x_0) / Secondderivative(X, S, x_0);
            if (Math.abs(x - x_0) > E /*|| Math.Abs(Function(x[i + 1]) - Function(x[i])) < constE*/) {
                z = true;
            } else
                break;
            x_0 = x;
        } while (z);
        X[0] = X[0] + x * S[0];
        X[1] = X[1] + x * S[1];
    }

    //gradient projections
    public static double [] GradientProjectionsEquality (double [] X, double Epsilon, int n, int m, int k, String EquationType) {
        double[][] A = A(X,m); //Шаг 4
        double[][] Tau = Limitations(X,m); //Шаг 5
        for (int j = 0; j < m; j++){
            Tau[0][j]=-Tau[0][j];
        }
        double[][] Delta2Xk = new double[1][n]; //Шаг 6
        Delta2Xk = MatrixMultiply(MatrixMultiply(MatrixTranspon(A), MatrixInversion(MatrixMultiply(A, MatrixTranspon(A)), m)),MatrixTranspon(Tau));
        double NormaDelta2Xk; //Шаг 7
        if (k>=1 & Objects.equals(EquationType, "linear")){
            NormaDelta2Xk = 0;
        }else {
            NormaDelta2Xk = Norma(Delta2Xk);
        }
        double[][] FunctionGradient = Gradient(X, n); //Шаг 8
        double[][] E = IdentityMatrix(n); //Шаг 9
        double[][] DeltaXk = new double[n][1];
        double[][] SlojnoeProizvedenie = MatrixMultiply(MatrixMultiply(MatrixTranspon(A),MatrixInversion(MatrixMultiply(A, MatrixTranspon(A)), m)),A);
        for (int i = 0; i < n; i++){
            for (int j = 0; j < n; j++){
                E[i][j]=-(E[i][j] - SlojnoeProizvedenie[i][j]);
            }
        }
        DeltaXk = MatrixMultiply(E,FunctionGradient);
        double NormaDeltaXk = Norma(DeltaXk);
        boolean Condition = false;
        if (NormaDeltaXk <= Epsilon & NormaDelta2Xk <= Epsilon){ //Шаг 10
            double[] newX = new double[n+1+m];
            for (int i = 0; i < n; i++){
                newX[i] = X[i];
            }
            newX[n] = 1;
            SlojnoeProizvedenie = MatrixMultiply(MatrixMultiply(MatrixInversion(MatrixMultiply(A, MatrixTranspon(A)),m),A), FunctionGradient);
            for (int i = 0; i < m; i++){
                newX[n+1+i] = -SlojnoeProizvedenie[i][0];
            }
            return newX;
        }else if (NormaDeltaXk > Epsilon & NormaDelta2Xk <= Epsilon){
            for (int i = 0; i < n; i++){
                Delta2Xk[i][0]=0;
            }
        }else if (NormaDeltaXk <= Epsilon){
            for (int i = 0; i < n; i++){
                DeltaXk[i][0]=0;
            }
            Condition = true;
        }
        if (!Condition) {
            double[] TransponDeltaXk = new double[n];
            for (int i = 0; i < n; i++) {
                TransponDeltaXk[i] = DeltaXk[i][0];
            }
            Newton(0.01, X, TransponDeltaXk); //Шаг 11, 12
        }
        double[] newX = new double[n + 1];
        for (int i = 0; i < n; i++){ //Шаг 13
            newX[i] = X[i]+Delta2Xk[i][0];
        }
        return newX;
    }

    public static double [] GradientProjectionsInequality (double [] X, double Epsilon1, double Epsilon, int n, int m, int k) {
        double[][] AllA = A(X,m);
        double[][] Tau = Limitations(X,m); //шаг 4
        boolean gCondition = false;
        boolean[] LimitationIsActive = new boolean[m];
        int iterator = 0 ;
        int NumberOfPassive =  m;
        for (int i = 0; i < m; i++){ //шаг 5
            if (Math.abs(Tau[0][i]) < 0.00001){
                LimitationIsActive[i] = true;
                iterator++;
                NumberOfPassive--;
            }
            if(Tau[0][i]>=Epsilon1 & Tau[0][i]<=0){
                gCondition = true;
            }
        }
        double[][] A = new double[iterator][n];
        int s = 0;
        for (int i=0; i<m; i++){ // считаем матрицу состоящую ТОЛЬКО из активных ограничений
            if (LimitationIsActive[i]) {
                for (int j = 0; j < n; j++) {
                    A[s][j] = AllA[i][j];
                }
                s++;
            }
        }
        double[][] FunctionGradient = Gradient(X, n); // Считаем градиент
        if (gCondition & Norma(FunctionGradient)==0 & k==0) { //шаг 5 а) при F(x0)=0
            return new double[]{Math.PI};
        }
        if (!gCondition){ //шаг 6 при невыполненении условий из шага 4
            for (int j = 0; j < m; j++){
                Tau[0][j]=-Tau[0][j];
            }
            double[][] SlojnoeProizvedenie = MatrixMultiply(MatrixMultiply(MatrixTranspon(A),MatrixInversion(MatrixMultiply(A, MatrixTranspon(A)), m)),Tau); // Шаг 6
            for (int i = 0; i < m; i++) {
                X[i]+= SlojnoeProizvedenie[i][0];
            }
        }
        double[][] DeltaXk = new double[n][1];
        double[] newX = new double[n+1+iterator];
        do {
            int DeleteRowIndex = -1;
            double Lambdamin = 0;
            double[][] E = IdentityMatrix(n); // Шаг 7
            double[][] SlojnoeProizvedenie = MatrixMultiply(MatrixMultiply(MatrixTranspon(A),MatrixInversion(MatrixMultiply(A, MatrixTranspon(A)), iterator)),A);
            for (int i = 0; i < n; i++){
                for (int j = 0; j < n; j++){
                    E[i][j]=-(E[i][j] - SlojnoeProizvedenie[i][j]);
                }
            }
            DeltaXk = MatrixMultiply(E,FunctionGradient); //Шаг 8
            double NormaDeltaXk = Norma(DeltaXk);
            if(NormaDeltaXk < Epsilon) { //Шаг 9
                double[] Lambda = new double[iterator];
                SlojnoeProizvedenie = MatrixMultiply(MatrixMultiply(MatrixInversion(MatrixMultiply(A, MatrixTranspon(A)), iterator), A), Gradient(X, n));
                for (int i = 0; i < iterator; i++) {
                    Lambda[i] = -SlojnoeProizvedenie[i][0];
                    if (Lambda[i]<Lambdamin){
                        Lambdamin = Lambda[i];
                        DeleteRowIndex=i;
                    }
                }
                if (DeleteRowIndex==-1){ //Вываливаемся если все лямбда положительны
                    for (int i = 0; i < n; i++){
                        newX[i] = X[i];
                    }
                    newX[n]=1;
                    for (int i = 0; i < iterator; i++){
                        newX[n+1+i] = Lambda[i];
                    }
                    return newX;
                }
                else {
                    iterator--;
                    AllA = A;
                    A = new double[iterator][n];
                    s=0;
                    for (int i=0; i<iterator; i++){ // Убираем одно из ограничений
                        if (i == DeleteRowIndex) {
                            s=1;
                        }
                        for (int j = 0; j < n; j++) {
                            A[i][j] = AllA[i+s][j];
                        }
                    }
                }
            }else {
                break;
            }
        }while (true);
        double[] t = new double[1+NumberOfPassive];
        double[] TransponDeltaXk = new  double[n];
        double[] TryX = new  double[n];
        for (int i = 0; i < n; i++) {
            TransponDeltaXk[i] = DeltaXk[i][0];
            TryX[i] = X[i];
        }
        Newton(0.01, TryX, TransponDeltaXk);
        for (int i = 0; i < n; i++) { //Найдем t-катую со звездой
            if(TryX[i]!=0){
                t[0]=(TryX[i]-X[i])/TransponDeltaXk[i];
            }
        }
        for (int j=0; j<NumberOfPassive; j++){ // Найжем t-катые из пасивных ограничений
            for (int i = 0; i < n; i++) {
                if (!LimitationIsActive[i]){
                    LimitationIsActive[i] = true;
                    if(i==0) {
                        t[j + 1] = Tk1(X,TransponDeltaXk);
                    } else if (i==1) {
                        t[j + 1] = Tk2(X,TransponDeltaXk);
                    }else {
                        t[j + 1] = Tk3(X,TransponDeltaXk);
                    }
                }
            }
        }
        double minTk = Double.MAX_VALUE;
        for (int i = 0; i < NumberOfPassive+1; i++) {
            if (t[i]>0 & t[i]<minTk){
                minTk=t[i];
            }
        }
        for (int i = 0; i < n; i++) {
            newX[i]=X[i]+minTk*TransponDeltaXk[i];
        }
        return newX;
    }

    public static void main(String[] args) {
        double[] X = {0,4};
        double Epsilon = 0.1;
        double Epsilon1 = -0.1;
        int M = 20;
        int n = 2;
        int m = 3;
        int k = 0;
        String EquationType = "notlinear";
        String EquatiunSymbol = "<";
        double[] Lambda = new double[m];
        if (EquatiunSymbol.equals("=")) {
            do {
                if (k >= M) {
                    double[][] SlojnoeProizvedenie = MatrixMultiply(MatrixMultiply(MatrixInversion(MatrixMultiply(A(X, m), MatrixTranspon(A(X, m))), m), A(X, m)), Gradient(X, n));
                    for (int i = 0; i < m; i++) {
                        Lambda[i] = SlojnoeProizvedenie[i][0];
                    }
                    break;
                }
                double[] newX = GradientProjectionsEquality(X, Epsilon, n, m, k, EquationType);
                for (int i = 0; i < n; i++) {
                    X[i] = newX[i];
                }
                if (newX[n] == 1) {
                    for (int i = 0; i < m; i++) {
                        Lambda[i] = newX[n + 1 + i];
                    }
                    break;
                }
                k++;
            } while (true);
            double[] Nevyazki = Lagranzh(X, Lambda);
            boolean TheNecessaryConditions = true;
            for (int i = 0; i < Nevyazki.length; i++) {
                if (Math.abs(Nevyazki[i]) > Epsilon) {
                    TheNecessaryConditions = false;
                    break;
                }
            }
            if (TheNecessaryConditions) {
                System.out.println("Необходимые условия минимума выполняются");
            } else {
                System.out.println("Необходимые условия минимума НЕ выполняются");
            }
            System.out.print("\nОптимальная точка x* = { ");
            for (int i = 0; i < n; i++) {
                System.out.printf("%.2f", X[i]);
                if (i != n - 1) {
                    System.out.print("; ");
                }
            }
            System.out.print("}");
        }
        else {
            do {
                if (k >= M) {
                    double[][] SlojnoeProizvedenie = MatrixMultiply(MatrixMultiply(MatrixInversion(MatrixMultiply(A(X, m), MatrixTranspon(A(X, m))), m), A(X, m)), Gradient(X, n));
                    for (int i = 0; i < m; i++) {
                        Lambda[i] = SlojnoeProizvedenie[i][0];
                    }
                    break;
                }
                double[] newX = GradientProjectionsInequality(X, Epsilon1, Epsilon, n, m, k);
                if(newX[0]==Math.PI){
                    System.out.print("Введите другую начальную точку");
                    exit(0);
                }
                for (int i = 0; i < n; i++) {
                    X[i] = newX[i];
                }
                if (newX[n] == 1) {
                    for (int i = 0; i < newX.length-n-1; i++) {
                        Lambda[i] = newX[n + 1 + i];
                        System.out.println("Lambda" + i + " = " + Math.abs(Lambda[i]));
                    }
                    System.out.print("\nОптимальная точка x* = { ");
                    for (int i = 0; i < n; i++) {
                        System.out.printf("%.2f", X[i]);
                        if (i != n - 1) {
                            System.out.print("; ");
                        }
                    }
                    System.out.print("}");
                    break;
                }
            }while (true);
        }
//        for (int i = 0; i < newA.length; i++) {
//            for (int j = 0; j < newA[0].length; j++) {
//                System.out.print(newA[i][j]+"\t");
//            }
//            System.out.println("");
//        }
    }

}