package mathlib;

import java.util.Random;

/**
 * Esta clase almacena un conjunto de funciones matemáticas útiles para la estimación
 * y pronóstico del SOC.  
 * @author Alex Díaz
 * @vesion 1.3
 */

public class Matlab {
	
	/**
	 * Método que retorna una matriz de ceros
	 * @param n = filas
	 * @param m = columnas 
	 * @return Matriz de ceros
	 */
	
	public static double[][] zeros(int rows, int columns) {
		double[][] matrixZeros = new double [rows][columns];
		for (int i = 0; i < rows; i++) { 			//filas
			for (int j = 0; j < columns; j++) {		//columnas
				matrixZeros[i][j] = 0;
			}
		}
		return matrixZeros;

	}
	
	/**
	 * Método que redondea un número
	 * @param numero double a redondear
	 * @param decimales int número de decimales a redondear 
	 * @return double redondeado
	 */
	public static double round( double number, int decimals ) {
	    return Math.round(number*Math.pow(10,decimals))/Math.pow(10,decimals);
	}
	
	/**
	 * Método que retorna una muestra de una matriz
	 * @param v double[][] matriz a muestrear
	 * @param i int inicio del muestreo
	 * @param cada int espacio entre cada muestra
	 * @param f int final del muestreo
	 * @return muestra de una matriz del ancho de la matriz v y largo de (f-i+1/cada)
	 */
	public static double[][] sample(double[][] _matrix, int initiation, int each, int fin) {
		float finalLength 		= (float) (fin-initiation+1/each);
		double[][] matrix	= new double [(int) finalLength][_matrix[0].length];
		int j 				= initiation;
		
		for(int l = 0; l<_matrix[0].length; l++){
			for(int k = 0; k<matrix.length; k++){
				matrix[k][l] = _matrix[j-1][l];
				j = j + each;
			}
			j = initiation;
		}
		return matrix;
	}
	
	/**
	 * Método que retorna la matriz transpuesta de una matriz de n x m
	 * @param mat double[][] matriz a transponer
	 * @return mat_t double[][] transpuesta de m x n de mat
	 */
	public static double[][] transpose(double[][] matrix) {
		double [][] transposedMatrix = new double [matrix[0].length][matrix.length];
		for (int i = 0 ; i < matrix.length ; i++){
            for (int j = 0 ; j < matrix[0].length ; j++){
            	transposedMatrix[j][i] = matrix[i][j];
            }
        }
		return transposedMatrix;
	}
	
	/**
	 * Método que retorna una matriz de n x m rellena de num
	 * @param num float número a rellenar
	 * @param n int número de filas
	 * @param m int número de columnas
	 * @return matriz double[][] de n x m
	 */
	public static double[][] repeatArray(float number, int rows, int columns) {
		double[][] matrix = new double [rows][columns];
		for(int i = 0; i<rows; i++){
			for(int j = 0; j<columns; j++){
				matrix[i][j] = round(number, 4);
			}
		}
		return matrix;
	}
	
	/**
	 * Método que retorna una matriz de n x m rellena con números aleatorios entre 0 y 1
	 * @param n int número de filas
	 * @param m int número de columnas
	 * @return mat_rand double[][] de n x m
	 */
	public static double[][] rand(int rows, int columns) {
		double[][] randomMatrix = new double [rows][columns];
		for (int i = 0; i<rows; i++){
			for (int j = 0; j < columns; j++) {
				randomMatrix[i][j] = Math.random();
			}
		}
		return randomMatrix;
	}
	
	/**
	 * Método que retorna la suma de una fila o columna de una matriz	
	 * @param mat double[][] matriz de entrada 
	 * @param k int indica fila o columna a sumar
	 * @param roc String indica si se quiere sumar una "fila" o "columna"
	 * @return suma double retorna la suma total
	 */
	public static double sum(double[][] matrix, int k, String roc) {
		double sum = 0;
		if(roc.equals("columna")){
			for (int i = 0; i<matrix.length; i++){ 
				sum = sum+matrix[i][k];
			}
		} else if(roc.equals("fila")){
			for (int i = 0; i<matrix[0].length; i++){
				sum = sum+matrix[k][i];
			}
		}
		return sum;
	}
	
	/**
	 * Método que retorna la suma del cuadrado de cada celda de un vector
	 * @param mat double[][] matriz de entrada
	 * @param k int posición de la fila o columna
	 * @param roc String indica fila o columna.
	 * @return suma double retorna la suma del cuadrado de cada celda de un vector
	 */
	public static double sum2(double[][] matrix, int k, String roc) {
		double sum = 0;
		if(roc.equals("fila")){
			for (int j = 0; j<matrix[0].length; j++){ 
				sum = sum+matrix[k][j]*matrix[k][j];
			}
		} else if(roc.equals("columna")){
			for (int i = 0; i<matrix.length; i++){ 
				sum = sum+matrix[i][k]*matrix[i][k];
			}
		}
		return sum;
	}
	
	/**
	 * Método que retorna un número aleatorio con distribución normal
	 * @param mu double media
	 * @param sigma double desviación estandar
	 * @return gaussians double número aleatorio
	 */
	public static double normrnd(double mu, double sigma){
		Random r = new Random();
        double gaussians = r.nextGaussian()*sigma+mu;
		return gaussians;
	}
	
	/**
	 * Función densidad de probabilidad Normal (pdf)
	 * @param x double valor a ser evaluado
	 * @param mu double media
	 * @param sigma double desviación estandar
	 * @return y double Retorna la pdf de distribución normal con media mu y desviación estandar
	 * sigma, evaluada con x
	 */
	public static double normpdf(double x, int mu, double sigma) {
		//BigDecimal num = new BigDecimal((1/(sigma*Math.sqrt(2*Math.PI)))*Math.exp((-Math.pow(x-mu, 2))/(2*Math.pow(sigma, 2))));
		//float aux1 = (float) (1/(sigma*Math.sqrt(2*Math.PI)));
		//float aux2 = (float) (-Math.pow(x-mu, 2));
		//float aux3 = (float) (2*Math.pow(sigma, 2));
		//float aux4 = aux2/aux3;
		//float aux2 = (float) (-Math.pow(x-mu, 2)/(2*Math.pow(sigma, 2)));
		//double aux5 = Math.exp(aux2);
		//double y = aux1*aux5;
		float y = (float) ((1/(sigma*Math.sqrt(2*Math.PI)))*Math.exp((-Math.pow(x-mu, 2)/(2*Math.pow(sigma, 2)))));
		return (double) y;
	}
	
	/**
	 * Crea una matriz de valores aleatorios con distribución normal. 
	 * @param rows int número de filas
	 * @param columns int número de columnas
	 * @return m double[][] matriz de valores aleatorios
	 */
	public static double[][] randn(int rows, int columns) {
		double[][] matrix = new double[rows][columns];
		java.util.Random r = new java.util.Random();
		for (int i = 0; i < rows; i++) {
			for (int j = 0; j < columns; j++) {
				double sign = Math.random();
				if (sign<0.5){
					matrix[i][j] = r.nextGaussian();
				} else {
					matrix[i][j] = -r.nextGaussian();
				}
			}
		}
		return matrix;
	}
	
	/**
	 * Devuelve un vector columna de muestras ponderadas tomadas con reemplazo si 
	 * "reemplazo" es true, o sin reemplazo si "reemplazo" es false. Se toman "ctd" 
	 * de muestras entre "ini" y "fin". Se utilizan pesos positivos donde "wk" es
	 * un vector de probabilidades. 
	 * @param ini int límite inferior de muestras
	 * @param fin int límite superior de muestras
	 * @param numberSamples int cantidad de muestras
	 * @param reemplazo boleean indica reemplazo
	 * @param wk double[][] vector de probabilidades o pesos
	 * @return v int[][] vector columna de muestras
	 */
	public static int[][] randsample(int ini, int fin, int numberSamples, double[][] wk) {
		int[][] matrix		= new int[numberSamples][1];
		double sumWeights 	= 0;
		double random 	 	= 0;
			
		for (int i = 0; i < numberSamples; i++){
			random 		= Math.random();
			sumWeights = 0;
			for (int k = 0; k < numberSamples; k++) {
				sumWeights += wk[k][0];
				if (random <= sumWeights) {
					matrix[i][0] = k+ini;
					break;
				}
			}
		}
		
		return matrix;
	}
	
	/**
	 * Método que entrega una matriz de n x m de unos
	 * @param n int número de filas
	 * @param m int número de columnas
	 * @return matriz_unos double[][] matriz de unos
	 */
	public static double[][] ones(int rows, int columns){
		double[][] matrixOnes = new double [rows][columns];
		for (int i = 0; i < rows; i++) { 			//filas
			for (int j = 0; j < columns; j++) {		//columnas
				matrixOnes[i][j] = 1;
			}
		}
		return matrixOnes;
	}
	
	/**
	 * Función que encuentra el máximo de un vector columna
	 * @param data double[][] vector columna 
	 * @return max double máximo valor del vector
	 */
	public static double max(double[][] data) {
		double max = data[0][0];
		for (int i = 1; i<data.length; i++){
			if(data[i][0]>max){
				max = data[i][0];
			}
		}
		return max;
	}
	
	/**
	 * Función que encuentra el mínimo de un vector columna
	 * @param vec double[][] vector columna 
	 * @return min double minimo valor del vector
	 */
	public static double min(double[][] data) {
		double min = data[0][0];
		for (int i = 1; i<data.length; i++){
			if(data[i][0]<min){
				min = data[i][0];
			}
		}
		return min;
	}
	
	/**
	 * Corta un vector columna entre los índices inicio y fin 
	 * @param vector double[][] vector columna de entrada
	 * @param inicio double
	 * @param fin
	 * @return
	 */
	public static double[][] cutVector(double[][] vector, int initiation, int fin) {
		double[][] r = new double[fin-initiation+1][1];
		for(int i = 0; i<fin-initiation+1; i++){
			r[i][0] = vector[initiation-1+i][0];
		}
		return r;
	}
	
	
	/*public static double[][] cut_col(double[][] m, int j) {
		double[][] r = new double[m.length][1];
		for(int i=0; i<m.length; i++) {
			r[i][0] = m[i][j-1];
		}
		return r;
	}*/
	
	/**
	 * Función que reemplaza el vector de entrada según los indices de entrada
	 * @param _matrix double[][] matriz de entrada de 2 filas
	 * @param index int[][] vector columna que contiene los indices
	 * @return matrix double[][] matriz de 2 filas
	 */
	public static double[][] extract(double[][] _matrix, int[][] index) {
		double[][] matrix = new double[_matrix.length][_matrix[0].length];
		for (int i = 0 ; i < _matrix[0].length ; i++){
			matrix[0][i] = _matrix[0][index[i][0]-1];
			matrix[1][i] = _matrix[1][index[i][0]-1];
		}
		return matrix;
	}
	
	/**
	 * Encuentra el primer índice de un vector fila menor que c  
	 * @param matrix double[][] vector fila de entrada
	 * @param threshold double umbral
	 * @return index int índice del primer valor menor que c
	 */
	public static int find(double[][] matrix, double threshold) {
		int index = -1;
		
		for (int k = 0; k<matrix[0].length; k++){
			if (matrix[0][k]<threshold){
				index = k;
				break;
			}
		}
		return index;
	}
	
	/**
	 * Encuentra el primer índice de un vector columna mayor o igual a threshold
	 * @param matrix double[][] vector columna de entrada
	 * @param threshold double umbral
	 * @return r int índice del primer valor mayor o igual a e
	 */
	public static int find_ma_i(double[][] matrix, double threshold) {
		int index = -1;
		//int largo = d.length;
		for (int k = 0; k<matrix.length; k++){
			if (matrix[k][0] >= threshold){
				index = k;
				break;
			}
		}
		return index;
	}
	
	/**
	 * Función que calcula el average de un vector columna entre un indice de inicio y final
	 * @param vector double[][] vector columna
	 * @param i int índice inicial
	 * @param f int índice final
	 * @return v double average
	 */
	public static double mean(double[][] vector, int initiation, int fin) {
		float average = 0;
		float numberData = fin-initiation+1;
		for (int k = initiation ; k < fin+1 ; k++){
			average = (float) (average + vector[k][0]);
		}
		average = average/numberData;		
		return (double) average;
	}
	
	/**
	 * Suma acumulativa de elementos de un vector
	 * @param vector double[][] vector de entrada
	 * @param k int posición
	 * @param roc String indica fila o columna
	 * @return m_cumsum double[][] suma acumulada del vector de entrada
	 */
	public static double[][] cumsum(double[][] vector, int k, String roc) {
		if(roc.equals("columna")){
			double[][] m_cumsum = new double [vector.length][1];
			//int columna = vector.length;
			for (int i = 0; i < vector.length; i++) {
				if(i==0){
					m_cumsum[i][0] = vector[i][k];
				}
				else {
					m_cumsum[i][0] = vector[i][k] + m_cumsum[i-1][0];
				}
			}
			return m_cumsum;
		}
		else if (roc.equals("fila")){
			double[][] m_cumsum = new double [1][vector[0].length];
			//int fila = vector.length;
			for (int i = 0; i < vector[0].length; i++) {
				if(i==0){
					m_cumsum[i][0] = vector[k][i];
				}
				else {
					m_cumsum[i][0] = vector[k][i] + m_cumsum[i-1][0];
				}
			}
			return m_cumsum;
		} else {
			return null;
		}
	}
	
	/**
	 * Función que crea una matriz de dos columnas
	 * @param X double[][] vector columna
	 * @param Y double[][] vector columna
	 * @return r double[][] matriz de dos columnas [X Y]
	 */
	public static double[][] createKernel(double[][] vectorX, double[][] vectorY) {
		double[][] matrix = new double[vectorX.length][2];
		for(int i = 0; i<vectorX.length; i++){
			for(int j=0; j<2; j++){
				if(j==0){
					matrix[i][j]	= vectorX[i][0];
				} else {
					matrix[i][j]	= vectorY[i][0];
				}
			}
		}
		return matrix;
	}
	
	/**
	 * Función que genera un vector con los elementos dentro de un rango
	 * @param ini double límite inferior del rango
	 * @param d double distancia entre cada elemento
	 * @param fin double límite superior del rango
	 * @return r double[][] retorna un vector columna
	 */
	public static double[][] createRange(double lowerLimit, double distance, double upperLimit) {
		float difference	= (float) ((upperLimit-lowerLimit)/distance);
		double[][] vector	= new double[(int) (difference+1)][1];
		for(int i = 0; i<difference+1; i++){
			vector[i][0] = Matlab.round(lowerLimit + i*distance, 2);
		}
		return vector;
	}
	
	/*public static double average (double[][] vector){
		double average = 0;
		for(int i=0; i<vector.length; i++){
			average = average + vector[i][0];
		}
		average = average / vector.length;
		return average;
		
	}*/
	
	/**
	 * Clustering K-means. 
	 * 'Start' - Selecciona k observaciones de X en forma aleatoria.
	 * 'EmptyAction' - Acción que toma si un cluster pierde todos sus miembros,
	 * 'singleton' - Crea un nuevo cluster consitente en la observación más lejana al centroide
	 * @param X double[][] vector de entrada
	 * @param K int número de clusters
	 * @return centroids double[][] retorna los K cluster con sus centroides
	 */
	public static double[][] kmeans(double[][] X, int K) {
		
		final int CLUSTERS 		= K;
		double[][] data 		= X;
		int length 				= data.length;
		double[][] sums 		= new double[CLUSTERS][length];
		//double[][] centroids 	= {{0, 0}, {Math.random()*10+1, Math.random()*10+1}};  //1° actual - 2° anterior
		//double[][] centroids 	= {{0, 0}, {1, 9}};
		//double[][] centroids 	= {{0, 0}, {average(X)+0.5, average(X)-0.5}};
		double[][] centroids 	= {{0, 0}, {X[(int) Math.floor(Math.random()*59)][0], X[(int) Math.floor(Math.random()*59)][0]}};  //1° actual - 2° anterior
		int[] count 			= new int[CLUSTERS];
		int i, j, k;
		double minimum, difference;
		boolean converged = false;
		float div;

		do {
			for(i = 0; i < CLUSTERS; i++) {			//Cluster Anterior por Actual
				centroids[0][i] = centroids[1][i];
				count[i] = 0;
				centroids[1][i] = 0;
			}
			
			for(i = 0; i < length; i++) {
				
				sums[0][i] = 0;
				minimum = centroids[0][0] > data[i][0] ? centroids[0][0] - data[i][0] : data[i][0] - centroids[0][0];
				k = 0;
				
				for(j = 1; j < CLUSTERS; j++) {
					
					sums[j][i] = 0;
					difference = centroids[0][j] > data[i][0] ? centroids[0][j] - data[i][0] : data[i][0] - centroids[0][j];
					
					if(difference < minimum) {
						
						minimum = difference;
						k = j;
					}
				}
				
				sums[k][i] = data[i][0];
				count[k]++;
			}
			
			converged = true;
			
			for(i = 0; i < CLUSTERS; i++) {
				
				if(count[i] == 0 && i==0){
					centroids[1][i] = Math.abs(centroids[0][0]-Matlab.min(data)) >= Math.abs(centroids[0][0]-Matlab.max(data)) ? Matlab.min(data) : Matlab.max(data);
				} else if(count[i] == 0 && i==1){
					centroids[1][i] = Math.abs(centroids[0][0]-Matlab.min(data)) >= Math.abs(centroids[0][0]-Matlab.max(data)) ? Matlab.min(data) : Matlab.max(data);
				} else if(count[i] > 0) {
					for(j = 0; j < length; j++) {
						div = (float) (sums[i][j] / count[i]);
						centroids[1][i] += div;
					}
				}
				
				centroids[0][i] = Matlab.round(centroids[0][i], 4);
				centroids[1][i] = Matlab.round(centroids[1][i], 4);
				converged &= centroids[0][i] == centroids[1][i];
			}
		}while(!converged);

		/*for(i = 0; i < CLUSTERS; i++) {
			
			System.out.println(centroids[1][i]);
		}*/
			
		return centroids;
	}
	
	/**
	 * Diferencia a través de la primera fila de un vector fila
	 * @param m double[][] vector fila de entrada con dimensión N
	 * @return r vector fila de entrada con dimensión N-1
	 */
	public static double[][] diff(double[][] _vector) {
		//double aux = m[0].length;
		double[][] vector 	= new double [1][_vector[0].length-1];
		for(int i=0; i<(_vector[0].length-1); i++){
			vector[0][i] 	= _vector[0][i+1]-_vector[0][i];
		}
		return vector;
	}
	
	/**
	 * Esta función retorna el average de un vector fila
	 * @param m double[][] vector fila
	 * @return v double average de los datos
	 */
	public static double mean(double[][] vector) {
		float mean	= 0;
		float numberData 	= vector[0].length;
		for (int i = 0 ; i < numberData ; i++){
			mean = (float) (mean + vector[0][i]);
		}
		mean = mean/numberData;
		return (double) mean;
	}	
	
	/**
	 * Media
	 * @param vector double[][] vector de entrada
	 * @param k int ínfice de la fila o columna
	 * @param roc String indica fila o columna
	 * @return average double media de la fila o columna seleccionada
	 */
	public static double mean (double[][] vector, int k, String roc){
		float average = 0;
		if(roc.equals("columna")){
			for(int i=0; i<vector.length; i++){
				average = (float) (average + vector[i][k]);
			}
			average = average / vector.length;
		} else if(roc.equals("fila")){
			for(int i=0; i<vector[0].length; i++){
				average = (float) (average + vector[k][i]);
			}
			average = average / vector[0].length;
		}
		return average;
	}
	
	/**
	 * Varianza
	 * @param X double[][] vector fila 
	 * @param m double media de X
	 * @return varianza double retorna la varianza de los valores en X
	 */
	public static double variance(double[][] X, double m) {
		double variance = 0;
		//double aux1 	= 0;
		for(int i = 0; i<X[0].length; i++){
			//aux1 		= x[0][i]-m;
			//aux1 		= Math.pow(aux1, 2);
			//variance 	= variance + aux1;
			variance 	= variance + Math.pow(X[0][i]-m, 2);
		}
		//aux1 		= x[0].length-1;
		float aux2 	= (float) (variance/(X[0].length-1));
		variance 	= aux2;
		return variance;
	}
	
	/**
	 * Función que retorna la desviación estandar
	 * @param x double[][] matriz
	 * @param k int posición de la fila o columna
	 * @param roc String indica fila o columna
	 * @return desv double retorna la desviación estandar de la fila o columna seleccionada
	 */
	public static double std(double[][] matrix, int k, String roc) {
		//double media 	= mean(x,k,roc);
		//double var 	= varianza(x,media);
		double deviation = Math.sqrt(variance(matrix, mean(matrix, k, roc)));
		return deviation;
	}
	
	/**
	 * Función que una fila o columna de una matriz
	 * @param m double[][] matriz
	 * @param k int posición de la fila o columna
	 * @param roc String indica fila o columna
	 * @return r double[][] vector correspondiente a la posición k de la matriz m
	 */
	public static double[][] cut(double[][] matrix, int k, String roc) {
		if(roc.equals("columna")){
			double[][] r = new double[matrix.length][1];
			for(int i=0; i<matrix.length; i++) {
				r[i][0] = matrix[i][k];
			}
			return r;
		} 
		else if(roc.equals("fila")){
			double[][] r = new double[1][matrix[0].length];
			for(int j=0; j<matrix[0].length; j++) {
				r[0][j] = matrix[k][j];
			}
			return r;
		}
		
		return null;
	}
	
	/**
	 * Función que cuenta números enteros
	 * @param inicio int número entero de inicio
	 * @param fin int número entero final
	 * @return r int[][] vector columna de valores entre inicio y fin
	 */
	public static int[][] count(int initiation, int fin) {
		int[][] vector = new int[fin-initiation][1];
		for(int i = 0; i<fin-initiation;i++){
			vector[i][0] = initiation+i+1;
		}
		return vector;
	}
	
	/**
	 * Adhiere un valor a un vector columna
	 * @param vector double[][] vector columna de entrada
	 * @param position int indica la posición donde será añadido el nuevo valor
	 * @param value double valor a ser añadido
	 * @return vector double[][] retorna un vector columna con una fila añadida
	 */
	public static double[][] addValue(double[][] _vector, int position, double value) {
		double[][] vector 	= new double[_vector.length+1][1];
		for(int i = 0; position-1>i; i++){
			vector[i][0]	= _vector[i][0];
		}
		vector[position-1][0] = value;
		for(int i = position; _vector.length+1>i; i++){
			vector[i][0]	= _vector[i-1][0];
		}
		return vector;
	}
	
	/**
	 * Método que retorna los indices de un vector dada la cuenta del histograma
	 * @param vector double[][] vector fila de entrada (muestra)
	 * @param edges [][] vector columna de los límites el cual debe contener valores monotónicamente no decrecientes 
	 * @return m int[][] devuelve vector de indices
	 */
	public static int[][] histc(double[][] vector, double[][] edges) {
		int [][] vectorIndex = new int [vector[0].length][1];
		for(int k = 0; k<edges.length; k++){
			for(int i = 0; i<vector[0].length; i++){
				if(edges[k][0]<=vector[0][i] && vector[0][i]<edges[k+1][0]){
					vectorIndex[i][0] = k+1;
				}
			}
		}
		return vectorIndex;
	}
	
	/**
	 * Función que crea un vector fila con un rango entre [ini,fin]
	 * @param ini float valor inicial del rango
	 * @param cada float valor que indica el espacio entre cada dato
	 * @param fin int valor que indica el fin del rango
	 * @return m double[][] vector fila con el rango creado
	 */
	public static double[][] sample(float initiation, float each, int fin) {
		float lfinal = (fin-initiation+each)/each;
		double[][] m = new double [1][(int) lfinal];
		for(int k = 0; k<(int) lfinal; k++){
			m[0][k] = initiation+k*each;
		}
		return m;
	}

}