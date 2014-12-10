package soc;

import mathlib.Matlab;

/**
 * 
 * @author Alex J. Díaz Millán
 *
 */
public class MarkovChain {

	private double[][] predictionCurrent;
	
	private double lowCurrent;
	private double highCurrent;
    private double P12_out;
    private double P21_out;
	
    /**
     * Inicializa vector de predicción de corriente
     * @param predictionSize Int Largo de la predicción
     */
	public void initializePredictionCurrent(int predictionSize){
		predictionCurrent	= Matlab.zeros(1, predictionSize + 1);
	}
	
	/**
	 * Se establece el valor de la predicción de corriente
	 * @param _predictionCurrent Double[][] Predicción de corriente
	 */
	public void setPredictionCurrent(double[][] _predictionCurrent){
		predictionCurrent = _predictionCurrent;
	}
	
	/**
	 * Se establece el valor en la posición i de la predicción de corriente
	 * @param i Int Posición a establecer el valor
	 * @param value Double Valor a ser establecido
	 */
	public void setPredictionCurrent(int i, double value) {
		predictionCurrent[0][i]	= value;
	}
	
	/**
	 * Entrega el valor de la predicción de corriente
	 * @return Double[][] Predicción de corriente
	 */
	public double[][] getPredictionCurrent() {
		return predictionCurrent;
	}
	
	/**
	 * Método que realiza la caracterización del perfil de uso de la batería
	 * @param currentSensor_k Double[][] Corriente sensada hasta la iteración k
	 * @param k Int Iteración k
	 * @param predictionStart Int Comienzo de la predicción
	 */
	public void CharacterizationProfile(double[][] currentSensor_k, int k, int predictionStart) {
		// Se definen los intervalos y datos a considerar
		double[][] discretizedCurrent	= Matlab.zeros(predictionStart, 1);		// Vector corriente discretizado
		double forgettingFactor    		= 0.65;                             	// Ponderador para ecuación EWMA (valor empírico)

		int largo_datos 	= predictionStart;						
		int delta       	= 60;                               				// Tamaño de cada intervalo (en segundos)
		int intervals  		= (int) Math.floor((float) (largo_datos/delta));	// Número de intervalos (en segundos)
		int rest			= (int) (largo_datos - intervals*delta);			// Datos de corriente no considerados en CM (se toman los primeros datos)

		double[][] minimumInterval  = Matlab.zeros(intervals, 1);           	// Limites superior e inferior en cada intervalo (Centroides de cada Cluster)
		double[][] maximumInterval	= Matlab.zeros(intervals, 1);    
		double[][] upperLimit       = Matlab.ones(delta*intervals, 1);       	// Limites superior e inferior en cada instante de tiempo
		double[][] lowerLimit      	= Matlab.ones(delta*intervals, 1);
		double[][] finalUpperLimit	= Matlab.ones(delta*intervals, 1);
		double[][] finalLowerLimit 	= Matlab.ones(delta*intervals, 1);
		
		// Probabilidades de la Cadena de Markov (2 estados)
		double p12 = 0.45;                                         				// Probabilidades iniciales (en base a paper)
		double p22 = 0.45;
		double p11 = 1 - p12;
		double p21 = 1 - p22;

		// Probabilidades de transición discretas finales "por intervalo"
			// Probabilidad P22[k] es calculable con los restantes
		double[][] P11 = Matlab.zeros(intervals, 1);								
		double[][] P12 = Matlab.zeros(intervals, 1);                       		
		double[][] P21 = Matlab.zeros(intervals, 1);

		// Cantidad de transiciones de estado "por intervalo" P_ij[k]
		double[][] counterP11  = Matlab.zeros(intervals, 1);     		
		double[][] counterP12  = Matlab.zeros(intervals, 1);
		double[][] counterP21  = Matlab.zeros(intervals, 1);

		// Estado (BAJO = 1 o ALTO = 2) de corriente para cada instante de muestreo
		double[][] stateCM 	= Matlab.zeros(delta,intervals);					
		
		// Discretización en 2 estados en cada intervalo y transiciones
		int fin = rest;
		int initial;
		double[][] datos_inter, C;
		int indice;
		float aux;

		for (int i=0; i<intervals; i++){

		    initial  				= fin + 1;											// Inicio desde la primera muestra considerada en la CM
		    fin     				= fin + delta;                          			// La próxima muestra será DESDE la primera del intervalo siguiente
		    datos_inter 			= Matlab.cutVector(currentSensor_k, initial, fin);
	        C 						= Matlab.kmeans(datos_inter,2);       				// K-Means Clustering (K=2) por cada intervalo para determinar estados
	        minimumInterval[i][0] 	= Math.min(C[0][0],C[0][1]);						// Valor de consumo BAJO (estado 1)
	        maximumInterval[i][0] 	= Math.max(C[0][0],C[0][1]);						// Valor de consumo ALTO (estado 2)		    
		    
		    for(int j = (int) (initial-1); j<fin; j++){
		        indice = (j - delta*(i+1-1));
		        if ( (currentSensor_k[j][0] - minimumInterval[i][0]) < (maximumInterval[i][0] - currentSensor_k[j][0]) ){	
		            stateCM[indice][i] 			= 1;								// Se discretiza al valor I BAJO
		            discretizedCurrent[j][0]	= minimumInterval[i][0];
		        } else {
		            stateCM[indice][i] 			= 2;           						// Se discretiza al valor I ALTO
		            discretizedCurrent[j][0] 	= maximumInterval[i][0];
		        }
		    }

		 // Contador de transiciones entre estados por intervalo
		    for (int j = 0; j<delta-1; j++){ 
		        if (stateCM[j][i] == 1){
		         // Se mantiene consumo BAJO
		            if (stateCM[j+1][i] == 1){
		                counterP11[i][0] = counterP11[i][0] + 1; 					
		            } 
		         // Existe transición de consumo BAJO a ALTO
		            else { 
		                counterP12[i][0] = counterP12[i][0] + 1;					
		            }
		        } 
		        else {
				 // Existe transición de consumo ALTO a BAJO
		            if (stateCM[j+1][i] == 1) {   									
		                counterP21[i][0] = counterP21[i][0] + 1;
		            }
		        }
		    }

		 // Determinación de probabilidades de transición por intervalo
		    if ((counterP11[i][0] + counterP12[i][0]) == 0){		
		        if (i == 1){
		            P11[i][0] = p11;
		            P12[i][0] = p12;
		        } else {
		            P11[i][0] = P11[i-1][0];
		            P12[i][0] = P12[i-1][0];
		        }
		        
		        if (counterP21[i][0] == 0){ 									// Uso contínuo de corriente ALTA
		            P21[i][0] = 0;                          					// P22[i][0] = 1
		        } else {
		        	aux 		= (float) (counterP21[i][0]/(delta - 1));
		            P21[i][0] 	= aux;                     						// P21[i][0] = 1 / (N°de muestras por intervalo - 1)
		        }
		    } else {
		    	aux = (float) (counterP11[i][0]/(counterP11[i][0] + counterP12[i][0]));
		        P11[i][0] = aux;                   								// Transiciones desde corriente BAJA
		        aux = (float) (counterP12[i][0]/(counterP11[i][0] + counterP12[i][0]));
		        P12[i][0] = aux;
		        if ((counterP11[i][0] + counterP12[i][0]) == (delta-1)){		// Uso contínuo de corriente BAJA
		            if (i == 1){
		                P21[i][0] = p21;
		            } else {
		                P21[i][0] = P21[i-1][0];
		            }
		        } else {
		        	aux = (float) (counterP21[i][0]/(delta - (counterP11[i][0] + counterP12[i][0]) - 1));
		            P21[i][0] = aux; 											// Transiciones desde corriente ALTA
		        }
		    }
		    
		    // Valores MIN y MAX sin aplicar EWMA a la CM
		    for(int j=(int) (initial-1-rest); j<fin-rest; j++){
		    	upperLimit[j][0] = maximumInterval[i][0];                        			
		    	lowerLimit[j][0] = minimumInterval[i][0];
		    }
		    
		    // Ponderación EWMA (Exponentially Weighted Moving Average)
		    if ( i!=0 ){
		        maximumInterval[i][0] = forgettingFactor*maximumInterval[i][0] + (1-forgettingFactor)*maximumInterval[i-1][0];
		        minimumInterval[i][0] = forgettingFactor*minimumInterval[i][0] + (1-forgettingFactor)*minimumInterval[i-1][0];

		        P11[i][0] = forgettingFactor*P11[i][0] + (1-forgettingFactor)*P11[i-1][0];
		        P12[i][0] = forgettingFactor*P12[i][0] + (1-forgettingFactor)*P12[i-1][0];
		        P21[i][0] = forgettingFactor*P21[i][0] + (1-forgettingFactor)*P21[i-1][0];
			} else {
		        P11[i][0] = forgettingFactor*P11[i][0] + (1-forgettingFactor)*p11;
		        P12[i][0] = forgettingFactor*P12[i][0] + (1-forgettingFactor)*p12;
		        P21[i][0] = forgettingFactor*P21[i][0] + (1-forgettingFactor)*p21;
		    }
		    
		    // Valores MIN y MAX aplicando EWMA (Ajuste de importancia a intervalos finales)
		    for(int j = (int) (initial-1-rest); j<fin-rest; j++){
		    	finalUpperLimit[j][0] = maximumInterval[i][0];
		    	finalLowerLimit[j][0] = minimumInterval[i][0];
		    }
		    
		}

		// Variables de salida (MAX y MIN ponderados por EWMA)
		lowCurrent	= minimumInterval[intervals-1][0];
		highCurrent	= maximumInterval[intervals-1][0];
		// Variables de salida (Prob. de transición entre estados con EWMA)
		P12_out  	= P12[intervals-1][0];                     	
		P21_out  	= P21[intervals-1][0];
		
	}

	/**
	 * Método que realiza la Cadena de Markov, establece la predicción de la corriente dado las probabilidades de transición y los estados calculados en la caracterización del perfil de uso de la batería 
	 * @param predictionSize Int Tamaño de la predicción
	 */
	public void RealizationMC(int predictionSize) {
		
		// Genera un número aleatorio entre 0 y 1 (Uniforme)
		double random   			= Math.random();		
		
		// Asignación inicial aleatoria del perfil
		if (random > P12_out){			// Probabilidad de mantener consumo bajo en el perfil P11					
			predictionCurrent[0][0]	= lowCurrent;
		} else {
			predictionCurrent[0][0] = highCurrent;
		}

		// Generación aleatoria del perfil de uso futuro
		for (int k = 0; k < predictionSize; k++){
			random   = Math.random();

		    if (predictionCurrent[0][k] == lowCurrent){
		    	// Probabilidad de mantener consumo bajo en la iteracion siguiente
		        if (random > P12_out){           						
		        	predictionCurrent[0][k+1] = predictionCurrent[0][k];	// P11(k) > P12(k)
		        } else{                                                                 
		        	predictionCurrent[0][k+1] = highCurrent;            	// P12(k) > P11(k)
		        }
		    } else {
		    	// Probabilidad de cambiar a un consumo bajo en la iteracion siguiente
		        if (random > (1 - P21_out)){
		        	predictionCurrent[0][k+1] = lowCurrent;					// P21(k) > P22(k) 
		        } else {
		        	predictionCurrent[0][k+1] = predictionCurrent[0][k];    // P22(k) > P21(k)
		        }
		    }
		}
		
	}
	
}