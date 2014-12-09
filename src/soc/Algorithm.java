package soc;

import mathlib.Matlab;
import data.DataManagement;
import data.ParameterizationBattery;
import data.Route;

public class Algorithm {
	
	private static int predictionStart, predictionSize;
	private static int realizationsMC		= 25;               // N° de realizaciones de Cadenas de Marvok o de prediccion en pronóstico

	public static void main(String[] args) throws Exception {
		
		// Inicio Tiempo de ejecución
		long startTime, finalTime;
		//startTime = System.currentTimeMillis();
		
		// Carga de valores Parametrización (OFF-LINE)
		ParameterizationBattery dataBattery = new ParameterizationBattery();
		double[][] voltageOcSen = dataBattery.getVoltageOcSen();
		double[][] socSen		= dataBattery.getSocSen();
		
		// Carga de mediciones sensor de voltaje y sensor de corriente (OFF-LINE)
		Route route = new Route();
		double[][] currentSensor = route.getCurrentSensor();
		double[][] voltageSensor = route.getVoltageSensor();
		
		startTime = System.currentTimeMillis();
		
		//Condición inicial de ejecución algoritmo (Voltaje inicial mayor al Vcut-off + 0.5[V])
		double cutoff = 30 + 0.5;				// Vcut-off programado en PCM 
		if (voltageSensor[0][0] <= cutoff){ 
		   System.out.print("Voltaje insuficiente para alimentar motor");
		}
		
		// Números de largo_iteraciones para filtrado
		int lengthPF   			= voltageSensor.length;		// Tiempo en segundos de estimación (largo_iteraciones del FP)
		predictionStart  		= 900;              		// Tiempo inicial de predicción medido en segundos
		int predictionHorizon   = lengthPF;          		// Horizonte de predicción medido en segundos
		//double criticalEnergy 		= 1065600;					// Energía nominal de la batería [V*A*seg](criticalEnergy_exp = 993135)
		double deltaT     		= 1;						// Tiempo de muestreo [seg]
		
		//double[][] internalResistancePrediction   	= Matlab.zeros(1, (predictionHorizon-predictionStart));	
		
		double characterizationProfile      		= 0;            											// Caracterización del perfil futuro ( =1 Ipromedio; ~=1 CM )
		
		// Modelo de transición de estados (X1=Resistencia Interna, X2=SOC) 
		//--------------------------MODELO-------------------------------------------------------------------------------------------------------
		//   x1(k+1) = x1(k) + w1(k)
		//   x2(k+1) = x2(k) + (v(k)*i(k)*delta_T)/criticalEnergy + w2(k)
		//   v(k) = vL + (v0-vL)*exp(gamma*(x2(k)-1)) + alfa*vL*(x2(k)-1) + (1-alfa)*vL*(exp(-beta)-exp(-beta*sqrt(x2(k)))) - i(k)*x1(k) + eta(k)
		//---------------------------------------------------------------------------------------------------------------------------------------
		Model model 			= new Model();
		model.initializeStates(lengthPF);							// PDF inicial estados
		model.initializeEstimatedVoltage(currentSensor[1][0], lengthPF);	// Se genera el vector de observación inicial
		
		// Filtro de Partículas
		ParticleFilter filter 	= new ParticleFilter();
		filter.initializeParticles(model);							// Inicializan las partículas
		filter.initializeWeights();									// Inicializan los pesos de las partículas
		filter.setReliability(0.95);								// Probabilidad considerada para formar el intervalo de confianza
		
		model.initializeEstimatedVoltageK(filter.getNumberParticles());
		
		// Cálculo del SOC real (OFF-LINE)
		double[][] realSOC	= CalculateRealSOC(lengthPF, voltageSensor, currentSensor, deltaT);	// Cálculo del SOC "real"
		
		// Definición de Ground Truth EOD asociado a Vcutoff
		int eodGroundTruth	= Matlab.find(voltageOcSen, cutoff);	// Instante de tiempo [seg] a pronosticar                                                      
		model.setSocCutOff(socSen[0][eodGroundTruth]); 		         	// SOC a pronosticar
																		// Cálculo HZ_SOC para DistributionPrediction.m
 		
		//----- Etapa de estimación y pronóstico (FP)--------------------------------------------------------------
		double[][] voltage 	= Matlab.transpose(Matlab.sample(voltageSensor, 1, 1, lengthPF));
		double[][] states, particles;
		
		for (int k = 1; k < lengthPF; k++){
			
			filter.setK(k);
		    states = model.getStates();
		    
		    // Ruido de proceso 2 adaptativo
		    if (k < 30) {
		    	model.setSigmaW2(0.016);
		    } else if ( Math.abs(states[1][k-1] - realSOC[k-1][0]) > Math.abs(states[1][k-2] - realSOC[k-2][0]) ) {
		    	model.setSigmaW2(0.016*1.5);
		    } else {
		    	model.setSigmaW2(0.016/5);	// REVISAR LA DIVISION
		    }
		   
		    // Proceso de Estimación
		    EstimationPF(model, filter, voltage[0][k], voltage[0][k-1], currentSensor[k][0]);
		   		
			// Voltaje estimado en iteración k según modelo
			model.setEstimatedVoltage(filter.getK(), currentSensor[k][0]);
			
			if (k == predictionStart) {
				
				PredictionPF(filter, model, Matlab.sample(currentSensor, 1, 1, k), predictionHorizon, characterizationProfile);
						 
				model.setPredictionVoltage(0, model.getEstimatedVoltage(k));
				particles = filter.getParticles();
				
				for (int i = 0; i<filter.getNumberParticles(); i++){
					model.setPredictionSoc(0, model.getPredictionSoc(0) + particles[1][i]*filter.getWeights(i));
				}
				
				for (int i = 0; i<predictionSize; i++){
					model.setPredictionSoc(i, model.getPredictionSoc(i) * 100);
				}
			}
		}
		
		double [][] estimatedSOC		= new double [1][lengthPF];
		double [][] internalResistance	= new double [1][lengthPF];
		states = model.getStates(); //Revisar... no es necesario pedirlo.
		for(int i=0; i<lengthPF; i++){
			internalResistance[0][i] 	= Math.abs(states[0][i]);	// Estado 1 estimado [Ohms] (valores esperados)
			estimatedSOC[0][i] 			= states[1][i]*100;			// Estado 2 estimado [%] (valores esperados)
		}
		
		for (int i = 0; i<predictionSize; i++) {
			model.setPredictionInternalResistance(i, internalResistance[0][predictionStart+i]);
		}
		
		//  Cálculo del error de estimación (RMSE)
		double[][] sampleVoltageSensor 	= Matlab.sample(voltageSensor, 1, 1, lengthPF);
		double[][] sampleRealSOC 		= Matlab.sample(realSOC, 1, 1, lengthPF);
		double[][] voltageError			= new double [lengthPF][1];
		double[][] socError 			= new double [lengthPF][1];
		
		for(int i = 0; i<lengthPF; i++){
			voltageError[i][0]	= Math.pow((sampleVoltageSensor[i][0] - model.getEstimatedVoltage(i)), 2);
			socError[i][0] 		= Math.pow((sampleRealSOC[i][0] - estimatedSOC[0][i]), 2); 
		}
		
		double averageVoltageError   	= Math.sqrt(Matlab.mean(voltageError, 0, "columna")); 
		double averageSocError			= Math.sqrt(Matlab.mean(socError, 0, "columna"));
		
		// Guardar Resultados obtenidos
		saveData(sampleRealSOC, estimatedSOC, model.getEstimatedVoltage(), model.getPredictionVoltage(), model.getPredictionSoc(), averageVoltageError, averageSocError);
		
		// Fin tiempo de ejecución
		finalTime = System.currentTimeMillis();
		System.out.println("La tarea ha tomado "+ ( finalTime - startTime ) +" milisegundos");

	}
	
	/**
	 * Algoritmo de Estimación basado en Filtro de Partículas
	 */	
	private static void EstimationPF(Model model, ParticleFilter filter, double voltageObserved, double voltage_1, double current) {
		
		// Estimación del estado por cada partícula
		filter.setParticles( model.calculationParticles(filter, voltage_1, current) );
		
		// Actualizacion de pesos (cuando se utiliza la PDF a priori)
		model.setVoltageEstimated_k(filter, current);
		filter.updateWeights(voltageObserved, model);
		
		// Normalización del vector de pesos
		filter.WeightsNormalization();
		
		// Calculo del índice de eficiencia
		double Neff = 1/Matlab.sum2(filter.getWeights(), 0, "columna");		//sum2 -> Suma de w^2 [REVISAR DIVISIÓN -> funcionó]
	
		// Condición de degeneración de particulas (Remuestreo o Resampling)
		double umbral 		   	= 0.85;							// Valor límite de entrada a Resampling
		double umbral_global   	= umbral * filter.getNumberParticles();	
		if (Neff < umbral_global){
			filter.Resampling();
		}
	
		// Se calcula la estimación del estado
		model.StatesEstimation(filter);
	}

	private static void PredictionPF(ParticleFilter filter, Model model, double[][] currentSensor_k, int predictionHorizon, double characterizationProfile) {
		
		/*
		 * Algoritmo de Predicción basado en Filtro de Partículas
		 */
		
		// Se promedian todas las predicciones para un pronóstico (Cada predicción conlleva 1 realización de CM)
		predictionSize  = predictionHorizon - filter.getK();	// Tamaño de predicción = Horizonte de predicción - Inicio de predicción
		
		// Matrices resultantes para cada realizacion
  		double[][] realizationTof   			= Matlab.zeros(realizationsMC, 1);
  		double[][] realizationMinIC				= Matlab.zeros(realizationsMC, 1);
  		double[][] realizationMaxIC   			= Matlab.zeros(realizationsMC, 1);
  		double[][] realizationPredictionSoc 	= Matlab.zeros(realizationsMC, predictionSize);
  		double[][] realizationPredictionVoltage	= Matlab.zeros(realizationsMC, predictionSize);
  		double[][] realizationPdfEod  			= Matlab.zeros(realizationsMC, predictionSize);
  		
  		// Número de realizaciones que no cruzan el umbral de falla (Todas las partículas)
  		double[][] realizationsNotCross 		= Matlab.zeros(realizationsMC, 1);
  		
  		// Vectores resultantes promedio
  		model.initializePrediction(predictionSize);
  		
  		// Cadena de Markov
  		MarkovChain mc	= new MarkovChain();
  		mc.initializePredictionCurrent(predictionSize);
		
		// Selección de enfoque de perfil futuro
	    if (characterizationProfile == 1){
	    	// Corriente predicha constante promediada
	        for(int i =0; i<predictionSize; i++) {
	        	mc.setPredictionCurrent(i, Matlab.mean(Matlab.sample(currentSensor_k, 1, 1, filter.getK()), 0, "fila") );
	    	}
	        realizationsMC	= 1;
	    } else {
	    	// Valores de los estados y matriz de transición (probabilidades) para la CM
	    	mc.CharacterizationProfile(currentSensor_k, filter.getK(), predictionStart);
	    }
  		
  		// Calculo de la distribucion de Epanechnikov para Prediccion_Reg.m
  		filter.CalculationEpanechnikov();
  		
		// Realizaciones de la Cadena de Markov
		for (int j=0; j<realizationsMC; j++){
			
			// Genera el perfil de uso futuro de la batería (una realización)
		    if (characterizationProfile != 1) {
		        mc.RealizationMC(predictionSize);
		    }
		    
		    // Se obtienen las predicciones, la PDF, el tof y los límites del intervalo de confianza
		    PredictionDistribution(filter, model, mc.getPredictionCurrent(), predictionSize );
		    
		    // Resultados SOC (Distribuciones FDP, intervalos, tof, vector de estado predicho)
		    for (int i=0; i<realizationPredictionVoltage[0].length; i++){
			    realizationPredictionVoltage[j][i]  = model.getPredictionVoltage(i);
			    realizationPredictionSoc[j][i]		= model.getPredictionSoc(i);                             
			    realizationPdfEod[j][i]  			= model.getPdfEOD(i);
		    }
		    
		    realizationTof[j][0]	= filter.getTof();
		    realizationMinIC[j][0]	= filter.getMin();
		    realizationMaxIC[j][0]  = filter.getMax();
		    
		    // Cruza o no cruza el umbral de falla
		    if (realizationTof[j][0] == 0){
		        realizationsNotCross[j][0] = 1;	 // 1 si no cruza?	
		    }
		    
		}
		
		// Se obtienen cuántas realizaciones no cruzaron el umbral de falla
		double realizationsNotCross_total = Matlab.sum(realizationsNotCross, 0, "columna");

		// Promedio de los resultados de predicción
		double total_SOC;
		if (realizationsMC == realizationsNotCross_total) {	// Cuando NO se entra a la Hazard Zona, no se consideran los resultados de esa realización
		    total_SOC = 1;
		} else {
		    total_SOC = realizationsMC - realizationsNotCross_total;
		}
		

		// Promedio de los resultados de prediccion CASO SI ENTRA EN HZ
		if (realizationsMC != 1){	// Al no cruzar el umbral de falla, se promedian los resultados de las realizaciones que si lograron entrar a la HZ
			// Promedios
			filter.setTof( (float) (Matlab.sum(realizationTof,0,"columna")/total_SOC) );
		    filter.setMin( (float) (Matlab.sum(realizationMinIC,0,"columna")/total_SOC) );
		    filter.setMax( (float) (Matlab.sum(realizationMaxIC,0,"columna")/total_SOC) );
		    
		    if (realizationsNotCross_total != 0){	// Se suman las predicciones y la PDF de las realizaciones "válidas"
		    
		        for (int j=0; j<realizationsMC; j++) {
		        	// j-ésima realización se suma a cada variable
		            if(realizationsNotCross[j][0] != 1) { 
		            	for (int k=0; k<predictionSize; k++){
			                model.setPredictionSoc(j, model.getPredictionSoc(j) + realizationPredictionSoc[j][k] );
			                model.setPredictionVoltage(j, model.getPredictionVoltage(j) + realizationPredictionVoltage[j][k] );
							model.setPdfEOD(j, model.getPdfEOD(j) + realizationPdfEod[j][k] );
		            	}
		            }
		        }
		        
		        // Se promedian las predicciones y la PDF
		        for(int i=0; i<predictionSize; i++) {
		        	model.setPredictionSoc(i, (float) (model.getPredictionSoc(i)/total_SOC));
		        	model.setPredictionVoltage(i, (float) (model.getPredictionVoltage(i)/total_SOC) );
				    model.setPdfEOD(i, (float) (model.getPdfEOD(i)/total_SOC) );
		        }

		    } 
		    else // Caso que cruza el umbral de falla (Resultado final)
		    {							
		    	//float aux = 0;
		    	for(int i=0; i<predictionSize; i++) {
		    		model.setPredictionSoc(i, (float) (Matlab.sum(realizationPredictionSoc,i,"columna")/realizationsMC) );
			    	model.setPredictionVoltage(i, (float) (Matlab.sum(realizationPredictionVoltage,i,"columna")/realizationsMC) );
			        model.setPdfEOD(i, (float) (Matlab.sum(realizationPdfEod,i,"columna")/realizationsMC) );
		    	}
		    }
		    
		} else {	// Caso 1 realización de predicción (1 realización de CM)
		    model.setPredictionVoltage(realizationPredictionVoltage);
		    model.setPredictionSoc(realizationPredictionSoc);
		    model.setPdfEOD( realizationPdfEod );
		    filter.setTof( (float) realizationTof[0][0] );
		    filter.setMin( (float) realizationMinIC[0][0] );
		    filter.setMax( (float) realizationMaxIC[0][0] );
		}
		
	}

	private static void PredictionDistribution(ParticleFilter filter, Model model, double[][] predictionCurrent, int predictionSize) {
		
		int numberParticles 			= filter.getNumberParticles();
		
		// Predicción (requiere particulas de cada estado)
		model.PredictionModel(filter, predictionCurrent, predictionSize);
		double[][] predictionVoltage_N	= model.getVoltage_N();
		double[][] predictionSoc_N		= model.getPredictionSoc_N();
		
		// Restablecer la predicción del voltaje y del soc
		model.clearPredictionVoltage(predictionSize);
		model.clearPredictionSOC(predictionSize);
		
		for(int j=0; j<predictionSize; j++){
			for(int i=0; i<numberParticles; i++){
				// Valor esperado del voltaje (voltaje predicho) con pesos constantes
				model.setPredictionVoltage(j, model.getPredictionVoltage(j) + filter.getWeights(i)*predictionVoltage_N[i][j]);
				// SOC predicho con pesos constantes
				model.setPredictionSoc(j, model.getPredictionSoc(j)+filter.getWeights(i)*predictionSoc_N[i][j]);
			}
		}
		
		// Generación de distribución acumulada
		double[][] cumulativeDistributionSoc 	= Matlab.zeros(1, predictionSize);
		double socDistribution;
		float socPercentage;
		
		// Se genera la distribucion acumulada por partícula
		for (int j=0; j<numberParticles; j++){
			// Se itera en toda la ventana de pronostico
		    for (int i=0; i<predictionSize; i++){
		    	
		        if (predictionSoc_N[j][i] <= model.getUpperLimitSoc() && predictionSoc_N[j][i] >= model.getLowerLimitSoc()){
		            socDistribution	= Math.abs(predictionSoc_N[j][i] - model.getUpperLimitSoc());		// Distancia o diferencia entre el valor de y(j,i) y el limite superior
		            socPercentage	= (float) (socDistribution / model.getWidthHazardZone());
		            cumulativeDistributionSoc[0][i] = cumulativeDistributionSoc[0][i] + filter.getWeights(j) * socPercentage;
		        } else {                                                       
		            if(predictionSoc_N[j][i] < model.getLowerLimitSoc()){                               	// Si pasa el limite se suma completamente
		                cumulativeDistributionSoc[0][i] = filter.getWeights(j) + cumulativeDistributionSoc[0][i];
		            }
		        }
		        
		    }
		}
		
		// Generación PDF (densidad) de probabilidad del EOD
		double[][] socDensity_aux = PdfEod(cumulativeDistributionSoc);	// Vector PDF del EOD
		
		//double[][] socDensity 	= new double[1][predictionSize];
		
		for(int i=0; i<predictionSize; i++){
			if(predictionSize-1-i == 0){	// Math.abs(predictionSize-1-i)<.0000001
				model.setPdfEOD(predictionSize-1, 0);
			} else {
				model.setPdfEOD(i, socDensity_aux[0][i]);
			}			
		}
		
		if (Matlab.sum(model.getPdfEOD(), 0, "fila") != 0){
			for(int i=0; i<predictionSize; i++){
				model.setPdfEOD(i, (float) (model.getPdfEOD(i) / Matlab.sum(model.getPdfEOD(), 0, "fila")));	// Normalizacion de la PDF EOD
			}
		}
		
		int[][] socIndex 	= Matlab.count(predictionStart, (int) (predictionStart + predictionSize));
		double ToF_SOC = 0;		
		for(int i=0; i<predictionSize; i++) {
			ToF_SOC			= ToF_SOC + socIndex[i][0]* model.getPdfEOD(i);
		}
		
		ConfidenceInterval(filter, model, model.getPdfEOD(), ToF_SOC, predictionStart);
				
	}

	private static double[][] PdfEod(double[][] cumulativeDistributionSoc) {
		// Condiciones de CDF (Función Densidad Acumulada) 
		double[][] pdf   = Matlab.diff(cumulativeDistributionSoc);
		int length = pdf[0].length;
		int aux   = 0;  												// Variable que encuentra el instante de tiempo donde pdf es positivo
		double m, n;

		if (Matlab.sum(pdf,0,"fila") != 0){								// CDF del SOC no alcanza el umbral de la HZ (cdf = 0)
		    
		    if (pdf[0][length-1] < 0){ 									// Último término negativo
		    	pdf[0][length-1] = 0;
		    }

		    if (pdf[0][0] < 0){ 										// Primer término negativo  
		        for (int j=1; j<length; j++){
	                if (pdf[0][j] >= 0){								// Encuentro la posición del siguiente valor positivo
	                    aux = j;
	                    break;
	                }
		        }

		        m = (pdf[0][aux]-0)/(aux-0);							// Ecuación lineal de valores positivos (Mapeo)
		        n = pdf[0][aux] - m*aux;
		           
		        for (int j=0; j<aux; j++){								// Hay que reemplazar los valores negativos
		             pdf[0][j] = m*j + n;
		        }
		    }  

		    for (int i=0; i<length-1; i++){								// Condición fuera de borde (valor negativo entremedio)
		        if (pdf[0][i+1] < 0){									// Encuentro la posición del valor positivo(o cero) anterior al valor negativo
		            for (int j=i+1; j<length; j++){
		                if (pdf[0][j] >= 0){ 							// Encuentro la posición del siguiente valor positivo
		                    aux = j;
		                    break;
		                }
		            }
		            
		            m = (float) ((pdf[0][aux] - pdf[0][i])/(aux - i));	// Hay que encontrar la ecuación de la recta entre los dos valores positivos
		            n = pdf[0][i] - m*i;
		           
		            for (int j=i+1; j<aux; j++){                    	// Hay que reemplazar los valores negativos
		                pdf[0][j] = m*j + n;
		            }
		        }
		    }
		}
		
		return pdf;
	}

	private static void ConfidenceInterval(ParticleFilter filter, Model model, double[][] pdf, double tofSOC, double predictionStart){
		// Inicialización de variables
		tofSOC       	= tofSOC - predictionStart;
		int min, max;
		double initialProbability, lowerLimit, upperLimit;
	
		// Cálculo del IC para cada instante (IC variable)
		if (tofSOC < 0) {
		    min = 0;
		    max = 0;
		} else {
		    tofSOC    			= Math.round(tofSOC);		// Número entero
		    initialProbability	= pdf[0][(int) tofSOC-1];	// Probabilidad de estar inicialmente en el ToF 
		    max        			= (int) tofSOC;            	// Ambos límites comienzan en la posición del ToF
		    min        			= (int) tofSOC;
		    lowerLimit 			= 1;                        // Instantes de tiempo donde comienza y termina pdf
		    upperLimit 			= pdf[0].length;
	
		    // Búsqueda desde el inicio
			for (int i=0; i<pdf[0].length; i++) {
			    if (pdf[0][i] != 0){
			    	lowerLimit = i;
			        break;
			    }
			}
			
			// Búsqueda desde el final
			for (int i = pdf[0].length-1;i>0; i--){	
			    if (pdf[0][i] != 0){
			    	upperLimit = i;
			        break;
			    }
			}
	
			// Determinación del IC
			for (int i = 0; i<pdf[0].length; i++) {	
			    if (initialProbability < filter.getReliability()){      
			        if(min > lowerLimit && max < upperLimit) {   
			            min      = min - 1;
			            max      = max + 1;         
			            initialProbability = initialProbability + pdf[0][min] + pdf[0][max];        
			        } else {
			            if(min > lowerLimit && max >= upperLimit) {    
			                min      = min - 1;
			                initialProbability = initialProbability + pdf[0][min];
			            } else {
			                if (min <= lowerLimit && max < upperLimit){    
			                    max      = max + 1;
			                    initialProbability = initialProbability + pdf[0][max];
			                } else {
			                    break;
			                }
			            }
			        }
			    } else {      
			        break;
			    }
			}
	
			min = (int) (min + predictionStart);
			max = (int) (max + predictionStart);
	
		}
		filter.setMin(min);
		filter.setMax(max);
	}




	
	private static double[][] CalculateRealSOC(int lengthPF, double[][] voltageSensor, double[][] currentSensor, double deltaT) {
		double realInitialEnergy 		= 0;
		double[][] realSOC				= Matlab.zeros(lengthPF, 1);
		
		for(int i=0; i<voltageSensor.length; i++){
			realInitialEnergy			= realInitialEnergy+voltageSensor[i][0]*currentSensor[i][0]*deltaT;	// Energía inicial en [V*A*seg]	
		}
		
		double criticalEnergyReal 		= realInitialEnergy;                   				// Energía crítica real según set de datos [V*A*seg]
		float initialEnergyNormalized	= (float) (realInitialEnergy/criticalEnergyReal); 	// Energía inicial normalizada [V*A*seg]

		realSOC[0][0] 					= initialEnergyNormalized;                   		// Primer valor realSOC igual o menor a 1
		
		float DTCER						= (float) (deltaT / criticalEnergyReal);
		
		// Evolución del SOC real según modelo (extracción de energía)
		for (int i = 0;i<lengthPF-1;i++){
		    realSOC[i+1][0] = realSOC[i][0]-voltageSensor[i][0]*currentSensor[i][0]*DTCER;	
		}
		for (int i = 0;i<lengthPF;i++){
			realSOC[i][0]	= realSOC[i][0]*100;											// SOC real porcentual
		}
		
		return realSOC;
	}
	
	private static void saveData(double[][] sampleRealSoc, double[][] estimatedSOC, double[][] estimatedVoltage, double[][] predictionVoltage, double[][] predictionSoc, double averageVoltageError, double averageSocError) {
		DataManagement.grabarDatosExel(sampleRealSoc, "sampleRealSoc.xls");								// SOC real
		DataManagement.grabarDatosExel(Matlab.transpose(estimatedSOC), "estimatedSOC.xls");				// SOC estimado
		DataManagement.grabarDatosExel(Matlab.transpose(estimatedVoltage), "estimatedVoltage.xls");		// Voltaje estimado
		DataManagement.grabarDatosExel(Matlab.transpose(predictionVoltage), "predictionVoltage.xls");	// Voltaje predicho
		DataManagement.grabarDatosExel(Matlab.transpose(predictionSoc), "predictionSOC.xls");			// SOC predicho
		//DataManagement.grabarDatosExel(averageVoltageError, "averageVoltageError.xls");					// Promedio del Error de Voltaje
		//DataManagement.grabarDatosExel(averageSocError, "averageSocError.xls");							// Promedio del Error del SOC

	}

}