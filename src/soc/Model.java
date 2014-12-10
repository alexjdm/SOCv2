package soc;

import mathlib.Matlab;

/**
 * 
 * @author Alex J. D�az Mill�n
 * <p>Modelo de transici�n de estados (X1=Resistencia Interna, X2=SOC)</p>
 * <p>x1(k+1) = x1(k) + w1(k)</p>
 * <p>x2(k+1) = x2(k) + (v(k)*i(k)*delta_T)/energiaCritica + w2(k)</p>
 * <p>v(k) = vL + (v0-vL)*exp(gamma*(x2(k)-1)) + alfa*vL*(x2(k)-1) + (1-alfa)*vL*(exp(-beta)-exp(-beta*sqrt(x2(k)))) - i(k)*x1(k) + eta(k)</p>
 */
public class Model {
	
	// Modelo de transici�n de estados (X1=Resistencia Interna, X2=SOC) 
	//--------------------------MODELO-------------------------------------------------------------------------------------------------------
	//   x1(k+1) = x1(k) + w1(k)
	//   x2(k+1) = x2(k) + (v(k)*i(k)*delta_T)/energiaCritica + w2(k)
	//   v(k) = vL + (v0-vL)*exp(gamma*(x2(k)-1)) + alfa*vL*(x2(k)-1) + (1-alfa)*vL*(exp(-beta)-exp(-beta*sqrt(x2(k)))) - i(k)*x1(k) + eta(k)
	//---------------------------------------------------------------------------------------------------------------------------------------
	
	private double[][] states;
	
	// Vectores de estimaci�n
	private double[][] voltageEstimated;
	private double[][] voltageEstimated_k;
	
	// Vectores de predicci�n
	private double[][] predictionSoc;
	private double[][] predictionSoc_N;
	private double[][] predictionVoltage;	// Vector en donde se guardan los valores de la predicci�n del voltaje
	private double[][] predictionVoltage_N;
	private double[][] predictionInternalResistance;
	private double[][] pdfEOD;
	
	private int numberStates		= 2;	// N�mero de Estados
	private int numberObservations	= 1;	// N�mero de Observaciones
	
	// Ruido de proceso y observaci�n (PDF de observacion)
	private double muW1    	= 0;				
	private double sigmaW1 	= 0.005;		// Desviaci�n est�ndar ruido X1
	private double muW2    	= 0;
	private double sigmaW2 	= 0.016; 		// ESTE VALOR SERA VARIABLE CON UN CONDICIONAL "IF" SEGUN RMSE DEL SOC (SI ES ALTO SE AUMENTA RUIDO Y VICEVERSA)  "valor empirico 0.016"

	private double sigmaV  	= 0.3;			// Desviaci�n est�ndar ruido observacion
	
	private double deltaT	= 1;
	private double E_crit 	= 1065600;
	float deltaT_Ecrit 		= (float) (deltaT/E_crit);
	
	// Valores iniciales de los estados
	private double states1_0	= 0.2408;
	private double states2_0	= 1.0;
	
	// Par�metros del modelo
	private double vL    	= 39.2;			// Par�metro asociado a zona lineal
	private double v0    	= 41.49;        // Par�metro asociado a voltaje inicial
	private double alpha 	= 0.14;         // Par�metro asociado a zona lineal    
	private double beta  	= 9.29;         // Par�metro asociado a decaimiento codo final
	private double gamma 	= 6.69;         // Par�metro asociado al codo inicial y zona lineal
	private double Nu    	= 0;
	
	// Definici�n de Ground Truth EOD asociado a VcutOff
	private double socCutOff;				// SOC a pronosticar
	private double delta 	= 0.005;		
	private double upperLimitSoc, lowerLimitSoc, widthHazardZone;
	
	/**
	 * Inicializaci�n de estados. Se genera el vector de estados inicial
	 * @param lengthPF Int Largo del Filtro de Part�culas
	 */
	public void initializeStates(int lengthPF) {
		states	 		= Matlab.zeros(numberStates, lengthPF);
		states[0][0] 	= states1_0;
		states[1][0]	= states2_0;
	}
	
	/**
	 * Inicializaci�n de la Estimaci�n del Voltaje
	 * @param current Double Corriente
	 * @param lengthPF Int Largo del Filtro de Part�culas
	 */
	public void initializeEstimatedVoltage(double current, int lengthPF) {
		voltageEstimated		= Matlab.zeros(numberObservations, lengthPF);
		voltageEstimated[0][0] 	= vL+(v0-vL)*Math.exp(gamma*(states2_0-1))+alpha*vL*(states2_0-1)+(1-alpha)*vL*(Math.exp(-beta)-Math.exp(-beta*Math.sqrt(states2_0)))-states1_0*current+Matlab.normrnd(0,sigmaV);
	}
	
	/**
	 * Inicializaci�n del Voltaje Estimado por cada part�cula
	 * @param numberParticle Int N�mero de part�culas
	 */
	public void initializeEstimatedVoltageK(int numberParticle){
		voltageEstimated_k = new double[1][numberParticle];
	}
	
	/**
	 * Se establecen los estados estimados
	 * @param states_k Double[][] Estados estimados en la iteraci�n k
	 */
	public void setStates(double[][] states_k)
	{
		states = states_k;
	}
	
	/**
	 * Se establecen los par�metros asociados al SOC l�mite de corte
	 * @param _socCutOff Double SOC de corte establecido
	 */
	public void setSocCutOff(double _socCutOff){
		socCutOff 		= _socCutOff;
		upperLimitSoc 	= socCutOff + delta;
		lowerLimitSoc 	= socCutOff - delta;
		widthHazardZone = upperLimitSoc - lowerLimitSoc;
	}
	
	/**
	 * Retorna el l�mite inferior del SOC de corte
	 * @return Double L�mite inferior del SOC de corte
	 */
	public double getLowerLimitSoc() {
		return lowerLimitSoc;
	}
	
	/**
	 * Retorna el l�mite superior del SOC de corte
	 * @return Double L�mite superior del SOC de corte
	 */
	public double getUpperLimitSoc() {
		return upperLimitSoc;
	}
	
	/**
	 * Retorna el ancho de la zona de riesgo asociada al SOC de corte
	 * @return Double Ancho de la zona de riesgo asociada al SOC de corte
	 */
	public double getWidthHazardZone() {
		return widthHazardZone;
	}
	
	/**
	 * Retorna los estados estimados
	 * @return Double[][] Estados estimados
	 */
	public double[][] getStates() {
		return states;
	}
	
	/**
	 * Retorna el valor de las part�culas dado los estados anteriores, voltaje anterior y la corriente
	 * @param filter Object Objeto de la clase ParticleFilter, contiene par�metros y funciones asociadas al filtro de part�culas
	 * @param voltage_1 Double Voltaje anterior observado por el sensor 
	 * @param current Double Corriente observada por el sensor
	 * @return Double[][] Valor actualizado de las part�culas
	 */
	public double[][] calculationParticles(ParticleFilter filter, double voltage_1, double current)
	{
		int numberParticles		= filter.getNumberParticles();
		double[][] particles	= new double[numberStates][numberParticles];
		double[][] particles_1	= filter.getParticles();
		
		double[][] random1 		= Matlab.randn(numberParticles, 1);
		double[][] random2 		= Matlab.randn(numberParticles, 1);
		
		double w1, w2;
		
		// Estimaci�n del estado por cada part�cula
		for(int i=0; i<numberParticles; i++)
		{	
			w1 		 		= muW1 + sigmaW1*random1[i][0];						
			particles[0][i] = particles_1[0][i] + w1;
			
			w2 		 		= muW2 + sigmaW2*random2[i][0];
			particles[1][i] = Math.abs(particles_1[1][i] - voltage_1 * current * deltaT_Ecrit + w2); // En caso de tener part�culas menores que cero
		}
		
		return particles;
	}
	
	/**
	 * M�todo que realiza la estimaci�n de los estados
	 * @param filter Object Objeto de la clase ParticleFilter, contiene par�metros y funciones asociadas al filtro de part�culas
	 */
	public void StatesEstimation(ParticleFilter filter) {
		double[][] weights 		= filter.getWeights();
		double[][] particles	= filter.getParticles();
		int k					= filter.getK();
		
		for(int i = 0; i<filter.getNumberParticles(); i++){
			states[0][k] = states[0][k] + particles[0][i] * weights[i][0];
			states[1][k] = states[1][k] + particles[1][i] * weights[i][0];
		}
	}
	
	/**
	 * Establece el voltaje estimado
	 * @param voltage Double[][] Voltaje
	 */
	/*public void setVoltageEstimated(double[][] voltage){
		voltageEstimated = voltage;
	}*/
	
	/**
	 * Retorna el voltaje estimado para cada part�cula dado las part�culas y la corriente
	 * @param filter Object Objeto de la clase ParticleFilter, contiene par�metros y funciones asociadas al filtro de part�culas
	 * @param current
	 */
	public void setVoltageEstimated_k(ParticleFilter filter, double current){		
		double[][] particles = filter.getParticles();
		for(int i=0; i<filter.getNumberParticles(); i++) {
			voltageEstimated_k[0][i] = vL+(v0-vL)*Math.exp(gamma*(particles[1][i]-1))+alpha*vL*(particles[1][i]-1)+(1-alpha)*vL*(Math.exp(-beta)-Math.exp(-beta*Math.sqrt(particles[1][i])))-particles[0][i]*current+Nu;
		}
	}
	
	/**
	 * Retorna el voltaje estimado de cada part�cula en la iteraci�n k
	 * @return Double[][] Voltaje estimado de cada part�cula
	 */
	public double[][] getEstimatedVoltageK(){
		return voltageEstimated_k;
	}
	
	/**
	 * Retorna el voltaje estimado de la iteraci�n k
	 * @param k Int Iteraci�n k
	 * @return Double Voltaje estimado de la iteraci�n k
	 */
	public double getEstimatedVoltage(int k) {
		return voltageEstimated[0][k];
	}
	
	/**
	 * Retorna el voltaje estimado
	 * @return Double[][] Voltaje estimado
	 */
	public double[][] getEstimatedVoltage() {
		return voltageEstimated;
	}
	
	/**
	 * Establece el voltaje estimado en la iteraci�n k dada la corriente
	 * @param k Int Iteraci�n k
	 * @param current Double Corriente observada por el sensor 
	 */
	public void setEstimatedVoltage(int k, double current) {
		voltageEstimated[0][k] = vL+(v0-vL)*Math.exp(gamma*(states[1][k]-1))+alpha*vL*(states[1][k]-1)+(1-alpha)*vL*(Math.exp(-beta)-Math.exp(-beta*Math.sqrt(states[1][k])))-states[0][k]*current+Nu;
	}
	
	/**
	 * Retorna el n�mero de estados del modelo
	 * @return Int N�mero de estados
	 */
	public int getNumberStates(){
		return numberStates;
	}
	
	/**
	 * Retorna el n�mero de observaciones del modelo
	 * @return Int N�mero de observaciones
	 */
	public int getNumberObservations(){
		return numberObservations;
	}
	
	/**
	 * Retorna la desviaci�n est�ndar asociada al ruido del estado X1
	 * @return Double Desviaci�n est�ndar asociada al ruido del estado X1
	 */
	public double getSigmaW1() {
		return sigmaW1;
	}
	
	/**
	 * Retorna la desviaci�n est�ndar asociada al ruido del estado X2
	 * @return Double Desviaci�n est�ndar asociada al ruido del estado X2
	 */
	public double getSigmaW2(){
		return sigmaW2;
	}
	
	/**
	 * Establece la desviaci�n est�ndar asociada al ruido del estado X2 
	 * @param sigma Double Desviaci�n est�ndar asociada al ruido del estado X2
	 */
	public void setSigmaW2(double sigma){
		sigmaW2	= sigma;
	}

	/**
	 * Retorna la desviaci�n est�ndar asociada al ruido de observaci�n
	 * @return Double Desviaci�n est�ndar asociada al ruido de observaci�n
	 */
	public double getSigmaV() {
		return sigmaV;
	}

	/**
	 * Retorna el valor inicial del estado X1
	 * @return Double Valor inicial del estado X1
	 */
	public double getStates1_0(){
		return states1_0;
	}
	
	/**
	 * Retorna el valor inicial del estado X1
	 * @return Double Valor inicial del estado X1
	 */
	public double getStates2_0(){
		return states2_0;
	}

	/**
	 * Inicializa los vectores de predicci�n (predictionSoc, predictionVoltage, predictionInternalResistance, pdfEOD)
	 * @param predictionSize Int Tama�o de la predicci�n
	 */
	public void initializePrediction(int predictionSize) {
		predictionSoc 					= Matlab.zeros(1, predictionSize);		
  		predictionVoltage 				= Matlab.zeros(1, predictionSize);
  		predictionInternalResistance	= Matlab.zeros(1, predictionSize);
  		pdfEOD 							= Matlab.zeros(1, predictionSize);
	}

	/**
	 * Establece la predicci�n del voltaje
	 * @param _predictionVoltage Double[][] Predicci�n del voltaje
	 */
	public void setPredictionVoltage(double[][] _predictionVoltage) {
		predictionVoltage = _predictionVoltage;
	}
	
	/**
	 * Establece la predicci�n del voltaje para la iteraci�n j
	 * @param j Int Iteraci�n j de la predicci�n
	 * @param value Double Valor del voltaje predicho
	 */
	public void setPredictionVoltage(int j, double value) {
		predictionVoltage[0][j] = value;
	}

	/**
	 * Establece la predicci�n del SOC
	 * @param _predictionSoc Double[][] Predicci�n del SOC
	 */
	public void setPredictionSoc(double[][] _predictionSoc) {
		predictionSoc	= _predictionSoc;
	}
	
	/**
	 * Establece la predicci�n del SOC para la iteraci�n j
	 * @param j Int Iteraci�n j de la predicci�n
	 * @param value Double Valor del SOC predicho
	 */
	public void setPredictionSoc(int j, double value) {
		predictionSoc[0][j] = value;
		
	}
	
	/**
	 * Retorna la predicci�n del voltaje en la iteraci�n j
	 * @param j Int Iteraci�n j de la predicci�n
	 * @return Double Predicci�n del voltaje en la iteraci�n j
	 */
	public double getPredictionVoltage(int j) {
		return predictionVoltage[0][j];
	}

	/**
	 * Retorna la predicci�n del SOC en la iteraci�n j
	 * @param j Int Iteraci�n j de la predicci�n
	 * @return Double Predicci�n del SOC en la iteraci�n j
	 */
	public double getPredictionSoc(int j) {
		return predictionSoc[0][j];
	}
	
	/**
	 * Retorna la predicci�n del SOC
	 * @return Double[][] Predicci�n del SOC
	 */
	public double[][] getPredictionSoc() {
		return predictionSoc;
	}

	/**
	 * Retorna la funci�n densidad de probabilidad del tiempo de descarga de la bater�a
	 * @return Double[][] Funci�n densidad de probabilidad del EOD
	 */
	public double[][] getPdfEOD() {
		return pdfEOD;
	}
	
	/**
	 * Retorna el valor de la funci�n densidad de probabilidad del tiempo de descarga de la bater�a en la iteraci�n j de la predicci�n
	 * @param j Int Iteraci�n de la predicci�n
	 * @return Double Valor de la funci�n densidad de probabilidad del tiempo de descarga de la bater�a
	 */
	public double getPdfEOD(int j) {
		return pdfEOD[0][j];
	}

	/**
	 * Establece el valor de la funci�n densidad de probabilidad del tiempo de descarga de la bater�a en la iteraci�n j de la predicci�n
	 * @param j Int Iteraci�n de la predicci�n
	 * @param value Double Valor a establecer en la funci�n densidad de probabilidad del tiempo de descarga de la bater�a
	 */
	public void setPdfEOD(int j, double value) {
		pdfEOD[0][j] = value;
	}

	/**
	 * Establece el valor de la funci�n densidad de probabilidad del tiempo de descarga de la bater�a
	 * @param _pdfEod Double[][] valor de la funci�n densidad de probabilidad del tiempo de descarga de la bater�a
	 */
	public void setPdfEOD(double[][] _pdfEod) {
		pdfEOD = _pdfEod;
	}
	
	/**
	 * M�todo que utiliza el modelo para realizar la predicci�n del voltaje y del SOC
	 * @param filter Object Objeto de la clase ParticleFilter, contiene par�metros y funciones asociadas al filtro de part�culas
	 * @param predictionCurrent Double[][] Predicci�n de la corriente
	 * @param predictionSize Int Tama�o de la predicci�n
	 */
	public void PredictionModel(ParticleFilter filter, double[][] predictionCurrent, int predictionSize) {
		// En la predicci�n el valor del estado 1 (x1) es una constante
		
		int numberParticles = filter.getNumberParticles();
		
		// Vector que contiene el valor del voltaje entregado por cada particula en cada instante
		predictionVoltage_N	= Matlab.zeros(numberParticles, predictionSize);
			
		// Vector que contiene el valor del SOC por cada particula en cada instante
		predictionSoc_N		= Matlab.zeros(numberParticles, predictionSize);

		// Predicci�n de estado SOC (3 Enfoques)
		//   1.- Enfoque VE en Modelo (Pesos constantes)  
		
		double[][] _voltage 	= Matlab.zeros(1, numberParticles);
		
		double[][] aux			= filter.getParticles();
		double[][] _particles 	= new double[aux.length][aux[0].length];
		double[][] _particles2;
		
		for(int i = 0; i<aux.length; i++){
			for(int j = 0; j<aux[0].length; j++){
				_particles[i][j]	=	aux[i][j];
			}
		}
		
		for (int t=0; t<(predictionSize - 1); t++){
			// Se regulariza el SOC con fuentes de incertidumbre
			_particles2 = RegularizationPrediction(_particles, filter); // Estado 2 regularizado
			
			for(int k = 0; k<numberParticles; k++){
				_particles[1][k]	=	_particles2[0][k];
			}

		    for (int k = 0; k<numberParticles; k++){
		    	_voltage[0][k]  	= vL+(v0-vL)*Math.exp(gamma*(_particles[1][k]-1))+alpha*vL*(_particles[1][k]-1)+(1-alpha)*vL*(Math.exp(-beta)-Math.exp(-beta*Math.sqrt(_particles[1][k]))) - _particles[0][k] * predictionCurrent[0][t];
		    	// Modelo de estado SOC para el instante t (con V observado en t)
		    	_particles[1][k] 	= _particles[1][k] - _voltage[0][k]*predictionCurrent[0][t]*deltaT_Ecrit;
		    }
		 
		    for (int n = 0; n<numberParticles; n++){
		        if (t > 1 && predictionSoc_N[n][t] <= socCutOff) {
		        	predictionSoc_N[n][t+1]	= socCutOff - delta;
		            _particles[1][n]    	= socCutOff - delta;
		        } else { // Si cruza el valor socCutOff, la part�cula se queda con un valor constante
		        	predictionSoc_N[n][t+1]	= _particles[1][n];
		        }
		    }
		    
		    for(int k=0; k<predictionVoltage_N.length; k++){
		    	predictionVoltage_N[k][t+1] = vL+(v0-vL)*Math.exp(gamma*(_particles[1][k]-1))+alpha*vL*(_particles[1][k]-1)+(1-alpha)*vL*(Math.exp(-beta)-Math.exp(-beta*Math.sqrt(_particles[1][k])))-_particles[0][k]*predictionCurrent[0][t];
		    }
		    
		}

		// 2.- Enfoque propagaci�n con 2 part�culas 
	}
	
	/**
	 * M�todo que regulariza las part�culas del estado X2 del modelo
	 * @param _particles Double[][] Part�culas a ser regularizadas
	 * @param filter Object Objeto de la clase ParticleFilter, contiene par�metros y funciones asociadas al filtro de part�culas
	 * @return Retorna las part�culas del estado X2 regularizadas
	 */
	public double[][] RegularizationPrediction( double[][] _particles, ParticleFilter filter) {

		// Uso de PDF de una fuente de incertidumbre (Epanechnikov)
		double[][] state2 	= new double[1][_particles[0].length];
		for(int i = 0; i<_particles[0].length; i++){
			state2[0][i] 	= _particles[1][i];
		}
		
		int numberStates 				= state2.length;					// N� de estados es 1 (estado SOC)
		double cNumberStates   			= 1;								// Volumen de esfera unitaria (1 dimension)
		double[][] state2Regularized  	= Matlab.zeros(numberStates, filter.getNumberParticles() );
		
		float aux 			= (float) 1/(numberStates + 4);
		double A			= Math.pow((8 * Math.pow(cNumberStates, -1) * (numberStates + 4) * Math.pow(2*Math.sqrt(Math.PI), numberStates)), aux);
		double H			= A*Math.pow(filter.getNumberParticles(), -aux);
		aux 				= (float) H/10;
		H         			= aux;											// Ajuste
		double deviation 	= Matlab.std(state2, 0, "fila");
		
		int index;
		double[][] kernel 	= filter.getKernel();
		
		// Se buscan los �ndices de los elementos con valor mayor o igual al n�mero aleatorio
		for (int j = 0; j<filter.getNumberParticles(); j++){     										
		    index            			= Matlab.find_ma_i(filter.getDistributionEpanechnikov(), Math.random());	// Se crea un valor aleatorio entre 0 y 1
		    state2Regularized[0][j] 	= state2[0][j] + H*deviation*kernel[index][0];
		}
		
		return state2Regularized;

	}

	/**
	 * Retorna el valor del voltaje entregado por cada particula (N) en cada instante de predicci�n
	 * @return Double[][] El valor del voltaje entregado por cada particula en cada instante de predicci�n
	 */
	public double[][] getVoltage_N() {
		return predictionVoltage_N;
	}

	/**
	 * Retorna el valor del SOC por cada particula (N) en cada instante de predicci�n
	 * @return Double[][] El valor del SOC por cada particula en cada instante de predicci�n
	 */
	public double[][] getPredictionSoc_N() {
		return predictionSoc_N;
	}

	/**
	 * Retorna la predicci�n del voltaje
	 * @return Double[][] Predicci�n del voltaje
	 */
	public double[][] getPredictionVoltage() {
		return predictionVoltage;
	}

	/**
	 * Establece la predicci�n del SOC en la iteraci�n k
	 * @param value Double Valor a ser establecido
	 * @param k Int Iteraci�n k de la predicci�n
	 */
	public void setPredictionSoc(double value, int k) {
		predictionSoc[0][k] = value;
		
	}

	/**
	 * Restablece los valores de la predicci�n del voltaje
	 * @param predictionSize Int Tama�o de la predicci�n
	 */
	public void clearPredictionVoltage(int predictionSize) {
		predictionVoltage = Matlab.zeros(1, predictionSize);
		
	}

	/**
	 * Restablece los valores de la predicci�n del SOC
	 * @param predictionSize Int Tama�o de la predicci�n
	 */
	public void clearPredictionSOC(int predictionSize) {
		predictionSoc = Matlab.zeros(1, predictionSize);
	}

	/**
	 * Establece el valor de la predicci�n de la resistencia interna del modelo en la iteraci�n j
	 * @param i Int Iteraci�n j de la predicci�n
	 * @param value Double Valor a ser establecido
	 */
	public void setPredictionInternalResistance(int j, double value) {
		predictionInternalResistance[0][j] = value;
	}

}