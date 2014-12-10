package soc;

import mathlib.Matlab;

/**
 * 
 * @author Alex J. Díaz Millán
 *
 */
public class ParticleFilter {
	
	private int numberParticles 	= 40;				// Número de Partículas
	private String strategy 		= "multinomial";	// Estrategia de remuestreo
	private double[][] weights, particles;
	private int k;
	
	private double[][] kernel;							// Núcleos definidos como (x,raiz(1-x^2)) con norma de x menor a 1
	private double[][] distributionEpanechnikov;		// Suma acumulada de los valores del Kernel
	
	private double reliability;							// Probabilidad considerada para formar el intervalo de confianza
	private double min;									// Mínimo del Intervalo de Confianza (medido en [seg])
	private double max;									// Máximo del Intervalo de Confianza (medido en [seg])
	private double tof;									// Instante de tiempo del EoD predicho (medido en [seg])
	
	
	/**
	 * Inicialización del Algoritmo de Filtro de Partículas
	 */
	public void initializeParticleFilter() {
		k = 1;
	}
	
	/**
	 * Inicialización de las partículas iniciales
	 * @param model Object Objeto de la clase Model, contiene los parámetros y funciones asociadas al modelo a utilizar
	 */
	public void initializeParticles(Model model) {
		particles			= new double[model.getNumberStates()][numberParticles];
		double[][] random1 	= Matlab.randn(1,numberParticles);
		double[][] random2 	= Matlab.rand(1,numberParticles);
		for(int i = 0; i<numberParticles; i++){
			particles[0][i] = model.getStates1_0() + model.getSigmaW1() * random1[0][i];
			particles[1][i] = 0.8 + 0.1 * model.getStates2_0() * random2[0][i]; 
		}
	}
	
	/**
	 * Inicialización de los pesos de las partículas
	 */
	public void initializeWeights() {		
		weights = new double[numberParticles][1];
		UniformWeights();
	}
	
	/**
	 * Retorna los pesos de las partículas
	 * @return Double[][] Pesos de las partículas
	 */
	public double[][] getWeights(){
		return weights;
	}
	
	/**
	 * Retorna el peso de la partícula i
	 * @param i Int Posición de la partícula
	 * @return Double Peso de la partícula
	 */
	public double getWeights(int i) {
		return weights[i][0];
	}
	
	/**
	 * Retorna los pesos actualizados dado el voltaje observado y voltaje estimado
	 * @param voltageObserved Double Voltaje observado
	 * @param model Object Objeto de la clase Model, contiene los parámetros y funciones asociadas al modelo a utilizar
	 */
	public void updateWeights(double voltageObserved, Model model)
	{
		double[][] voltageEstimated	= model.getEstimatedVoltageK();
		double sigmaV				= model.getSigmaV();
		double error, pdfNu;
		
		for(int i=0; i<numberParticles; i++)
		{
			error			= voltageObserved - voltageEstimated[0][i];
			pdfNu  			= Matlab.normpdf(error, 0, sigmaV);
			weights[i][0]	= weights[i][0] * pdfNu;
		}
	}
	
	/**
	 * Retorna el número de partículas
	 * @return Int Número de partículas
	 */
	public int getNumberParticles(){
		return numberParticles;
	}

	/**
	 * Retorna la estrategia de remuestreo a utilizar 
	 * @return String Estrategia de remuestreo
	 */
	public String getStrategy() {
		return strategy;
	}
	
	/**
	 * Retorna la iteración k del algoritmo
	 * @return Iteración del algoritmo
	 */
	public int getK() {	
		return k;
	}
	
	/**
	 * Método que realiza el remuestreo de las partículas
	 */
	public void Resampling() {
		
		int[][] idx = new int[1][numberParticles];
		boolean ok	= false;
		
		switch (strategy){
			case "multinomial":
				// Vector columna de particulas (muestras) con reemplazo y peso nuevo
				idx = Matlab.randsample(1, numberParticles, numberParticles, weights);
				ok	= true;
				break;
			case "sistematico":
				double[][] edges 		 	= Matlab.cumsum(weights, 0, "columna");
				edges 					  	= Matlab.addValue(edges, 1, 0.0);
				edges[edges.length-1][0] 	= 1;
				float u1 					= (float) (Math.random()/numberParticles);
				float u2					= 1/((float) numberParticles);
				idx 						= Matlab.histc(Matlab.sample(u1, u2, 1), edges);
				ok	= true;
				break;
		   default:
			   System.out.print("Error: Estrategia desconocida");
			   ok	= false;
			   break;
		}
		
		if(ok){
			particles = Matlab.extract(particles,idx);	// Extrae nuevas particulas (reemplazadas por degeneración)
			UniformWeights();
		}
		
	}
	
	/**
	 * Establece el peso de las partículas en forma uniforme y normalizado
	 */
	public void UniformWeights() {
		float div 			= 1/(float) numberParticles;
		weights 			= Matlab.transpose(Matlab.repeatArray(div, 1, numberParticles));
	}

	/**
	 * Retorna los pesos de las partículas normalizados
	 * @return Double[][] Pesos de las partículas normalizados
	 */
	public double[][] WeightsNormalization() {
		double sumaWk 	  	= Matlab.sum(weights, 0, "columna");
		float aux3 		  	= 0;
		for(int i = 0; i<weights.length; i++){
			aux3 	 		= (float) (weights[i][0]/sumaWk);
			weights[i][0] 	= aux3;
		}
		return weights;
	}

	/**
	 * Establece el valor de la iteración del algoritmo
	 * @param _k Int Iteración del algoritmo
	 */
	public void setK(int _k) {
		k = _k;
	}

	/**
	 * Retorna las partículas
	 * @return Double[][] Partículas
	 */
	public double[][] getParticles() {
		return particles;
	}
	
	/*public double[][] getParticles(int i){
		return Matlab.cut(particles, i-1, "fila");
	}*/

	/**
	 * Establece el valor de los pesos de las partículas
	 * @param _weights Double[][] Pesos de las partículas
	 */
	public void setWeights(double[][] _weights) {
		weights = _weights;
	}

	/**
	 * Establece el valor de las partículas
	 * @param _particles Double[][] Partículas
	 */
	public void setParticles(double[][] _particles) {
		particles = _particles;
	}

	/**
	 * Método que calcula la distribución de Epanechnikov
	 */
	public void CalculationEpanechnikov() {
		double[][] xRange			= Matlab.createRange(-1,0.01,1);
		double[][] yRange			= new double[xRange.length][1];
		for(int k = 0; k<xRange.length; k++){
			yRange[k][0]			= Math.sqrt(1-Math.pow(xRange[k][0], 2));
		}
		kernel						= Matlab.createKernel(xRange, yRange);		// Núcleos definidos como (x,raiz(1-x^2)) con norma de x menor a 1
		distributionEpanechnikov	= Matlab.cumsum(yRange, 0, "columna");		// Suma acumulada de los valores del Kernel
		
		// Se normaliza la distribucion de Epanechnikov
		for(int k = 0; k<distributionEpanechnikov.length; k++){
			distributionEpanechnikov[k][0]		= (float) (distributionEpanechnikov[k][0]/distributionEpanechnikov[distributionEpanechnikov.length-1][0]);	
		}
	}

	/**
	 * Retorna el kernel generado en CalculationEpanechnikov()
	 * @return Double[][] Kernel
	 */
	public double[][] getKernel() {
		return kernel;
	}

	/**
	 * Retorna la distribución de Epanechnikov
	 * @return Double[][] Distribución de Epanechnikov
	 */
	public double[][] getDistributionEpanechnikov() {
		return distributionEpanechnikov;
	}

	/**
	 * Establece el nivel de confianza
	 * @param _reliability Double Confianza
	 */
	public void setReliability(double _reliability) {
		reliability = _reliability; 
	}
	
	/**
	 * Retorna el nivel de confianza
	 * @return Double Confianza
	 */
	public double getReliability() {
		return reliability;
	}

	/**
	 * Establece el valor mínimo del intervalo de confianza [seg]
	 * @param _min Mínimo del intervalo de confianza
	 */
	public void setMin(double _min) {
		min = _min;
	}
	
	/**
	 * Retorna el valor del mínimo del intervalo de confianza [seg]
	 * @return Double Mínimo del intervalo de confianza
	 */
	public double getMin() {
		return min;
	}
	
	/**
	 * Establece el valor máximo del intervalo de confianza [seg]
	 * @param _max Máximo del intervalo de confianza
	 */
	public void setMax(double _max) {
		max = _max;
	}
	
	/**
	 *  Retorna el valor del máximo del intervalo de confianza [seg]
	 * @return Máximo del intervalo de confianza
	 */
	public double getMax() {
		return max;
	}
	
	/**
	 * Establece el valor del tiempo de descarga de la batería
	 * @param _tof Double Tiempo de descarga de la batería
	 */
	public void setTof(double _tof) {
		tof = _tof;
	}
	
	/**
	 * Retorna el valor del tiempo de descarga de la batería
	 * @return Double Tiempo de descarga de la batería
	 */
	public double getTof() {
		return tof;
	}
	
}