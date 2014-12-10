package soc;

import mathlib.Matlab;

/**
 * 
 * @author Alex J. D�az Mill�n
 *
 */
public class ParticleFilter {
	
	private int numberParticles 	= 40;				// N�mero de Part�culas
	private String strategy 		= "multinomial";	// Estrategia de remuestreo
	private double[][] weights, particles;
	private int k;
	
	private double[][] kernel;							// N�cleos definidos como (x,raiz(1-x^2)) con norma de x menor a 1
	private double[][] distributionEpanechnikov;		// Suma acumulada de los valores del Kernel
	
	private double reliability;							// Probabilidad considerada para formar el intervalo de confianza
	private double min;									// M�nimo del Intervalo de Confianza (medido en [seg])
	private double max;									// M�ximo del Intervalo de Confianza (medido en [seg])
	private double tof;									// Instante de tiempo del EoD predicho (medido en [seg])
	
	
	/**
	 * Inicializaci�n del Algoritmo de Filtro de Part�culas
	 */
	public void initializeParticleFilter() {
		k = 1;
	}
	
	/**
	 * Inicializaci�n de las part�culas iniciales
	 * @param model Object Objeto de la clase Model, contiene los par�metros y funciones asociadas al modelo a utilizar
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
	 * Inicializaci�n de los pesos de las part�culas
	 */
	public void initializeWeights() {		
		weights = new double[numberParticles][1];
		UniformWeights();
	}
	
	/**
	 * Retorna los pesos de las part�culas
	 * @return Double[][] Pesos de las part�culas
	 */
	public double[][] getWeights(){
		return weights;
	}
	
	/**
	 * Retorna el peso de la part�cula i
	 * @param i Int Posici�n de la part�cula
	 * @return Double Peso de la part�cula
	 */
	public double getWeights(int i) {
		return weights[i][0];
	}
	
	/**
	 * Retorna los pesos actualizados dado el voltaje observado y voltaje estimado
	 * @param voltageObserved Double Voltaje observado
	 * @param model Object Objeto de la clase Model, contiene los par�metros y funciones asociadas al modelo a utilizar
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
	 * Retorna el n�mero de part�culas
	 * @return Int N�mero de part�culas
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
	 * Retorna la iteraci�n k del algoritmo
	 * @return Iteraci�n del algoritmo
	 */
	public int getK() {	
		return k;
	}
	
	/**
	 * M�todo que realiza el remuestreo de las part�culas
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
			particles = Matlab.extract(particles,idx);	// Extrae nuevas particulas (reemplazadas por degeneraci�n)
			UniformWeights();
		}
		
	}
	
	/**
	 * Establece el peso de las part�culas en forma uniforme y normalizado
	 */
	public void UniformWeights() {
		float div 			= 1/(float) numberParticles;
		weights 			= Matlab.transpose(Matlab.repeatArray(div, 1, numberParticles));
	}

	/**
	 * Retorna los pesos de las part�culas normalizados
	 * @return Double[][] Pesos de las part�culas normalizados
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
	 * Establece el valor de la iteraci�n del algoritmo
	 * @param _k Int Iteraci�n del algoritmo
	 */
	public void setK(int _k) {
		k = _k;
	}

	/**
	 * Retorna las part�culas
	 * @return Double[][] Part�culas
	 */
	public double[][] getParticles() {
		return particles;
	}
	
	/*public double[][] getParticles(int i){
		return Matlab.cut(particles, i-1, "fila");
	}*/

	/**
	 * Establece el valor de los pesos de las part�culas
	 * @param _weights Double[][] Pesos de las part�culas
	 */
	public void setWeights(double[][] _weights) {
		weights = _weights;
	}

	/**
	 * Establece el valor de las part�culas
	 * @param _particles Double[][] Part�culas
	 */
	public void setParticles(double[][] _particles) {
		particles = _particles;
	}

	/**
	 * M�todo que calcula la distribuci�n de Epanechnikov
	 */
	public void CalculationEpanechnikov() {
		double[][] xRange			= Matlab.createRange(-1,0.01,1);
		double[][] yRange			= new double[xRange.length][1];
		for(int k = 0; k<xRange.length; k++){
			yRange[k][0]			= Math.sqrt(1-Math.pow(xRange[k][0], 2));
		}
		kernel						= Matlab.createKernel(xRange, yRange);		// N�cleos definidos como (x,raiz(1-x^2)) con norma de x menor a 1
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
	 * Retorna la distribuci�n de Epanechnikov
	 * @return Double[][] Distribuci�n de Epanechnikov
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
	 * Establece el valor m�nimo del intervalo de confianza [seg]
	 * @param _min M�nimo del intervalo de confianza
	 */
	public void setMin(double _min) {
		min = _min;
	}
	
	/**
	 * Retorna el valor del m�nimo del intervalo de confianza [seg]
	 * @return Double M�nimo del intervalo de confianza
	 */
	public double getMin() {
		return min;
	}
	
	/**
	 * Establece el valor m�ximo del intervalo de confianza [seg]
	 * @param _max M�ximo del intervalo de confianza
	 */
	public void setMax(double _max) {
		max = _max;
	}
	
	/**
	 *  Retorna el valor del m�ximo del intervalo de confianza [seg]
	 * @return M�ximo del intervalo de confianza
	 */
	public double getMax() {
		return max;
	}
	
	/**
	 * Establece el valor del tiempo de descarga de la bater�a
	 * @param _tof Double Tiempo de descarga de la bater�a
	 */
	public void setTof(double _tof) {
		tof = _tof;
	}
	
	/**
	 * Retorna el valor del tiempo de descarga de la bater�a
	 * @return Double Tiempo de descarga de la bater�a
	 */
	public double getTof() {
		return tof;
	}
	
}