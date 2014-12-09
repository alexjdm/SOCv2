package soc;

import mathlib.Matlab;

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
	
	
	// Inicialización del Algoritmo de Filtro de Partículas
	public void initializeParticleFilter() {
		k = 1;
	}
	
	// Inicialización de las partículas iniciales
	public void initializeParticles(Model model) {
		particles			= new double[model.getNumberStates()][numberParticles];
		double[][] random1 	= Matlab.randn(1,numberParticles);
		double[][] random2 	= Matlab.rand(1,numberParticles);
		for(int i = 0; i<numberParticles; i++){
			particles[0][i] = model.getStates1_0() + model.getSigmaW1() * random1[0][i];
			particles[1][i] = 0.8 + 0.1 * model.getStates2_0() * random2[0][i]; 
		}
	}
	
	// Inicialización de los pesos de las partículas
	public void initializeWeights() {		
		weights = new double[numberParticles][1];
		UniformWeights();	// Genera todas las partículas de peso uniforme y normalizado para Resampling.m
	}
	
	// Retorna los pesos
	public double[][] getWeights(){
		return weights;
	}
	
	public double getWeights(int i) {
		return weights[i][0];
	}
	
	// Retorna los pesos actualizados dado el voltaje observado y voltaje estimado
	public void updateWeights(double voltageObserved, Model model)
	{
		double[][] voltageEstimated	= model.getVoltageEstimated_k();	
		double sigmaV				= model.getSigmaV();
		double error, pdfNu;
		
		for(int i=0; i<numberParticles; i++)
		{
			error			= voltageObserved - voltageEstimated[0][i];
			pdfNu  			= Matlab.normpdf(error, 0, sigmaV);
			weights[i][0]	= weights[i][0] * pdfNu;
		}
	}
	
	public int getNumberParticles(){
		return numberParticles;
	}

	public String getStrategy() {
		return strategy;
	}
	
	public int getK() {	
		return k;
	}
	
	public void Resampling() {
		
		int[][] idx = new int[1][numberParticles];
		
		switch (strategy){
			case "multinomial":
				// Vector columna de particulas (muestras) con reemplazo y peso nuevo
				idx = Matlab.randsample(1, numberParticles, numberParticles, weights);	
				break;
			case "sistematico":
				double[][] edges 		 	= Matlab.cumsum(weights, 0, "columna");
				edges 					  	= Matlab.addValue(edges, 1, 0.0);
				edges[edges.length-1][0] 	= 1;
				float u1 					= (float) (Math.random()/numberParticles);
				float u2					= 1/((float) numberParticles);
				idx 						= Matlab.histc(Matlab.sample(u1, u2, 1), edges);
				break;
		   default:
			   System.out.print("Estrategia no implementada o desconocida");
			   break;
		}
		
		particles = Matlab.extract(particles,idx);	// Extrae nuevas particulas (reemplazadas por degeneración)
		UniformWeights();
	}
	
	public void UniformWeights() {
		float div 				= 1/(float) numberParticles;
		weights 				= Matlab.transpose(Matlab.repeatArray(div, 1, numberParticles));
	}

	public double[][] WeightsNormalization() {
		double sumaWk 	  	= Matlab.sum(weights, 0, "columna");
		float aux3 		  	= 0;
		for(int i = 0; i<weights.length; i++){
			aux3 	 		= (float) (weights[i][0]/sumaWk);
			weights[i][0] 	= aux3;
		}
		return weights;
	}

	public void setK(int k2) {
		k = k2;
	}

	public double[][] getParticles() {
		return particles;
	}
	
	public double[][] getParticles(int i){
		return Matlab.cut(particles, i-1, "fila");
	}

	
	public void setWeights(double[][] _weights) {
		weights = _weights;
	}

	public void setParticles(double[][] _particles) {
		particles = _particles;
	}

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


	public double[][] getKernel() {
		return kernel;
	}

	public double[][] getDistributionEpanechnikov() {
		return distributionEpanechnikov;
	}

	public void setReliability(double _reliability) {
		reliability = _reliability; 
		
	}

	public double getReliability() {
		return reliability;
	}

	public void setMin(double _min) {
		min = _min;
	}
	
	public double getMin() {
		return min;
	}
	
	public void setMax(double _max) {
		max = _max;
	}
	
	public double getMax() {
		return max;
	}
	
	public void setTof(double _tof) {
		tof = _tof;
	}
	
	public double getTof() {
		return tof;
	}
	
}