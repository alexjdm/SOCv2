package data;

public class ParameterizationBattery {
	
	private double[][] voltageOcSen, socSen;
	
	public ParameterizationBattery() throws Exception{
		double[][] parameterizationBattery 	= DataManagement.main("Parametrizacion.xls");
		voltageOcSen 						= new double [1][parameterizationBattery.length];
		socSen 								= new double [1][parameterizationBattery.length];
		for (int i = 0; i<parameterizationBattery.length; i++){
			for(int j = 0; j<2; j++){
				if(j==0){
					voltageOcSen[0][i] 	= parameterizationBattery[i][j];
				}
				else if(j==1){
					socSen[0][i] 		= parameterizationBattery[i][j];
				}
			}
		}
	}
	
	public double[][] getVoltageOcSen(){
		return voltageOcSen;
	}
	
	public double[][] getSocSen(){
		return socSen;
	}

}
