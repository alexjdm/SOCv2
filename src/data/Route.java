package data;

public class Route {
	
	private double[][] currentSensor, voltageSensor;
	
	public Route() throws Exception{
		double[][] route	= DataManagement.main("ruta5.xls");
		currentSensor 		= new double [route.length][1];
		voltageSensor 		= new double [route.length][1];
		for (int i = 0; i<route.length; i++){
			for(int j = 0; j<2; j++){
				if(j==0){
					currentSensor[i][0] = route[i][j];
				}
				else if(j==1){
					voltageSensor[i][0] = route[i][j];
				}
			}
		}
	}
	
	public double[][] getCurrentSensor(){
		return currentSensor;
	}
	
	public double[][] getVoltageSensor(){
		return voltageSensor;
	}

}