package data;

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException; 
import java.util.ArrayList; 
import java.util.Iterator; 
import java.util.List; 

import org.apache.poi.hssf.usermodel.HSSFCell; 
import org.apache.poi.hssf.usermodel.HSSFRow; 
import org.apache.poi.hssf.usermodel.HSSFSheet;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.ss.usermodel.Cell;

/**
 * Clase encargada del manejo de datos en hojas de cálculo Excel
 * @author Alexhp
 *
 */
public class DataManagement {
	
	public static double[][] main(String nombre) throws Exception {
		String nombre_archivo = nombre;
		List ListaDatos = new ArrayList();
		ListaDatos = leerDatosExel(nombre_archivo);
		double[][] Datos = listtodouble(ListaDatos);
		return Datos;
	 }
	
	/**
	 * Función que guarda los datos ingresados a una hoja de cálculo Excel
	 * @param matriz1 double[][] matriz de entrada
	 * @param nombre String nombre del archivo en el que son guardados los datos
	 */
	public static void grabarDatosExel(double[][] matriz1, String nombre) {
		String nombre_archivo = nombre;
		HSSFWorkbook libro = new HSSFWorkbook();
		HSSFSheet hoja1 = libro.createSheet("Hoja1");
		//System.out.println(matriz1.length);
		//System.out.println(matriz1[0].length);
		for(int i = 0; i<matriz1.length; i++){
			//for(int j = 0; j<matriz[0].length; j++){
				HSSFRow fila = hoja1.createRow(i);	
				HSSFCell celda = fila.createCell(0);
				celda.setCellValue(matriz1[i][0]);	
			//}
		}
		try{
			FileOutputStream arch = new FileOutputStream(nombre_archivo);
			libro.write(arch);
			arch.close();
		}catch(Exception e){
			System.out.println("No sirve");
		}
	}
	
	private static List leerDatosExel(String n_sheet) {
		String filename = n_sheet;
		List sheetData = new ArrayList();
		FileInputStream fis = null;	 
		try {
			fis = new FileInputStream(filename);
		 	HSSFWorkbook workbook = new HSSFWorkbook(fis);
		 	HSSFSheet sheet = workbook.getSheetAt(0);		 	
		 	Iterator rows = sheet.rowIterator();
		 	while (rows.hasNext()) {
		 		HSSFRow row = (HSSFRow) rows.next();
		 		Iterator cells = row.cellIterator();
		 		List data = new ArrayList();
		 		while (cells.hasNext()) {
		 			HSSFCell cell = (HSSFCell) cells.next();
		 			data.add(cell);
		 		}
		 		sheetData.add(data);
		 	}
		 } catch (IOException e) {
			 e.printStackTrace();
		 }
		return sheetData;
	}
	
	private static double[][] listtodouble(List ListaDatos){
		double[][] Datos = new double [ListaDatos.size()][2];
		for (int i = 0; i < ListaDatos.size(); i++) {
			List list = (List) ListaDatos.get(i);
			for (int j = 0; j < 2; j++) { //j < list.size()
				Cell cell = (Cell) list.get(j);
				if (cell.getCellType() == Cell.CELL_TYPE_NUMERIC) {
					Datos [i][j]= (double) cell.getNumericCellValue();
					//System.out.print(Datos [i][j]);
				} 
				if (j < list.size() - 1) {
					//System.out.print(", ");
				}
			}
			//System.out.println("");
		 }
		return Datos;
	}
	
	static void visualizarDatos(double[][] Datos){
		for (int i = 0; i < Datos.length; i++) { 	//filas
			for (int j = 0; j < Datos[0].length; j++) {			//j < Datos[0].length columnas
				System.out.print(Datos [i][j]); 
				if (j < 2 - 1) {
					System.out.print(", ");
				}
			}
			System.out.println("");
		 }
	}
	
}